""" This script is to calculate the varying pupil of HET """
from shapely.geometry import Polygon, Point
from shapely.affinity import translate, rotate
from shapely.ops import transform

import numpy as np
import matplotlib.pyplot as plt

from . import HETparams

class PrimaryMirror(object):
    """ This object is HET's Primary mirror """
    def __init__(self,BadMirrorSegmentCoords=None):
        """ Create the primary mirror of HET.
        Input:
              BadMirrorSegmentCoords: (optional) See HETparams.py for default value and format of how to list bad primary mirror segments.
                                     Example: [(+3,7),(+1,7)]
        """
        # all the coordinates below are on the curved specrical surface of the primary mirrior.
        # so it is a 2D manifold embeded in the 3d world.

        self.SegmentMirrorApothem = HETparams.SegmentMirrorApothem

        # Gaps between mirros for creating the Primary mirror Honey comb grid
        self.Gap = HETparams.Gap
        self.ColumnGap = self.Gap*np.sqrt(3)/2.


        self.SegmentMirror = self.create_single_segment_mirror(self.SegmentMirrorApothem)

        self.ColumnPositions = (3*self.SegmentMirrorApothem/np.sqrt(3)+self.ColumnGap) * np.arange(6) # Only positive coordinates

        # Create the Primary Mirror Grid
        self.PrimaryMirrorGrid = self.create_primary_mirror_grid(ColumnPositions = self.ColumnPositions)

        # Create the primary mirro by shifting the Segment mirrors on to the grid, and the dictionary to number the segments
        self.PrimaryMirrorSegmentsList, self.SegMirrorIndexDic = self.create_primarymirror_segments_and_gridindx(self.PrimaryMirrorGrid)
        
        # Optinally create a list of Bad Mirror segments which has dropped out of the allignemnt
        if BadMirrorSegmentCoords is None:
            self.BadMirrorSegmentCoords = HETparams.BadMirrorSegmentCoords
        else:
            self.BadMirrorSegmentCoords = BadMirrorSegmentCoords

        #Union the segments to create a full single Primary mirror of HET
        self.PrimaryMirror = self.generate_combined_primary_mirror(self.PrimaryMirrorSegmentsList,self.BadMirrorSegmentCoords)

    def create_single_segment_mirror(self,apothem=None):
        """ Returns one Hexagonal segment mirror of primary centered at (0,0) """
        apothem = apothem or self.SegmentMirrorApothem
        # Hexagon Segment mirror in primary (Hex Polygon centered on (0,0))
        SegmentMirror = Polygon([(apothem/np.sqrt(3),apothem),
                                 (2*apothem/np.sqrt(3),0),
                                 (apothem/np.sqrt(3),-apothem),
                                 (-apothem/np.sqrt(3),-apothem),
                                 (-2*apothem/np.sqrt(3),0),
                                 (-apothem/np.sqrt(3),apothem)])

        print('Area of one Hex Segment mirror of primary = {0} m^2'.format(SegmentMirror.area))
        return SegmentMirror

    def create_primary_mirror_grid(self,ColumnPositions=None, Gap=None):
        """ Returns the Primary Grid of the Primary mirror, on which the segment mirros are to be palced """
        if ColumnPositions is None:
            ColumnPositions =  self.ColumnPositions
        if Gap is None:
            Gap = self.Gap

        PrimaryMirrorGrid = {}
        for j,cposition in enumerate(ColumnPositions):
            if j == 0:
                # Central column 
                PrimaryMirrorGrid[cposition] = np.arange(-5,6)*(2*self.SegmentMirrorApothem + Gap)
            elif j == 1:
                # Next odd column
                PrimaryMirrorGrid[cposition] = (self.SegmentMirrorApothem + Gap/2.) + np.arange(-5,5)*(2*self.SegmentMirrorApothem + Gap)      
            elif j%2 == 0:
                # Next even column
                PrimaryMirrorGrid[cposition] = PrimaryMirrorGrid[ColumnPositions[j-2]][1:-1]
            elif j%2 == 1:
                # Nex odd column
                PrimaryMirrorGrid[cposition] = PrimaryMirrorGrid[ColumnPositions[j-2]][1:-1]

            # Add the symetrrically opposite -ve column position
            PrimaryMirrorGrid[-cposition] = PrimaryMirrorGrid[cposition]

        return PrimaryMirrorGrid

    def create_primarymirror_segments_and_gridindx(self, PrimaryMirrorGrid = None):
        """ Creates the primary mirrors by shifting the Segment mirrors on to the PrimaryMirrorGrid.
        Returns: the list of segment mirros as well as a dictionary to index the mirror."""
        if PrimaryMirrorGrid is None:
            PrimaryMirrorGrid = self.PrimaryMirrorGrid

        PrimaryMirrorSegmentsList = []
        SegMirrorIndexDic = {} # dictionary to number the segments
        i = 0
        for cposition in PrimaryMirrorGrid:
            for gridpos in PrimaryMirrorGrid[cposition]:
                PrimaryMirrorSegmentsList.append( translate(self.SegmentMirror,xoff=cposition,yoff=gridpos) )
                SegMirrorIndexDic[cposition,gridpos] = i
                i += 1

        return PrimaryMirrorSegmentsList, SegMirrorIndexDic

    def generate_combined_primary_mirror(self,PrimaryMirrorSegmentsList=None,BadMirrorSegmentCoords=None):
        """ Creates the combined Primary mirror by union of segments in 
        PrimaryMirrorSegmentsList and remove bad segments in BadMirrorSegmentCoords.
        Returns: the combined Primary mirror of HET.
        """
        if PrimaryMirrorSegmentsList is None:
            PrimaryMirrorSegmentsList = self.PrimaryMirrorSegmentsList 
        if BadMirrorSegmentCoords is None :
            BadMirrorSegmentCoords = self.BadMirrorSegmentCoords

        # Take union of all segments to create Primary Mirror
        print('No of Segments in Primary Mirror = {0}'.format(len(PrimaryMirrorSegmentsList)))

        PrimaryMirror  = PrimaryMirrorSegmentsList[0].union(PrimaryMirrorSegmentsList[1])
        for i in range(2,len(PrimaryMirrorSegmentsList)):
            PrimaryMirror = PrimaryMirror.union(PrimaryMirrorSegmentsList[i])

        # Remove Bad mirror segments (if any) from the grid

        for cindx,rindx in BadMirrorSegmentCoords:
            try:
                cposition = np.sign(cindx)*self.ColumnPositions[np.abs(cindx)]
                gridpos = self.PrimaryMirrorGrid[cposition][rindx]
            except (KeyError,IndexError) as e:
                print('Error: In coordinates of Bad Segment mirror {0}'.format((cindx,rindx)))
                print('Give correct coordintes..')
                continue
            BadMirrorIndx = self.SegMirrorIndexDic[cposition,gridpos]
            BadMirror = PrimaryMirrorSegmentsList[BadMirrorIndx]
            # Remove this Bad mirror from Primary mirror
            PrimaryMirror = PrimaryMirror.difference(BadMirror)


        print('Area of Primary Mirror = {0} m^2'.format(PrimaryMirror.area))
        return PrimaryMirror


class Pupil(object):
    """ Illumintaed Active Pupil of HET """
    def __init__(self):
        self.PupilDia = HETparams.PupilDia
        self.CentralShadowDiaWFC = HETparams.CentralShadowDiaWFC
        self.BeamPosition = HETparams.BeamPosition
        self.BeamWidth = HETparams.BeamWidth

        self.FullPupilShape = self.create_illuminated_pupil()

    def create_illuminated_pupil(self,PupilDia=None,CentralShadowDiaWFC =None,BeamPosition= None, BeamWidth = None):
        """ Returns the full cirvular illuminated pupil with all the secondary obscurations """
        PupilDia = PupilDia or self.PupilDia 
        CentralShadowDiaWFC = CentralShadowDiaWFC or self.CentralShadowDiaWFC
        BeamPosition = BeamPosition or self.BeamPosition
        BeamWidth = BeamWidth or self.BeamWidth

        # Create the Pupil
        FullPupil = Point(0,0).buffer(PupilDia/2.0)
        # Central obsceuration due to WFC
        WFCshadow = Point(0,0).buffer(CentralShadowDiaWFC/2.0)
        PupilShape = FullPupil.difference(WFCshadow)  # Add shadow of the tracker

        # Shadow of the two Tracker beams 
        # Approximate sice from figure 7 in Wolf et. al. on MARS
        BeamLeft = Polygon([(-BeamPosition,10),(-BeamPosition,-10),(-BeamPosition-BeamWidth,-10),(-BeamPosition-BeamWidth,10)])
        BeamRight = Polygon([(BeamPosition,10),(BeamPosition,-10),(BeamPosition+BeamWidth,-10),(BeamPosition+BeamWidth,10)])
        SupportBeams = BeamLeft.union(BeamRight)

        #Update Effective Pupil shadow
        PupilShape = PupilShape.difference(SupportBeams)  # Add shadow of the tracker
        return PupilShape

    def projected_mirror_shape(self, Mirror, xcenter, ycenter, radius_curv = HETparams.RadiusOfCurvatureHETprimary):
        """ Returns the projected shape of the Mirror sitting on the spehere of radius of curvature radius_curv.
        The projection plane is the plane prependicular to line connecteding (xcenter,ycenter) and center of curvature of mirror
        """
        def transform_func(x,y):
            """ Transformation function for projection """
            r = np.sqrt((x-xcenter)**2 + (y-ycenter)**2) # radial distance from center od the projection
            r_t =  radius_curv * np.sin(r/radius_curv)  # radial distance in projected plane.
            
            x_t = r_t * x/r   # new projected x and y distance in plane
            y_t = r_t * y/r   

            return (x_t, y_t)

        # Apply the transformation function to the Mirror shape and return it.
        return transform(transform_func, Mirror)

    def EffectiveActivePupil(self,PrimaryMirror, xoff=0,yoff=0, PupilShape=None, doMirrorCurvProjection = True):
        """ 
        Input: 
             xoff : X position offset of pupil on the Primarry mirror (unit in meters)
             yoff : Y position offset of pupil on the Primarry mirror (unit in meters)
             Note - xoff and yoff are not same as tracker position. There is a multiplicative factor (~ 2) since it is half way between center of curvature and mirror.

             PupilShape: The pupil projected on to primary mirror. (by default uses self.FullPupilShape)
             doMirrorCurvProjection = True , project the primary mirror to the prependicular plane seen by the star for calculating effective area

        Returns the effective Active Pupil shape for a given pupil position by shifting and intersecting with Primary mirror"""
        if PupilShape is None:
            PupilShape = self.FullPupilShape

        ShiftedPupil = translate(PupilShape,xoff=xoff,yoff=yoff)
        
        if doMirrorCurvProjection: # Apply the mirror curvature projection
            PrimaryMirror = self.projected_mirror_shape(PrimaryMirror, xoff, yoff)

        EffectivePupil = ShiftedPupil.intersection(PrimaryMirror)
        return EffectivePupil

###################################
if __name__ == "__main__":
    xshifts = np.linspace(-3.9/2,3.9/2,1000)
    AreaDrift = [EffectivePupil(xoff=x).area for x in xshifts]
    plt.plot(xshifts,AreaDrift)
    plt.show()
    

