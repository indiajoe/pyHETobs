""" This script is to calculate the varying pupil of HET """
from shapely.geometry import Polygon, Point
from shapely.affinity import translate, rotate
import numpy as np
import matplotlib.pyplot as plt

SegmentMirrorApothem = 0.5  # meters

# Hexagon Segment mirror in primary (Hex Polygon centered on (0,0))
apothem = SegmentMirrorApothem
SegmentMirror = Polygon([(apothem/np.sqrt(3),apothem),
                         (2*apothem/np.sqrt(3),0),
                         (apothem/np.sqrt(3),-apothem),
                         (-apothem/np.sqrt(3),-apothem),
                         (-2*apothem/np.sqrt(3),0),
                         (-apothem/np.sqrt(3),apothem)])

print('Area of one Hex Segment mirror of primary = {0} m^2'.format(SegmentMirror.area))

# Create the Primary mirror Honey comb grid
Gap = 0.011  # 6.2 to 15.8 mm # 11 mm on average
ColumnGap = Gap*np.sqrt(3)/2.

ColumnPositions = (3*SegmentMirrorApothem/np.sqrt(3)+ColumnGap) * np.arange(6) # Only positive coordinates
PrimaryMirrorGrid = {}
for j,cposition in enumerate(ColumnPositions):
    if j == 0:
        # Central column 
        PrimaryMirrorGrid[cposition] = np.arange(-5,6)*(2*SegmentMirrorApothem + Gap)
    elif j == 1:
        # Next odd column
        PrimaryMirrorGrid[cposition] = (SegmentMirrorApothem + Gap/2.) + np.arange(-5,5)*(2*SegmentMirrorApothem + Gap)      
    elif j%2 == 0:
        # Next even column
        PrimaryMirrorGrid[cposition] = PrimaryMirrorGrid[ColumnPositions[j-2]][1:-1]
    elif j%2 == 1:
        # Nex odd column
        PrimaryMirrorGrid[cposition] = PrimaryMirrorGrid[ColumnPositions[j-2]][1:-1]

    # Add the symetrrically opposite -ve column position
    PrimaryMirrorGrid[-cposition] = PrimaryMirrorGrid[cposition]

# Create the primary mirro by shifting the Segment mirrors on to the grid

PrimaryMirrorSegmentsList = []
SegMirrorIndexDic = {} # dictionary to number the segments
i = 0
for cposition in PrimaryMirrorGrid:
    for gridpos in PrimaryMirrorGrid[cposition]:
        PrimaryMirrorSegmentsList.append( translate(SegmentMirror,xoff=cposition,yoff=gridpos) )
        SegMirrorIndexDic[cposition,gridpos] = i
        i += 1

# Take union of all segments to create Primary Mirror
print('No of Segments in Primary Mirror = {0}'.format(len(PrimaryMirrorSegmentsList)))

PrimarryMirror  = PrimaryMirrorSegmentsList[0].union(PrimaryMirrorSegmentsList[1])
for i in range(2,len(PrimaryMirrorSegmentsList)):
    PrimarryMirror = PrimarryMirror.union(PrimaryMirrorSegmentsList[i])

# Remove Bad mirror segments (if any) from the grid
BadMirrorSegmentCoords = [(+3,7),(+1,7)]  # Format: (+/- column index, row index)  # where row index 0 is the bottom most segment
for cindx,rindx in BadMirrorSegmentCoords:
    try:
        cposition = np.sign(cindx)*ColumnPositions[np.abs(cindx)]
        gridpos = PrimaryMirrorGrid[cposition][rindx]
    except (KeyError,IndexError) as e:
        print('Error: In coordinates of Bad Segment mirror {0}'.format((cindx,rindx)))
        print('Give correct coordintes..')
        continue
    BadMirrorIndx = SegMirrorIndexDic[cposition,gridpos]
    BadMirror = PrimaryMirrorSegmentsList[BadMirrorIndx]
    # Remove this Bad mirror from Primary mirror
    PrimarryMirror = PrimarryMirror.difference(BadMirror)


print('Area of Primarry Mirror = {0} m^2'.format(PrimarryMirror.area))


# Create the Pupil
PupilDia = 10.0 # meters diameter
FullPupil = Point(0,0).buffer(PupilDia/2.0)
# Central obsceuration due to WFC
CentralShadowDiaWFC = 4.5 # meter dia
WFCshadow = Point(0,0).buffer(CentralShadowDiaWFC/2.0)
PupilShape = FullPupil.difference(WFCshadow)  # Add shadow of the tracker
# Shadow of the two Tracker beams 
# Approximate sice from figure 7 in Wolf et. al. on MARS
BeamPosition = 1.8 # meters on both side
BeamWidth = 0.64 # meters
BeamLeft = Polygon([(-BeamPosition,10),(-BeamPosition,-10),(-BeamPosition-BeamWidth,-10),(-BeamPosition-BeamWidth,10)])
BeamRight = Polygon([(BeamPosition,10),(BeamPosition,-10),(BeamPosition+BeamWidth,-10),(BeamPosition+BeamWidth,10)])
SupportBeams = BeamLeft.union(BeamRight)

#Update Effective Pupil shadow
PupilShape = PupilShape.difference(SupportBeams)  # Add shadow of the tracker

def EffectivePupil(xoff=0,yoff=0, PupilShape =PupilShape, PrimaryMirror=PrimarryMirror):
    """ Returns the effective Pupil shape after shifting and intersecting with Primary mirror"""
    ShiftedPupil = translate(PupilShape,xoff=xoff,yoff=yoff)
    EffectivePupil = ShiftedPupil.intersection(PrimaryMirror)
    return EffectivePupil

xshifts = np.linspace(-3.9/2,3.9/2,1000)
AreaDrift = [EffectivePupil(xoff=x).area for x in xshifts]
plt.plot(xshifts,AreaDrift)

    

