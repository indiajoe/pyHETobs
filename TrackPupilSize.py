""" This script contains functions to calculate the variation of effective HET pupil size as a function of time """

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, FK5
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

import HETparams
import HET_Tracker 
import HET_pupil

class StarNotObservableError(ValueError):
    """ To be raised when the star is not observable or not identified """

class ObservationBlock(object):
    """ One observation Block on HET """
    def __init__(self,StarName,east_track=True,otime=None):
        """ Define observation block for the input StarName
        Input:
              StarName: (str or SkyCoords) Name of the Star to reolve from Simbad or Direct SkyCoord object.
              east_track : True to observe during east track, False for observe in west track.
              otime : (optional) Day of observation. 
        """
        # Find the Coordinates of the user input Star Name if it is a string
        if isinstance(StarName,str):
            self.StarName = StarName
            self.StarCoo =  SkyCoord.from_name(StarName)
            print('Identified Star coordinates as {0}'.format(self.StarCoo))
        elif isinstance(StarName,SkyCoord):
            self.StarName = 'your star'
            self.StarCoo =  StarName
        else:
            print('Error: Unknown object provided as StarName: {0}'.format(StarName))

        # East or West track
        self.east_track = east_track

        # Approximate time of observation user wanted.
        if otime is None:
            self.ApproxTime = Time.now()
        else:
            self.ApproxTime = otime

        # Calculte Zenith crossing time
        self.ZenithCrossTime = self.calculate_ZenithCrossTime(self.ApproxTime)

        # Calculating optimal Azimuth and transit time..
        self.Optimal_Azimuth, self.Transit_time = self.calculate_optimalAz_TrackTransit(self.StarCoo,
                                                                                        self.ZenithCrossTime)
        self.TelescopePark_AltAz = AltAz(az=self.Optimal_Azimuth, alt=HETparams.HET_FixedAlt)
        
        #Calculating Track begin and End time..
        self.Track_StartTime, self.Track_EndTime = HET_Tracker.start_and_end_of_tracktime(self.Transit_time,self.StarCoo,self.TelescopePark_AltAz)

        print('Duration of Track length = {0.sec} seconds'.format(self.Track_EndTime - self.Track_StartTime))
        #Calculating the pupil's offsef on Primary  tracks..
        self.pXoff_calculator, self.pYoff_calculator = HET_Tracker.pupil_Xoff_Yoff_function(self.Transit_time, self.Track_StartTime, self.Track_EndTime, self.TelescopePark_AltAz, self.StarCoo)
        # save the max star and end time of the trackes in seconds w.r.t Track transit
        self.max_end_time = self.pXoff_calculator.max_end_time
        self.max_start_time = self.pXoff_calculator.max_start_time

        # Now create the primary mirror and pupil of HET
        print ('Initialising HET Primary Mirror and Pupil..')
        self.HETPrimaryMirror = HET_pupil.PrimaryMirror()
        self.HETPupil = HET_pupil.Pupil()


    def calculate_ZenithCrossTime(self,ApproxTime,StarCoo = None):
        """ Returns the time at which star will cross Zenith """
        StarCoo = StarCoo or self.StarCoo
        HA_toStar = ApproxTime.sidereal_time('apparent',longitude=HETparams.McDonaldObservatory.longitude) - StarCoo.transform_to(FK5(equinox=ApproxTime)).ra
        ZenithCrossTime = ApproxTime - HA_toStar.hour *u.hour 
        print('{0} will cross zenith at {1}'.format(self.StarName,ZenithCrossTime))

        return ZenithCrossTime

    def calculate_optimalAz_TrackTransit(self,StarCoo=None,ZenithCrossTime=None,find_east_track=None):
        """ Calculates the optimal azimuth to park the telescope during this observation block """
        StarCoo = StarCoo or self.StarCoo
        ZenithCrossTime = ZenithCrossTime or self.ZenithCrossTime
        find_east_track = find_east_track or self.east_track

        print('Calculating optimal Azimuth and transit time..')
        Optimal_Azimuth, min_alt_delta, Transit_time = HET_Tracker.find_HET_optimal_azimuth(StarCoo,ZenithCrossTime,find_east_track)

        # If the minimium altitude difference is beyoud tracker radius, exit
        if min_alt_delta > HETparams.Tracker_Radius.value :
            print('ERROR: I am sorry to break this sad news to you..')
            print('This star {0} is never observable with HET'.format(self.StarName))
            print('Closest approch to HET elevation is {0}'.format(min_alt_delta))
            raise StarNotObservableError('Star {0} is outside HET observable region in sky'.format(self.StarName))

        TrackName = 'East' if find_east_track else 'West'
        print('Optimal Azimuth to park telescope is at: {0}'.format(Optimal_Azimuth))
        print('Star will transit the {0} track window at {1}'.format(TrackName, Transit_time))

        return Optimal_Azimuth, Transit_time

    def EffectiveHETapperture(self,epochs):
        """ Returns an array of effective HET telescope area during the eposure time at inpot epochs.
        Input:
             epochs: (float array) Array of seconds epochs during this observation (unit: seconds)
                    Track Transit_time defines the zero epoch of this input.
        Returns:
             EffArea: (numpy array) Array of the effective telescope apperuture
        """
        
        EffArea = []
        for t in epochs:
            pXoff = self.pXoff_calculator(t)
            pYoff = self.pYoff_calculator(t)
            area = self.HETPupil.EffectiveActivePupil(self.HETPrimaryMirror.PrimaryMirror, xoff=pXoff,yoff=pYoff).area
            EffArea.append(area)

        return np.array(EffArea)

    def plot_effectiveHETapperture(self,epoch,alpha=0.4,color="g",inp_ax=None,outputfile=None):
        """ Plots the effective apperture of the telscope at the input epoch """
        pXoff = self.pXoff_calculator(epoch)
        pYoff = self.pYoff_calculator(epoch)
        EffPupil = self.HETPupil.EffectiveActivePupil(self.HETPrimaryMirror.PrimaryMirror, xoff=pXoff,yoff=pYoff)
        
        if inp_ax is None:
            fig, ax = plt.subplots(figsize=(8, 8))
        else:
            ax = inp_ax

        for polygon in EffPupil:
            mpl_poly = Polygon(np.array(polygon.exterior), facecolor=color, lw=0, alpha=alpha)
            ax.add_patch(mpl_poly)

        ax.set_xlim(-6,6)
        ax.set_ylim(-6,6)
        ax.set_aspect('equal')
        if (inp_ax is None) and (outputfile is not None):
            fig.savefig(outputfile)
        return ax









#######################################################
if __name__ == "__main__":
    ################################################################################
    # Plot of how the pupil area changes during the track at various dec values
    ################################################################################
    DecRange = np.linspace(-15,75,45)
    StarList = [SkyCoord(0.0*u.deg, dec*u.deg, frame='icrs') for dec in DecRange]
    
    EffPupilAreaDic = {} # Dictionary to store results

    for dec,StarCoo in zip(DecRange,StarList):
        print('Declination : {0}'.format(dec))
        try:
            StarObsBlock = ObservationBlock(StarCoo)
        except StarNotObservableError:
            pass
        else:
            # Calculate apperture for 2000 points in the total track
            t = np.linspace(StarObsBlock.max_start_time,StarObsBlock.max_end_time,2000)
            aper = StarObsBlock.EffectiveHETapperture(t)
            EffPupilAreaDic[dec] = {'t':t, 'aper':aper }
    
    for dec in sorted(EffPupilAreaDic.keys()):
        plt.plot(EffPupilAreaDic[dec]['t'], EffPupilAreaDic[dec]['aper'],label=str(round(dec,0)),alpha=0.8)
    plt.legend()
    plt.grid(True)
    plt.xlabel('Track time centered on transit (seconds)')
    plt.ylabel('Effective Aperture (m^2)')
    plt.savefig('EffectiveHETAreaAcrossDec.png')
    plt.show()
