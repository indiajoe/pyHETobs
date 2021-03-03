#!/usr/bin/env python
""" This script is to calculate the movement of HET tracker during a source observation """

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import numpy as np
from scipy import optimize
import sys
from scipy import interpolate
import matplotlib.pyplot as plt
from astroplan import Observer
from astroplan import FixedTarget

from . import HETparams


def find_HET_optimal_azimuth(StarCoo,ZenithCrossTime,find_east_track=True):
    """ 
    Input:
         StarCoo : Star coordinate objects
         ZenithCrossTime: Time as which Star will cross the Zenith
         find_east_track:  True will find the eastern side track, False will find western side track.

    Returns:
         AltAz_ofStar.az: the optimal azimuth to park the telescope to observe star, 
         distance_to_HETAlt_min: minimum delta distance to HET's altitude
         time_of_minimum_distance: time of closes approach to HET's altitude
    """
    McDobserver = Observer.at_site('McDonald Observatory')
    Target = FixedTarget(name='Star',coord=StarCoo)
    try:
        alt55crosstime = McDobserver.target_rise_time(ZenithCrossTime,Target,horizon=HETparams.HET_FixedAlt,which='nearest')
    except TypeError:
        minimal_time = 0 *u.hour
    else:
        if alt55crosstime.value < 0:
            minimal_time = 0 *u.hour
        else:
            minimal_time = abs(ZenithCrossTime - alt55crosstime)

    if find_east_track:
        time_of_minimum_distance = ZenithCrossTime - minimal_time
    else:
        time_of_minimum_distance = ZenithCrossTime + minimal_time

    AltAz_ofStar =  StarCoo.transform_to(AltAz(obstime=time_of_minimum_distance,location=HETparams.McDonaldObservatory))
    distance_to_HETAlt = np.abs((AltAz_ofStar.alt - HETparams.HET_FixedAlt).value)
    return AltAz_ofStar.az, distance_to_HETAlt, time_of_minimum_distance




def start_and_end_of_tracktime(TransitTime,StarCoo,TelescopeAltAz,TrackerRad=None):
    """ Returns the starting time and ending time corresponding to star entering and leaving the Tracker Radius circle """
    TrackerRad = TrackerRad or HETparams.Tracker_Radius  # Set default to global HETparams.Tracker_Radius

    def distance_to_TelescopeTrackerCircle(deltatime,finding_start=True):
        """ deltatime : (+ve float) in units of hours. finding_start: True to find the startinng time, False for finding end time """
        if finding_start :
            time = TransitTime - deltatime*u.hour
        else:
            time = TransitTime + deltatime*u.hour

        AltAz_ofStar =  StarCoo.transform_to(AltAz(obstime=time,location=HETparams.McDonaldObservatory))

        distance_to_HETAltAz = TelescopeAltAz.separation(AltAz_ofStar)

        return np.abs(np.abs(distance_to_HETAltAz.value) - TrackerRad.value)

    # First find starting interestion in timewith the tracker circle
    start_timedelta = optimize.minimize(distance_to_TelescopeTrackerCircle, 0.5,args=(True,),bounds=[(0,12)],tol=2./3600.)
    # Find end of track time
    end_timedelta = optimize.minimize(distance_to_TelescopeTrackerCircle, 0.5,args=(False,),bounds=[(0,12)],tol=2./3600.)

    Track_StartTime = TransitTime - start_timedelta.x[0]*u.hour
    Track_EndTime = TransitTime + end_timedelta.x[0]*u.hour

    return Track_StartTime, Track_EndTime


def Tracker_Xoff_Yoff_function(Transit_time, Track_StartTime, Track_EndTime, TelescopePark_AltAz, StarCoo):
    """ 
    Returns: two functions of time, namely Xoff_calculator and Yoff_calculator 
    They return Xoff and Yoff in units of meters for and input offset time from TransitTime.
    input time should be in seconds to thie function.
    Transit_time defines the zero epoch of time for the functions.
    """
    Deg2Meter_Scale = 1/4.37997  # Taken from operaz3.f code


    # No we sample the time into No_of_epochs time epochs for generating an interpolation function
    No_of_epochs = 100
    Track_Duration = Track_EndTime - Track_StartTime
    Epochs_array = Track_StartTime + np.linspace(0,1,No_of_epochs) * Track_Duration
    
    frames_DuringTheTrack = AltAz(obstime=Epochs_array, location=HETparams.McDonaldObservatory)
    AltAzs_ofStarDuringTrack =  StarCoo.transform_to(frames_DuringTheTrack)

    Xoff_values = []
    for altaz in AltAzs_ofStarDuringTrack:
        # PoleCorrection to subtract and remove any sudden jump at the north pole direction when azimuth of the star suddenly change from 0 to 360 degree in a track.
        if (altaz.az > 270*u.deg) and (TelescopePark_AltAz.az < 90*u.deg): 
            PoleCorrection = -360*u.deg
        elif (altaz.az < 90*u.deg) and (TelescopePark_AltAz.az > 270*u.deg): 
            PoleCorrection = 360*u.deg
        else:
            PoleCorrection = 0*u.deg

        Xoff_values.append( ((altaz.az + PoleCorrection - TelescopePark_AltAz.az)*np.cos(altaz.alt)).value * Deg2Meter_Scale )

    Yoff_values = [(altaz.alt - TelescopePark_AltAz.alt).value * Deg2Meter_Scale for altaz in AltAzs_ofStarDuringTrack]
    time_seconds = [timedelta.sec for timedelta in  Epochs_array - Transit_time]
    
    Xoff_calculator = interpolate.interp1d(time_seconds, Xoff_values, kind='cubic',bounds_error=True)
    Yoff_calculator = interpolate.interp1d(time_seconds, Yoff_values, kind='cubic',bounds_error=True)

    return Xoff_calculator, Yoff_calculator

def pupil_Xoff_Yoff_function(Transit_time, Track_StartTime, Track_EndTime, TelescopePark_AltAz, StarCoo):
    """ 
    Returns: two functions of time, namely pXoff_calculator and pYoff_calculator 
    They return illuminated pupil's Xoff and Yoff positions on primry mirror in units of meters for and input offset time from TransitTime.
    input time should be in seconds to this function.
    Transit_time defines the zero epoch of time for the functions.
    """
    # The optical axis line of the WFC on tracker, centre of the curvature of primary and the star, should always fall in a straight line. 
    # The pupil footprint simply shift by the radius of curvature of primary mirror multiplied by the offset angle of the star!

    # No we sample the time into No_of_epochs time epochs for generating an interpolation function
    No_of_epochs = 100
    Track_Duration = Track_EndTime - Track_StartTime
    Epochs_array = Track_StartTime + np.linspace(0,1,No_of_epochs) * Track_Duration
    
    frames_DuringTheTrack = AltAz(obstime=Epochs_array, location=HETparams.McDonaldObservatory)
    AltAzs_ofStarDuringTrack =  StarCoo.transform_to(frames_DuringTheTrack)

    pXoff_values, pYoff_values = zip(*[pupil_Xoff_Yoff_foraltaz(altaz,TelescopePark_AltAz) for altaz in AltAzs_ofStarDuringTrack])

    time_seconds = [timedelta.sec for timedelta in  Epochs_array - Transit_time]
    
    pXoff_calculator = interpolate.interp1d(time_seconds, pXoff_values, kind='cubic',bounds_error=True)
    pYoff_calculator = interpolate.interp1d(time_seconds, pYoff_values, kind='cubic',bounds_error=True)
    #Save the max star and end time in seconds of the interpolation function boundry
    pXoff_calculator.max_end_time = time_seconds[-1]
    pXoff_calculator.max_start_time = time_seconds[0]
    return pXoff_calculator, pYoff_calculator

def pupil_Xoff_Yoff_foraltaz(altaz,TelescopePark_AltAz):
    """
    Returns the instantanious Xoff and Yoff in pupil while telescope parked at `TelescopePark_AltAz` is staring at `altaz`
    """
    # The optical axis line of the WFC on tracker, centre of the curvature of primary and the star, should always fall in a straight line. 
    # The pupil footprint simply shift by the radius of curvature of primary mirror multiplied by the offset angle of the star!


    # PoleCorrection to subtract and remove any sudden jump at the north pole direction when azimuth of the star suddenly change from 0 to 360 degree in a track.
    if (altaz.az > 270*u.deg) and (TelescopePark_AltAz.az < 90*u.deg): 
        PoleCorrection = -360*u.deg
    elif (altaz.az < 90*u.deg) and (TelescopePark_AltAz.az > 270*u.deg): 
        PoleCorrection = 360*u.deg
    else:
        PoleCorrection = 0*u.deg
    pXoff = ((altaz.az + PoleCorrection - TelescopePark_AltAz.az)*np.cos(altaz.alt)).radian * HETparams.RadiusOfCurvatureHETprimary
    pYoff = -(altaz.alt - TelescopePark_AltAz.alt).radian * HETparams.RadiusOfCurvatureHETprimary
    return pXoff, pYoff
    

#######################################################
if __name__ == "__main__":
    StarName = '2MASS J23062928-0502285'  # TRAPPIST-1
    StarCoo = SkyCoord.from_name(StarName)

    CurrentTime = Time.now()
    # First calculate the Zenith corssing time of the star
    HA_toStar = CurrentTime.sidereal_time('apparent',longitude=HETparams.McDonaldObservatory.lon) - StarCoo.ra
    ZenithCrossTime = CurrentTime - HA_toStar.hour *u.hour 

    print('{0} will cross zenith at {1}'.format(StarName,ZenithCrossTime))

    print('Calculating optimal Azimuth and transit time..')
    Optimal_Azimuth, min_alt_delta, Transit_time = find_HET_optimal_azimuth(StarCoo,ZenithCrossTime)

    # If the minimium altitude difference is beyoud tracker radius, exit
    if min_alt_delta > HETparams.Tracker_Radius.value :
        print('This star {0} is never observable with HET'.format(StarName))
        sys.exit(1)

    print('Optimal Azimuth to park telescope is at: {0}'.format(Optimal_Azimuth))


    TelescopePark_AltAz = AltAz(az=Optimal_Azimuth, alt=HETparams.HET_FixedAlt)

    print('Calculating Track begin and End time..')
    Track_StartTime, Track_EndTime = start_and_end_of_tracktime(Transit_time,StarCoo,TelescopePark_AltAz,HETparams.Tracker_Radius)
    print('Duration of Track length = {0.sec} seconds'.format(Track_EndTime - Track_StartTime))

    print("Calculating the Tracker's tracks..")
    Xoff_calculator, Yoff_calculator = Tracker_Xoff_Yoff_function(Transit_time, Track_StartTime, Track_EndTime, TelescopePark_AltAz, StarCoo)


    ################################################################################
    # Ploting of toilet seat of target accessibility locus graph
    ################################################################################

    DecRange = np.concatenate([np.arange(-15,5,0.5), np.linspace(5,55,50/2), np.arange(55,75,0.5)])
    StarList = [SkyCoord(0.0*u.deg, dec*u.deg, frame='icrs') for dec in DecRange]

    Dec2plotList = []
    HAstart2plotList = []
    HAend2plotList = []


    for dec,StarCoo in zip(DecRange,StarList):
        # First calculate the Zenith corssing time of the star
        HA_toStar = CurrentTime.sidereal_time('apparent',longitude=HETparams.McDonaldObservatory.longitude) - StarCoo.ra
        ZenithCrossTime = CurrentTime - HA_toStar.hour *u.hour 

        print('Calculating optimal Azimuth and transit time..')
        Optimal_Azimuth, min_alt_delta, Transit_time = find_HET_optimal_azimuth(StarCoo,ZenithCrossTime)

        # If the minimium altitude difference is beyond tracker radius, skip
        if min_alt_delta > HETparams.Tracker_Radius.value :
            print('This star {0} is never observable with HET'.format(StarName))
            continue

        print('Optimal Azimuth to park telescope is at: {0}'.format(Optimal_Azimuth))


        TelescopePark_AltAz = AltAz(az=Optimal_Azimuth, alt=HETparams.HET_FixedAlt)

        print('Calculating Track begin and End time..')
        Track_StartTime, Track_EndTime = start_and_end_of_tracktime(Transit_time,StarCoo,TelescopePark_AltAz,HETparams.Tracker_Radius)
        print('Duration of Track length = {0.sec} seconds'.format(Track_EndTime - Track_StartTime))

        HAstart = Track_StartTime.sidereal_time('apparent',longitude=HETparams.McDonaldObservatory.longitude) - StarCoo.ra
        HAend = Track_EndTime.sidereal_time('apparent',longitude=HETparams.McDonaldObservatory.longitude) - StarCoo.ra

        Dec2plotList.append(dec)
        HAstart2plotList.append(np.abs(HAstart.value))
        HAend2plotList.append(np.abs(HAend.value))


    #Mask out bad points in the edges of dec in EndTrack which is not valid 
    Mask = np.array(HAend2plotList) > 12
    # Also Mask out initial few points which becme 24 due to error/noise
    Mask[0:4] = False
    plt.plot(24-np.array(HAend2plotList)[Mask],np.array(Dec2plotList)[Mask],'k')
    plt.plot(24-np.array(HAstart2plotList),Dec2plotList,'k')
    plt.plot(np.array(HAend2plotList)[Mask]-24,np.array(Dec2plotList)[Mask],'k')
    plt.plot(np.array(HAstart2plotList)-24,Dec2plotList,'k')
    plt.xlabel('HA (hours)')
    plt.ylabel('Declination (deg)')
    plt.grid()
    plt.savefig('HET_target_accessibility_locus.png')
    plt.show()
    ################################################################################
