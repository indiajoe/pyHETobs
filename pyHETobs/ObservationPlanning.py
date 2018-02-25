#!/usr/bin/env python
""" This script contains tools to plan observation with HET """
import os
from astropy.time import Time, TimeDelta
from astroplan import Observer
from astroplan import FixedTarget
from astropy.coordinates import SkyCoord
from pyHETobs import ObservationBlock
from TrackPupilSize import StarNotObservableError

import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
def get_target_list(targetliststring):
    """ Return a list of SkyCoordinates Object for each target in the file """
    TargetList = []
    with open(targetliststring) as targetlistfile:
        for line in targetlistfile:
            line = line.rstrip()
            if line[0] == '#':  # Skip commented out lines
                continue
            star = FixedTarget.from_name(line.split('#')[0])
            print('Identified Star as {0}'.format(star))
            TargetList.append(star)
    return TargetList


#######################################################
if __name__ == "__main__":

    McD = Observer.at_site('McDonald Observatory')

    NightOfObs = Time('2017-11-08T05:00:00')

    Targets = get_target_list('TargetList.txt')

    ObservationBlockSets = {}
    for star in Targets:
        try:
            EastObs = ObservationBlock(star.coord,east_track=True,otime=NightOfObs)
            EastObs.StarName = star.name
        except StarNotObservableError:
            EastObs = None

        try:
            WestObs = ObservationBlock(star.coord,east_track=False,otime=NightOfObs)
            WestObs.StarName = star.name
        except StarNotObservableError:
            WestObs = None

        if (EastObs is not None) or (WestObs is not None):
            ObservationBlockSets[star.name] = (EastObs,WestObs)
        else:
            print('Target :{0} not observable in East or West track'.format(star.name))

    DicOfStartAndEndTimes = {}
    # Now plot the apertures at night
    for starname in ObservationBlockSets:
        DicOfStartAndEndTimes[starname] = []
        for ObsBlock in ObservationBlockSets[starname]:
            if ObsBlock is None:
                continue
            # Calculate apperture for 100 points in the total track
            Tarray = np.linspace(ObsBlock.max_start_time,ObsBlock.max_end_time,100)
            aper = ObsBlock.EffectiveHETapperture(Tarray)
            TimeList = [ObsBlock.Transit_time + TimeDelta(t*u.second) for t in Tarray]
            DicOfStartAndEndTimes[starname].append((TimeList[0],TimeList[-1]))
            plt.plot_date([t.plot_date for t in TimeList],aper,'.',color='k',alpha=0.5)
            plt.text(ObsBlock.Transit_time.plot_date, max(aper),
                     starname+'({0})'.format(int(round(ObsBlock.Optimal_Azimuth.value))))

    
    # Print the table of observable windows
    print('='*10)
    for star in DicOfStartAndEndTimes:
        print('{0} : {1}'.format(star,DicOfStartAndEndTimes[star]))
    print('='*10)
    plt.axvline(McD.sun_set_time(NightOfObs,which='previous').plot_date,color='r')
    plt.axvline(McD.sun_rise_time(NightOfObs,which='next').plot_date,linestyle='dashed', color='r')
    plt.axvline(McD.moon_set_time(NightOfObs,which='next').plot_date,color='yellow')
    plt.axvline(McD.moon_rise_time(NightOfObs,which='previous').plot_date,linestyle='dashed',color='yellow')
    plt.xlabel('Date MM-DD HR (UT)')
    plt.ylabel(r'Pupil Area ($m^2$)')
    plt.show()
