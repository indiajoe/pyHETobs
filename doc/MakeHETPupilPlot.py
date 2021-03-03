import numpy as np
import pandas as pd
from astropy.time import Time
import astropy.units as u
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

import pyHETobs

def HETPupilPlot(Target, ObsTime, Obslocation, ParkAzimuth, ObjAltAz,
            InsertAltitude=False, InsertPupilArea=False):
    """
    INPUTS:
        Target: Astropy SkyCoord object with RA/Dec
        ObsTime: Astropy Time Object for Observation time
        Obslocation: Astropy Earth location
        ParkAzimuth: PArking Azimuth of telescope from 'STRTAZ' keyword. With u.deg unit
        ObjAltAz: AltAz AStropy object for ObsTime
        InsertAltitude: Insert Altitude plot. Default is False
        InsertPupilArea: Insert PupilArea vs Time plot. Default is False
    """
    # ParkAzimuth = Slope_Header['STRTAZ'] * u.deg
    HET = pyHETobs.TrackPupilSize.HET_Telescope(park_azimuth=ParkAzimuth)

    fig = plt.figure(constrained_layout=True)
    gs = gridspec.GridSpec(2, 3, wspace=0.7, hspace=0.3,)
    ax1 = fig.add_subplot(gs[0:2, 0:2])
    FullPupil = HET.HETPrimaryMirror.generate_combined_primary_mirror()
    EffPupil = HET.get_effective_pupil(ObjAltAz)
    PupilArea = EffPupil.area
    Coordinates = EffPupil.centroid.coords.xy
    ax1 = HET.plot_pupil(FullPupil, inp_ax=ax1, color='grey', alpha=0.2)
    ax1 = HET.plot_pupil(EffPupil, inp_ax=ax1, alpha=0.5)
    ax1.plot(Coordinates[0][0], Coordinates[1][0], 'x', markersize=12, color='r')
    ax1.xaxis.set_ticks(np.arange(-6, 7, 2))
    ax1.yaxis.set_ticks(np.arange(-6, 7, 2))
    ax1.set_xlim(-6.2, 6.2)
    ax1.set_xlabel("Mirror X (m)", size=15,  fontweight ='bold')
    ax1.set_ylim(-6.2, 6.2)
    ax1.xaxis.set_tick_params(labelsize=15)
    ax1.yaxis.set_tick_params(labelsize=15)
    ax1.set_ylabel("Mirror Y (m)", size=15,  fontweight ='bold')


    TBase = np.floor(ObsTime.jd)
    xboundary = [0.4, 1]

    if InsertAltitude:
        # from utils import find_utc_offset

        from astroplan import Observer

        delta = np.linspace(num=50, *xboundary)*24*u.hour
        times_obs = Time(TBase, format='jd')+delta
        frame_obsnight = AltAz(obstime = times_obs, location=Obslocation)
        targetaltaz_obsnight = Target.transform_to(frame_obsnight)

        ax2 = fig.add_subplot(gs[0, 2], aspect='auto')

        ax2.plot(times_obs.jd - TBase, targetaltaz_obsnight.alt ,'k', alpha=0.3)
        ax2.scatter(ObsTime.jd - TBase, ObjAltAz.alt.value, color='r', s=40)
        ax2.set_ylabel("Altitude $^{\circ}$", size=12,  fontweight ='bold')
        ax2.xaxis.set_tick_params(labelsize=12)
        ax2.yaxis.set_tick_params(labelsize=12)
        ax2.tick_params(top=False, right=False)
        ax2.xaxis.set_ticks(np.linspace(xboundary[0], xboundary[1], 3))
        ax2.set_ylim(40, 70)
        ax2.set_xlim(*xboundary)


    if InsertPupilArea:
        ax3 = fig.add_subplot(gs[1, 2])
        ax3.scatter(ObsTime.jd - TBase, PupilArea, color='r', s=40)

        ax3.set_ylim(30, 51)
        ax3.set_xlim(*xboundary)
        ax3.xaxis.set_tick_params(labelsize=12)
        ax3.yaxis.set_tick_params(labelsize=12)
        ax3.set_ylabel("Pupil Area (m$^2$)", size=12,  fontweight ='bold')
        ax3.set_xlabel("JD", size=12,  fontweight ='bold')
        ax3.tick_params(top=False, right=False)
        ax3.xaxis.set_ticks(np.linspace(xboundary[0], xboundary[1], 3))

    plt.tight_layout()

    return fig, ax1, Coordinates, PupilArea
