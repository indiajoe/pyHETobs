""" This is all the static parameter values for the HET telscope model """

#################################
# Related to HET primary mirror
#################################

RadiusOfCurvatureHETprimary = 26.164 # meters  Radius of curvature of HET primary mirror

SegmentMirrorApothem = 0.5  # meters # Apothem of hexagonal segment mirror 

# Gaps between mirrors for creating the Primary mirror Honey comb grid
Gap = 0.011  # 6.2 to 15.8 mm # 11 mm on average

# list of Bad Mirror segments which has dropped out of the allignemnt
BadMirrorSegmentCoords = [] #[(+3,7),(+1,7)]  # Format: (+/- column index, row index)  # where row index 0 is the bottom most segment


#################################
# Related to WFC and pupil image
#################################

PupilDia = 10.0 # meters # Outer diameter of the WFC pupil on Primary mirror
CentralShadowDiaWFC = 4.5 # meter # Diameter of the central obscuration.
# Support beams holding the WFC on top.
BeamPosition = 1.3 # meters on both side # position of the two beam's shadow from the center of the pupil
BeamWidth = 0.29 # meters  # width of the beam shadow


#################################
# Related to X-Y tracker
#################################
import astropy.units as u
## offset radius upto which Tracker can track
Tracker_Radius = 8.4 *u.deg

#################################
# Observatory Location and Telescope altitude
#################################
from astropy.coordinates import EarthLocation

McDonaldObservatory = EarthLocation.of_site(u'mcdonald')
HET_FixedAlt = 55. *u.deg


##############################
# CCAS Tower shadow numbers
##############################
CCAS_azimuth = 68.593 *u.deg

TelescopePierToCC = 16.433  #meters, Horizontal distance from telescope pier to CC

CCASDia = 4.191  #meters, Diamter of CCAS dome

CCAScenter = 0, -1.518  #meters, Position of the Dome circle center with respect to CC

CCASwidth = 2.134   #meters,  Width of the CCAS tower cylinder
