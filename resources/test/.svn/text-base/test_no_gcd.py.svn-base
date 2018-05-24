#!/usr/bin/env python

#This whole thing is a slight modification of the millipede test fit to work for StartingTrackVeto


import sys,os.path,numpy
from I3Tray import *
from icecube import dataio,dataclasses,photonics_service,StartingTrackVeto,phys_services

def get_pxs():
    base = os.path.expandvars("$I3_DATA/photon-tables/splines/InfBareMu_mie_%s_z20a10.fits")
    if not os.path.exists(base % "abs"):
        raise errno.ENOENT((base % "abs") + " does not exist!")
    if not os.path.exists(base % "prob"):
        raise errno.ENOENT((base % "prob") + " does not exist!")
    
    return photonics_service.I3PhotoSplineService(base % "abs", base % "prob", 0.)

try:
    pxs = get_pxs()
except:
    print("Can't find full-size spline tables, skipping test")
    sys.exit(0)

testdir = os.environ["I3_TESTDATA"]
files = [ #"sim/GeoCalibDetectorStatus_2012.56063_V0.i3.gz",
         "sim/Level3_nugen_numu_IC86.2012.011069.000000_20events.i3.bz2"]
filelist = [os.path.join(testdir, f) for f in files]

pulses="InIcePulses"
fit="SplineMPE"

tray = I3Tray()

tray.AddModule("I3Reader","reader",
               FileNameList  = filelist
               )

def pulli3omgeo(frame):
    global i3omgeo
    i3omgeo=frame["I3Geometry"].omgeo
tray.Add(pulli3omgeo,"soitonlyhappensonce",Streams=[icetray.I3Frame.Geometry])

def pullbadDOMList(frame):
    global BadOMs
    print frame["BadDomsList"],frame["BadDomsListSLC"]
    BadOMs=frame["BadDomsList"]
    BadOMs.extend(frame["BadDomsListSLC"])
tray.Add(pullbadDOMList,"soitonlyhappensonce2",Streams=[icetray.I3Frame.DetectorStatus])
           
def make_n_segment_vector(frame,fit,n=1):
    if n%2==0:
        print "n=",n,"is even! Change this!"
        sys.exit(910)
    try:
        basep=frame[fit]
    except:
        print "I don't see what you're looking for"
        return True
    #shift to closest approach to 0,0,0
    origin_cap=phys_services.I3Calculator.closest_approach_position(basep,dataclasses.I3Position(0,0,0))
    #print origin_cap
    basep_shift_d=numpy.sign(origin_cap.z - basep.pos.z) *\
                  numpy.sign(basep.dir.z) *\
                  (origin_cap-basep.pos).magnitude
    #print basep_shift_d
    basep_shift_pos=basep.pos+basep.dir*basep_shift_d#basep.shift_along_track(basep_shift_d)
    #print basep_shift_pos
    basep_shift_t=basep_shift_d/basep.speed
    #print basep_shift_t
    basep.pos=basep_shift_pos
    basep.time=basep.time+basep_shift_t
    segments=[]
    segment_length=1950./n
    for idx in range(n):
        dshift=segment_length*(idx-((n-1)/2.))
        particle=dataclasses.I3Particle()
        particle.time=basep.time+(dshift/basep.speed)
        particle.pos=basep.pos+basep.dir*dshift#basep.shift_along_track(dshift)
        particle.dir=basep.dir
        particle.energy=0.01
        if n==1:
            particle.shape=particle.shape.InfiniteTrack
            particle.length=0
        else:
            particle.shape=particle.shape.ContainedTrack
            particle.length=segment_length         
        segments.append(particle)
    del frame[fit+"_"+str(n)+"_segments"]
    frame[fit+"_"+str(n)+"_segments"]=dataclasses.I3VectorI3Particle(segments)
    del segments

tray.Add(make_n_segment_vector,'make_n_segment_vector_'+fit+"_1",fit=fit,n=1)
tray.Add("StartingTrackVeto","STV",Pulses=pulses,Photonics_Service=pxs,
                Miss_Prob_Thresh=1,Fit=fit,Particle_Segments=fit+"_1_segments",
                Distance_Along_Track_Type="cherdat",Supress_Stochastics=False,Min_CAD_Dist=300)

q1=0
q2=0
def test(frame): 
    global q1,q2
    all_obs_q_map=frame[pulses+"_"+fit+"_allObsQs_1"]
    all_obs_q=0
    for map_omk in all_obs_q_map.keys():
        all_obs_q+=all_obs_q_map[map_omk]
    q1+=all_obs_q
    ps=dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulses)
    qtot=0
    i3omgeo=frame["I3Geometry"].omgeo
    f=frame[fit]
    BadOMs=frame["BadDomsList"]
    BadOMs.extend(frame["BadDomsListSLC"])
    for omk in ps.keys():
        cad = phys_services.I3Calculator.closest_approach_distance(f,i3omgeo[omk].position)
        if cad<300 and omk not in BadOMs:
            qtot+=sum([this.charge for this in ps[omk]])
    q2+=qtot
tray.AddModule(test)

try:
    tray.Execute()
except Exception:
    pass
else:
    raise Exception('should have raised exception')


