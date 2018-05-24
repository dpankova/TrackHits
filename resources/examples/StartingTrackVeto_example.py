#!/usr/bin/env python
from os import uname
from os import path as path
from os import system as system
from datetime import datetime
import sys, numpy,argparse,glob,time,os,gzip,pickle
print uname()
desc="icetray script to extract info"
parser = argparse.ArgumentParser(description=desc)
parser.add_argument('-i','--infiles',dest='infiles',
				type=str,default=[],nargs='+',
				help="[I]nfiles with frames")
parser.add_argument('-o','--outfilebase',dest='outfilebase',type=str,
				default="",help='base name for [o]utfiles')
parser.add_argument('-v','--verbose',dest='verbose',type=int,
				default=0,help='verbose? 1=True 0=False')
parser.add_argument('-p','--pulsesname',dest='pulsesname',type=str,
				default="SplitInIcePulses",help='Name of pulse series to use.')
args = parser.parse_args()


infiles=args.infiles
pulsesname=args.pulsesname
outfilebase=args.outfilebase
verbose=args.verbose
print infiles,len(infiles)
print outfilebase

from icecube import dataclasses,phys_services,dataio,icetray,photonics_service,StartingTrackVeto
from I3Tray import *

#You will need to point this to the correct spline tables
PhotonicsDir = "/cvmfs/icecube.opensciencegrid.org/data/photon-tables/"
PhotonicsSplineDirectory = os.path.join(PhotonicsDir, "splines")

inf_muon_service = photonics_service.I3PhotoSplineService(
				   amplitudetable = os.path.join(PhotonicsSplineDirectory ,"InfBareMu_mie_abs_z20a10.fits"),  ## Amplitude tables 
				   timingtable = os.path.join(PhotonicsSplineDirectory ,"InfBareMu_mie_prob_z20a10.fits"),	## Timing tables
				   timingSigma  = 0.0)
				   #maxRadius	= 125.0)

seg_muon_service = photonics_service.I3PhotoSplineService(
				   amplitudetable = os.path.join(PhotonicsSplineDirectory ,"ZeroLengthMieMuons_250_z20_a10.abs.fits"),  ## Amplitude tables 
				   timingtable = os.path.join(PhotonicsSplineDirectory ,"ZeroLengthMieMuons_250_z20_a10.prob.fits"),	## Timing tables
				   timingSigma  = 0.0)
				   #maxRadius	= 125.0)
tray = I3Tray()

tray.Add('I3Reader','reader', FilenameList=infiles)	


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
	basep_shift_pos=basep.shift_along_track(basep_shift_d)
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
		particle.pos=basep.shift_along_track(dshift)
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
	print "Put", fit+"_"+str(n)+"_segments", "in the frame"
	del segments
	 
				
def fit_check_cut(tray,name,n=1,fits_to_try=[],phot_serv=inf_muon_service,dtype="cher_dat",cad_d=125):
	for fn in fits_to_try:
		tray.Add(make_n_segment_vector,'make_n_segment_vector_'+fn+"_"+str(n),fit=fn,n=n)
		tray.Add('StartingTrackVeto','STV_'+fn+"_"+str(n),Pulses=pulsesname,Photonics_Service=phot_serv,
				Miss_Prob_Thresh=1,Fit=fn,Particle_Segments=fn+"_"+str(n)+"_segments",
				Distance_Along_Track_Type=dtype,Supress_Stochastics=True,Min_CAD_Dist=cad_d)

tray.AddSegment(fit_check_cut,"fit_check_cut_1",n=1,phot_serv=inf_muon_service,fits_to_try=["SPEFit2"],
				dtype="cher_dat")
#You can use this, but be warned I3PhotoSpline will be upset and very loud about it. It complains you didn't choose where the read table read should happen even though it is. I need to dig into this more.
#tray.AddSegment(fit_check_cut,"fit_check_cut_ns_131",n=131,phot_serv=seg_muon_service,fits_to_try=["SPEFit2"],
#				dtype="contrib_dat")

if len(outfilebase)==0:
	outfilebase="test"
tray.AddModule('I3Writer', 'writer', Filename=outfilebase+".i3.bz2",Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],DropOrphanStreams=[icetray.I3Frame.DAQ])

tray.Execute()

