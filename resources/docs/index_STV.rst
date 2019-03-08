.. The Icecube Collaboration
..
.. $Id$
..
.. @version $Revision$
.. @date $LastChangedDate$
.. @author Kyle Jero <kjero@icecube.wisc.edu> $LastChangedBy$

.. _TrackHits:


=================

.. toctree::
   :maxdepth: 1

   code_review
   release_notes
   
This is an Icetray Module designed to preform a look-back veto using the photon tables. An overview of the veto can be found in Note 1.  Basically, this module is designed to preform the following

* Find the first hit associated with the track you provided.
* Get the scale of the muon's light yeild with respect to a minimum ionizing muon using the track information after the first hit
* Calculate the probability of the muon recording 0 PE on the DOMs before the first hit, if there are any.

The main output is this probability, which represents the probability of missing the event if it was an incoming muon of the scale calculated. This is often referred to as p_miss. However, the information which was used in calulating the probability is also booked. This includes:

* p_miss
* The per DOM observed charges
* The per DOM observed charges that are coincident with the track
* The un-normalized per DOM expected charges
* The normalized per DOM expected charges
* The normalization constant
* The normalization llh value
* The normalization reduced llh value
* The per DOM closest approach distance (cad)
* The per DOM distance along the track (dat) of the closest approach distances
* The min and max cad dat
* The min and max cad dat DOM keys
* The distance from min to max cad dat
* The min and max cad dat standard deviation
* The per DOM Cherenkov distance (cher)
* The per DOM distance along the track (dat) of the Cherenkov distances
* The min and max cher dat
* The distance from min to max cher dat
* The min and max cher dat DOM keys
* The min and max cher dat standard deviation
* The per DOM Yeild weighted average (contributed/contrib) distance
* The per DOM distance along the track (dat) of the contributed distances
* The min and max contrib dat
* The distance from min to max contrib dat
* The min and max contrib dat DOM keys
* The min and max contrib dat standard deviation
* A map of DOM key to per DOM list index

Some configuration options are available to the user

* Parameter:"Pulses"->Name of pulse series to use (Default:SplitInIcePulses)
* Parameter:"Fit"->Name of the fit to use (Default:None, if the fit is missing the module ignore the frame.)
* Parameter:"Photonics_Service"->Photon service to use (Default:None, you have to pass the photonics service to the module so it uses the right global one, you need this to make predictions.)
* Parameter:"Particle_Segments"->Name of the particle segments to use (Default:None, if you are using segmented track photon tables you need to provide a I3VectorI3Particle object of the segments you want to use, this is asking for the objects name in the frame.)
* Parameter:"Time_Edge_Min"->Minimum of the time edges (Default:0, in order to use the photon tables you have to request the number of PE in some time bin. This is the lower edge of that bin range.)
* Parameter:"Time_Edge_Max"->Maximum of the time edges (Default:10100, in order to use the photon tables you have to request the number of PE in some time bin. This is the upper edge of that bin range.)
* Parameter:"Time_Edge_NSteps"->Number of steps for the time edges (Default:10100/100, in order to use the photon tables you have to request the number of PE in some time bin. This is the number of steps in that bin range.) **Keep in mind the more bins you reqest, the longer the calculation takes**
* Parameter:"Min_CAD_Dist"->Minimum closest approach distance considered for calculation (Default:125, DOMs within this distance of the track contribute to the calculations, outside is ignored. If you think you need a larger radius than 125 you likely should try a different track. Again this is largely a time optimization.)
* Parameter:"Supress_Stochastics"->Apply stochastic supression? (Default:false, At high energies you can remove the effect of stochastics by removing large numbers of photons too far away, and too close. See Note 1 below.)
* Parameter:"Miss_Prob_Thresh"->Miss probability threshold for rejection (Default:1, the outcome of the calculation is a probability representing the probability of missing this event if it were an incoming muon, see Note 1 for a definition. The module will cut events above this value.)
* Parameter:"Distance_Along_Track_Type"->Options:

  * "cher_dat: Distance along track from Cherenkov angle "
  * "cad_dat: Distance along track from closest approach position "
  * "contrib_dat: Distance along track from Cherenkov Angles contributing from segments ", (Default:cher_dat, the calculation finds the hit which maps back the earliest along the Fit provided. There are a few options about how to calculate this with the cher_dat and cad_dat being their normal definition. The contrib_dat is only for vetos using the segmented track table. It computes the charge weighted average and std. dev. from each segment of the Cherenkov distance to that track. The DOM which as the earliest avg. position-2(std. dev.) is the earliest. This helps with regions where scattering in the ice means none of the photons are direct and the Cherenkov angle approx. breaks down.)

* Parameter:"Geometry"->Name of geometry object to use (Default:I3DefaultName<I3Geometry>::value(), this is here for completeness, show not need to be changed.)
* Parameter:"BadDOMs"->Name of BadDOMs object to use (Default:"BadDomsList", Info from bad DOMs is ignored, makes the veto aware of holes in the detector. I assume there is a corresponding SLC list.) 

Note 1:
https://docushare.icecube.wisc.edu/dsweb/Get/Document-77688/ESTES_Update_Muon_Call_Kyle_Jero_6_6_16.pdf

An example of an IceTray script using the StartingTrackVeto can be found in the examples directory, if you'd like to test the output without building this you can borrow my copy. ::

/data/user/kjero/IceCubeSoftware/personalbuilds/ESTES_combo/build/env-shell.sh python examples/StartingTrackVeto_example.py -i /cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_2012.56063_V1.i3.gz /data/sim/IceCube/2012/filtered/level2/CORSIKA-in-ice/11499/00000-00999/Level2_IC86.2012_corsika.011499.000000.i3.bz2 -o test/data/user/kjero/IceCubeSoftware/personalbuilds/ESTES_combo/build/env-shell.sh python examples/StartingTrackVeto_example.py -i /cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_2012.56063_V1.i3.gz /data/sim/IceCube/2012/filtered/level2/CORSIKA-in-ice/11499/00000-00999/Level2_IC86.2012_corsika.011499.000000.i3.bz2 -o test
