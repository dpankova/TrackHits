/*
 * copyright  (C) 2016
 * Kyle Jero
 * The Icecube Collaboration: http://www.icecube.wisc.edu
 *
 * $Id$
 *
 * @version $Revision: ?????? $
 * @date $LastChangedDate$$
 * @author $LastChangedBy$
 */


#include <icetray/I3ConditionalModule.h>
#include <icetray/I3Bool.h>
#include <icetray/OMKey.h>
#include <dataclasses/geometry/I3Geometry.h>
#include <dataclasses/status/I3DetectorStatus.h>
#include <dataclasses/I3Double.h>
#include <dataclasses/I3Vector.h>
#include <photonics-service/I3PhotonicsService.h>

#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>
#include <math.h>
#include <TrackHits/TrackHits.h>
#include <TrackHits/TrackHitsUtils.h>

I3_MODULE(TrackHits);

    TrackHits::TrackHits(const I3Context &context)
: I3ConditionalModule(context)
{
    AddParameter("Pulses", 
		 "Name of pulse series to use", 
		 "SplitInIcePulses");
    AddParameter("Fit", 
		 "Name of the fit to use",
		 "");
    AddParameter("Photonics_Service", 
		 "Photon service to use");

    AddParameter("Particle_Segments", 
		 "Name of the particle segments to use",
		 "");
    AddParameter("Time_Edge_Min", 
		 "Minimum of the time edges", 
		 0);
    AddParameter("Time_Edge_Max", 
		 "Maximum of the time edges", 
		 10100 );
    AddParameter("Time_Edge_NSteps", 
		 "Number of steps for the time edges", 
		 10100/100 );
    AddParameter("Min_CAD_Dist",
		 "Minimum closest approach distance \
                  considered for calculation", 
		 125);
    AddParameter("Percent", 
		 "Percent of max for time interval", 
		 0.01);
    AddParameter("Geometry", 
		 "Name of geometry object to use",
		 I3DefaultName<I3Geometry>::value());
    AddParameter("DeepCore",
                 "Look at Deepcore DOMS or only IC DOMS",
                 false);
    AddParameter("BadDOMs",
		 "Name of BadDOMs object to use\
                 (I assume there is a corresponding SLC list)",
		 "BadDomsList");
}

void TrackHits::Configure()
{
    GetParameter("Pulses", pulsesName_);
    GetParameter("Fit", fitName_);
    GetParameter("Photonics_Service", photonicsService_);
    GetParameter("Particle_Segments", particleSegmentsName_);
    GetParameter("Time_Edge_Min", timeEdgeMin_);
    GetParameter("Time_Edge_Max", timeEdgeMax_);
    GetParameter("Time_Edge_NSteps", timeEdgeNSteps_);
    GetParameter("Min_CAD_Dist", minCADDist_);
    GetParameter("Percent", percent_);
    GetParameter("DeepCore", deepcore_);
    GetParameter("Geometry", geoName_);
    GetParameter("BadDOMs", badDOMsName_);
}

void TrackHits::Geometry(I3FramePtr frame)
{
    geo_ = frame->Get<I3GeometryConstPtr>(geoName_);
    if (!geo_)
    {
        log_fatal("No I3Geometry object found with name \"%s\"",
		  geoName_.c_str());
    }
    PushFrame(frame);
}

void TrackHits::DetectorStatus(I3FramePtr frame)
{
    if (badDOMsName_ != "")
    {
        std::string badDOMsNameSLC=badDOMsName_+"SLC";
        I3VectorOMKeyConstPtr badDOMsAll = 
	  frame->Get<I3VectorOMKeyConstPtr>(badDOMsName_);
        I3VectorOMKeyConstPtr badDOMsSLC = 
	  frame->Get<I3VectorOMKeyConstPtr>(badDOMsNameSLC);

        if (!badDOMsAll)
	{
            log_fatal("No object found with name \"%s\"",
		      badDOMsName_.c_str());
	}
        if (!badDOMsSLC)
	{
	  log_fatal("No object found with name \"%s\"",
                    badDOMsNameSLC.c_str());
	}
        for (I3VectorOMKey::const_iterator i = 
	       badDOMsAll->begin(); i != badDOMsAll->end(); i++)
	{
            badDOMs_.push_back(*i);
        }
        for (I3VectorOMKey::const_iterator i = 
	       badDOMsSLC->begin(); i != badDOMsSLC->end(); i++)
        {
            badDOMs_.push_back(*i);
        }
    }
    PushFrame(frame);
}

void TrackHits::Physics(I3FramePtr frame)
{
    log_debug("%s",fitName_.c_str());
    if (!frame->Has(fitName_)) 
    { 
        PushFrame(frame); 
    }
    else
    {
        I3VectorI3ParticleConstPtr segments =
	    frame->Get<I3VectorI3ParticleConstPtr>(particleSegmentsName_);

        std::vector<double> timeEdges;
        for (double step = timeEdgeMin_; step <= timeEdgeMax_; 
                step += (timeEdgeMax_ - timeEdgeMin_)/timeEdgeNSteps_)
        {
            timeEdges.push_back(step);
        }


        //Given a track, photonics service, and pulse series, find the per DOM expectations and observations. 
        //Where the expectations and observations match define the region of the observed muon.
        //Make the DOM -> Track mappings for later use
        TrackHitsUtils::ObtainEventInfo(frame, 
					segments, 
					photonicsService_,
					timeEdges, 
					pulsesName_, 
					fitName_,
					minCADDist_, 
					percent_,
					deepcore_,
					geo_, 
					boost::make_shared<I3VectorOMKey>(badDOMs_) );

	PushFrame(frame);
    }
}
