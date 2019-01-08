/*
 * Copyright  (C) 2016
 * Kyle Jero
 * The Icecube Collaboration: http://www.icecube.wisc.edu
 *
 * $Id$
 *
 * @version $Revision: ?????? $
 * @date $LastChangedDate$$
 * @author $LastChangedBy$
 */

#include <TrackHits/TrackHitsUtils.h>

#include <dataclasses/I3Vector.h>
#include <dataclasses/I3Double.h>
#include <dataclasses/I3Map.h>
#include <dataclasses/I3Matrix.h>
#include <dataclasses/I3Position.h>
#include <dataclasses/I3Constants.h>
#include <dataclasses/I3TimeWindow.h>
#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/physics/I3RecoPulse.h>
#include <phys-services/I3Calculator.h>
#include <photonics-service/I3PhotonicsService.h>

#include <boost/foreach.hpp>
#include <boost/math/special_functions/sign.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/lexical_cast.hpp>
#include <vector>
#include <boost/assign/list_of.hpp>
#include <boost/make_shared.hpp>
#include <math.h>


//const double INDEX_OF_REFRACTION = I3Constants::n_ice;
const double CHER_ANGLE = I3Constants::theta_cherenkov;
const double SIN_CHER_ANGLE = sin( CHER_ANGLE );

namespace TrackHitsUtils 
{
    void ObtainEventInfo(I3FramePtr frame, 
			 I3VectorI3ParticleConstPtr segments,
			 I3PhotonicsServicePtr photonService,
			 std::vector<double> & timeEdges, 
			 std::string pulsesName,
			 std::string fitName, 
			 double mincadDist, 
			 double percent,
			 I3GeometryConstPtr geo, 
			 I3VectorOMKeyPtr badDOMs)
    {
        //Get pulses and fit fromthe frame
        I3RecoPulseSeriesMapConstPtr pulses =
	  frame->Get<I3RecoPulseSeriesMapConstPtr>(pulsesName);
	I3ParticleConstPtr fit = frame->Get<I3ParticleConstPtr>(fitName);

	
	//Get calibration errors
	I3TimeWindowSeriesMapConstPtr calibrationErrata;
	bool has_calibrationErrata = false;
	if (frame->Has("CalibrationErrata"))
	{ 
	    calibrationErrata = 
	      frame->Get<I3TimeWindowSeriesMapConstPtr>("CalibrationErrata"); 
	    has_calibrationErrata = true;
	}


	//Create Lists that will be saved
	I3MapKeyUInt keyToIndexMap;
	I3MapKeyDouble cads;//Closest Approach Distances
	I3MapKeyVectorDouble coincObsQsList;//Coincident Observed Charges
	I3MapKeyVectorDouble coincObsProbsList;//Coincident Observed Charges
	I3MapKeyDouble unnormedExQsList;//Unnormalized Expected Charges
	I3MapKeyDouble expectedBgsList;//Expected Background Charges
	unsigned int index = 0;


	//Loop over all DOMs, calculate track to DOM metrics, 
	//extract per DOM observed and expected charge information
	for (I3RecoPulseSeriesMap::const_iterator k = (*pulses).begin();
	     k != (*pulses).end(); k++)
	{
	    log_debug("-------------");
	    const OMKey& omkey = (*k).first;
	    //DC strings
	    if (omkey.GetString()>78 && omkey.GetOM()>10) 
	    { 
	        continue; 
	    }
	    //DC icecube strings
	    if (((omkey.GetString() == 35) || (omkey.GetString() == 36)
		 || (omkey.GetString() == 37)) && (omkey.GetOM()>38)) 
	    { 
	        continue; 
	    }
	    if (((omkey.GetString() == 26) || (omkey.GetString() == 27))
		&& (omkey.GetOM()>38)) 
	    { 
	        continue; 
	    }
	    if (((omkey.GetString() == 45) || (omkey.GetString() == 46))
		&& (omkey.GetOM()>38)) 
	    { 
	        continue; 
	    }
	    //Skip if it's IceTop DOM
	    if (omkey.GetString()<87 && omkey.GetOM()>60) 
	    { 
	        continue; 
	    }
	    //Gen2 Compatibility
	    if (omkey.GetString()>87 && omkey.GetOM()>80) 
	    { 
	        continue; 
	    }

	
	    //If it's a bad DOMS, skip
	    bool badOM = false;
	    for (I3VectorOMKey::const_iterator badK = badDOMs->begin();
		 badK != badDOMs->end(); badK++)
	    {
	        if ((*badK) == omkey) 
	        { 
		    badOM = true; 
		}
	    }
	    if (has_calibrationErrata)
	    {
	      for (I3TimeWindowSeriesMap::const_iterator badK = calibrationErrata->begin(); 
		   badK != calibrationErrata->end(); badK++)
	      {
                  if (badK->first == omkey) 
		  { 
		      badOM = true; 
		      break; 
		  }
	      }
	    }
	    log_debug("badOM=%d", badOM);
	    if (badOM) 
	    { 
	        continue; 
	    }
	

	    //distance metric
	    //I3Position kPos = (*k).second.position;
	    I3Position kPos = (geo->omgeo.find(omkey))->second.position;
	 
	    double cad = I3Calculator::ClosestApproachDistance((*fit), kPos);
	    log_debug("cad=%f",cad);
	    if (cad > mincadDist || std::isnan(cad)) 
	    { 
	        continue; 
	    }
	    keyToIndexMap[omkey] = index;
	    cads[omkey] = cad;

	    
	    //DOM observed charge
	    double qTot = 0;
	    std::vector<double> ts;
	    std::vector<double> qs;
	    BOOST_FOREACH(const I3RecoPulse &pulse, (*k).second)
	    {
	      ts.push_back(I3Calculator::TimeResidual((*fit), 
						      kPos, 
						      pulse.GetTime()));
	      qs.push_back(pulse.GetCharge());
	      qTot += pulse.GetCharge();
	    }
	     
	   
	    log_debug("String,OM:%d,%d",omkey.GetString(),omkey.GetOM());
	    log_debug("qTot= %f",qTot);
	

	    //DOM position for photon service
	    double kPos_x = kPos.GetX();
	    double kPos_y = kPos.GetY();
	    double kPos_z = kPos.GetZ();
	    photonService->SelectModuleCoordinates(kPos_x, 
						   kPos_y, 
						   kPos_z);


	    //Get charge info for every segment
	    double quantileSum=0;
	    I3VectorDouble sumQuantilesPerSegment;
	    I3VectorDouble sumQuantilesPerTime(timeEdges.size());
	    for (I3VectorI3Particle::const_iterator s = segments->begin();
		  s != segments->end(); s++)
	    {

	        //Set segment as a source
	        double meanPEs, ePointDistance, geoTime;
	        PhotonicsSource ps(*s);
		photonService->SelectSource(meanPEs, 
					    NULL, 
					    ePointDistance, 
					    geoTime, 
					    ps, 
					    true);
		log_debug("mean pes=%f ePointDistance=%f geoTime=%f",
			      meanPEs,ePointDistance, geoTime);

		//Get probabity quntiles for segment and DOM
		double *quantilesMuon = new double[timeEdges.size()];		
		if (meanPEs > 0)
		{
		    if ((*s).GetLength() > 0) 
		    { 
		        meanPEs *= (*s).GetLength(); 
		    }
		    photonService->GetProbabilityQuantiles(& timeEdges[0],
							   0, 
							   quantilesMuon, 
							   timeEdges.size());
		}


		//check for table oddities and build the per time 
		//yeilds from the segments sometimes the splines behave oddly,
		//you should never see a PE yield below 0 or above 1000
		//if you do something has gone arry and you should 
		//just take 0 as the PE yield to be safe
		if (meanPEs < 0 || meanPEs > 1000)
		{
		    for ( unsigned int i = 0; i < timeEdges.size(); i++ )
		    {
		        quantilesMuon[i] = 0;
		    }
		}
		else
		{
		    for ( unsigned int i = 0; i < timeEdges.size(); i++ )
		    {
		        quantilesMuon[i] = quantilesMuon[i]*meanPEs;
		    }
		}


		//build the per segment and per time yeilds
		double sumQuantilesPerSegmentTmp=0;
		for (unsigned int i = 0; i < timeEdges.size(); i++)
		{
  		    sumQuantilesPerSegmentTmp += quantilesMuon[i];
		    sumQuantilesPerTime[i] += quantilesMuon[i];
		    log_debug("muon qunatiles=%f",quantilesMuon[i]);
		}
		quantileSum+=sumQuantilesPerSegmentTmp;
		log_debug("quantile sum=%f",quantileSum);
		sumQuantilesPerSegment.push_back( sumQuantilesPerSegmentTmp );
		delete[] quantilesMuon;
	    }


	    //Find peak yeild in time so you can use it to set the 
	    //time window for comparing observed and expected charges
	    double maxQuantilePerTime=0;
	    for (unsigned int i = 0; i < timeEdges.size(); i++)
	    {
	        if (sumQuantilesPerTime[i] > maxQuantilePerTime)
		{ 
		    maxQuantilePerTime = sumQuantilesPerTime[i]; 
		}
	    }


	    
	    //Find the observed PEs which match with the expected PE 
	    //(look for yeild to drop below 1% of peak, thus the .01 below)
	    std::vector<unsigned int> gtThreshIdxs;        
	    for (unsigned int i = 0; i < timeEdges.size(); i++)
	    {
	        log_debug("sum_quantiles_pertime=%f maxQuantilePerTime=%f",
			  sumQuantilesPerTime[i],maxQuantilePerTime);
		if (sumQuantilesPerTime[i] > percent * maxQuantilePerTime)
		{ 
		    gtThreshIdxs.push_back(i); 
		}
	    }
	    log_debug("index=%d",index);
	    double minTime=0;
	    double maxTime=0;
	    if (gtThreshIdxs.size() > 0)
	    {
	        minTime=timeEdges[gtThreshIdxs[0]];
		maxTime=timeEdges[gtThreshIdxs[gtThreshIdxs.size()-1]+1];
	    }
	    log_debug("minTime=%f maxTime=%f",minTime,maxTime);


	    //Find charge in that time window
	    //double coincObsQsTmp=0;
	    std::vector<double> coincObsQsTmp;
	    std::vector<double> coincObsProbsTmp;
	    log_debug("q size=%lu",qs.size());
	    for (unsigned int  i = 0; i< qs.size(); i++)
	    {
	        log_debug("t=%f q=%f",ts[i],qs[i]);
		if (ts[i] > minTime && ts[i] < maxTime)
		{ 
		  //coincObsQsTmp+=qs[i]; 
		    coincObsQsTmp.push_back(qs[i]);
		    coincObsProbsTmp.push_back(sumQuantilesPerTime[i]);
		}
	    }
	    log_debug("coinc obs q= %f",coincObsQsTmp);
	    log_debug("coinc obs p= %f",coincObsPsTmp);
	

	    //If for some reason you encounter a nan value, 
	    //indicating something bad, record 0 and keep moving
	    if (std::isnan(quantileSum)) 
	    { 
	      std::vector<double> v =  boost::assign::list_of(0);
	      expectedBgsList[omkey] = 0; 
	      unnormedExQsList[omkey] = 0;
	      coincObsQsList[omkey] = v; 
	      coincObsProbsList[omkey] = v;
	    }
	    else 
	    { 
	      //The calculation below inserts 100 Hz of noise into every time bin
	      expectedBgsList[omkey] = (timeEdges[1] - timeEdges[0])*100*pow(10, -9);
	      unnormedExQsList[omkey] = quantileSum;
	      coincObsQsList[omkey] = coincObsQsTmp; 
	      coincObsProbsList[omkey] = coincObsProbsTmp;
	    }
	    index++;
	}


	//Book calculation pieces to the frame
	std::string namePrefix ("TrackHits_"+fitName+"_"+pulsesName);
	std::string nameSufixes[5] = { "coincObsQsList", 
				       "coincObsProbsList", 
				       "expectedQsList", 
				       "expectedBgsList",
				       "closestApDists"};
	log_debug("Names set");
	std::string sls = boost::lexical_cast<std::string>(segments->size());
	std::vector<std::string> names;
	for (unsigned int i = 0; i < 5; i++)
	{
	    std::string n = namePrefix + "_" + nameSufixes[i] + "_" +sls;
	    if (frame->Has(n))
	    { 
	        frame->Delete(n); 
	    }
	    names.push_back(n);
	    log_debug("%s",n.c_str());
	}

	log_debug("Begin booking");
	int nameIdx=0;
	frame -> Put(names[nameIdx], boost::make_shared<I3MapKeyVectorDouble>(coincObsQsList)); nameIdx++;
	frame -> Put(names[nameIdx], boost::make_shared<I3MapKeyVectorDouble>(coincObsProbsList)); nameIdx++;
	frame -> Put(names[nameIdx], boost::make_shared<I3MapKeyDouble>(unnormedExQsList)); nameIdx++;
	frame -> Put(names[nameIdx], boost::make_shared<I3MapKeyDouble>(expectedBgsList)); nameIdx++;
	frame -> Put(names[nameIdx], boost::make_shared<I3MapKeyDouble>(cads)); nameIdx++;
	  
    }
}

