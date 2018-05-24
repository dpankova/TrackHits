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

#include <StartingTrackVeto/StartingTrackVetoUtils.h>

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
#include <boost/make_shared.hpp>
#include <math.h>

//const double INDEX_OF_REFRACTION = I3Constants::n_ice;
const double CHER_ANGLE = I3Constants::theta_cherenkov;
const double SIN_CHER_ANGLE = sin( CHER_ANGLE );

namespace StartingTrackVetoUtils {
void ObtainEventInfo( I3FramePtr frame, I3VectorI3ParticleConstPtr segments,
        I3PhotonicsServicePtr photonService,
        std::vector<double> & timeEdges, std::string pulsesName,
        std::string fitName, double mincadDist,
        I3GeometryConstPtr geo, I3VectorOMKeyPtr badDOMs )
{
    I3RecoPulseSeriesMapConstPtr pulses =
        frame->Get<I3RecoPulseSeriesMapConstPtr>(pulsesName);
    I3ParticleConstPtr fit = frame->Get<I3ParticleConstPtr>(fitName);
    I3TimeWindowSeriesMapConstPtr calibrationErrata;
    bool has_calibrationErrata = false;
    if ( frame->Has( "CalibrationErrata" ) ){ calibrationErrata = frame->Get<I3TimeWindowSeriesMapConstPtr>("CalibrationErrata"); has_calibrationErrata = true;}
    I3MapKeyUInt keyToIndexMap;
    I3MapKeyDouble allObsQs;//All Observed Charges
    I3MapKeyDouble cads;//Closest Approach Distances
    I3MapKeyDouble caddats;//Closest Approach Distance Distance Along Tracks
    I3MapKeyDouble caddatsStdDevs;
    I3MapKeyDouble cherdats;//Cherenkov Distance Along Tracks
    I3MapKeyDouble cherdatsStdDevs;
    I3MapKeyDouble contribdats;//Average Contributed Distance Along Tracks
    I3MapKeyDouble contribdatsStdDevs;
    I3VectorDouble coincObsQsList;//Coincident Observed Charges
    I3VectorDouble unnormedExQsList;//Unnormalized Expected Charges
    I3VectorDouble expectedBgsList;//Expected Background Charges
    unsigned int index = 0;
    //Loop over all DOMs, calculate track to DOM metrics, extract per DOM observed and expected charge information
    for ( I3OMGeoMap::const_iterator k = geo->omgeo.begin();
            k != geo->omgeo.end(); k++ )
    {
        log_debug("-------------");
        const OMKey& omkey = k->first;
        I3VectorDouble sumQuantilesPerSegment;
        I3VectorDouble sumQuantilesPerTime ( timeEdges.size() );
        if ( omkey.GetString()<87 && omkey.GetOM()>60 ) { continue; }
        if ( omkey.GetString()>87 && omkey.GetOM()>80 ) { continue; }//Gen2 Compatibility
        bool badOM = false;
        for ( I3VectorOMKey::const_iterator badK = badDOMs->begin();
                badK != badDOMs->end(); badK++ )
        {
            if ( (*badK) == omkey ) { badOM = true; }
        }
        if ( has_calibrationErrata )
        {
            for ( I3TimeWindowSeriesMap::const_iterator badK = calibrationErrata->begin();
                    badK != calibrationErrata->end(); badK++ )
            {
                if ( badK->first == omkey ) { badOM = true; break; }
            }
        }
        log_debug( "badOM=%d", badOM );
        if ( badOM ) { continue; }
        //distance metric
        I3Position kPos = (*k).second.position;
        double cad = I3Calculator::ClosestApproachDistance( (*fit), kPos );
        log_debug("cad=%f",cad);
        if ( cad > mincadDist || std::isnan( cad ) ) { continue; }
        keyToIndexMap[omkey] = index;
        cads[omkey] = cad;
        //DOM observed charge
        double qTot = 0;
        std::vector<double> ts;
        std::vector<double> qs;
        for (I3RecoPulseSeriesMap::const_iterator i = (*pulses).begin();
                i != (*pulses).end(); i++)
        {
            if ( (*i).first == omkey )
            {
                BOOST_FOREACH(const I3RecoPulse &pulse, i->second)
                {
                    ts.push_back(I3Calculator::TimeResidual( (*fit), kPos, pulse.GetTime() ) );
                    qs.push_back(pulse.GetCharge());
                    qTot += pulse.GetCharge();
                }
            }
        }
        allObsQs[omkey] = qTot;
        log_debug("String,OM:%d,%d",omkey.GetString(),omkey.GetOM());
        log_debug("qTot= %f",qTot);
        //distance metric
        I3Position cap = I3Calculator::ClosestApproachPosition( (*fit), kPos );
        double caddat = ( (*fit).GetPos() - cap ).Magnitude();
        caddat *= boost::math::sign( cap.GetZ() - (*fit).GetPos().GetZ() ) * 
            boost::math::sign( (*fit).GetDir().GetZ() );
        caddats[omkey] = caddat;
        double cherdat = caddat - ( cad / tan( CHER_ANGLE ));
        cherdats[omkey] = cherdat;
        //DOM expected charge
        double kPos_x = kPos.GetX();
        double kPos_y = kPos.GetY();
        double kPos_z = kPos.GetZ();
        photonService->SelectModuleCoordinates( kPos_x, kPos_y, kPos_z );
        std::vector<double> contribdatsTmp;
        double quantileSum=0;
        for ( I3VectorI3Particle::const_iterator s = segments->begin();
                s != segments->end(); s++ )
        {
            PhotonicsSource ps(*s);
            double meanPEs, ePointDistance, geoTime;
            photonService->SelectSource ( meanPEs, NULL, ePointDistance, geoTime, ps, true);
            double *quantilesMuon = new double[timeEdges.size()];
            if ( meanPEs > 0 )
            {
                log_debug("mean pes=%f ePointDistance=%f geoTime=%f",meanPEs,ePointDistance, geoTime);
                if ( (*s).GetLength() > 0 ) { meanPEs *= (*s).GetLength(); }
                photonService->GetProbabilityQuantiles( & timeEdges[0],
                        0, quantilesMuon, timeEdges.size() );
            }
            I3Position thisPos( (*s).GetPos().GetX(), (*s).GetPos().GetY(), (*s).GetPos().GetZ() );
            double segdat = ( (*fit).GetPos() - thisPos ).Magnitude();
            log_debug("segdat makeup s=%f,%f,%f fit=%f,%f,%f",(*s).GetPos().GetX(), (*s).GetPos().GetY(), (*s).GetPos().GetZ(),(*fit).GetPos().GetX(),(*fit).GetPos().GetY(),(*fit).GetPos().GetZ() );
            log_debug("segdat1=%f",segdat);
            segdat *= boost::math::sign( thisPos.GetZ() - (*fit).GetPos().GetZ() ) *
                boost::math::sign( (*fit).GetDir().GetZ() );
            //distance metric (building segments for weighted averaging)
            contribdatsTmp.push_back( segdat );
            log_debug("segdat2=%f",segdat);
            //check for table oddities and build the per time yeilds from the segments
            //sometimes the splines behave oddly, you should never see a PE yield below 0 or above 1000
            //if you do something has gone arry and you should just take 0 as the PE yield to be safe
            if ( meanPEs < 0 || meanPEs > 1000 )
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
            for ( unsigned int i = 0; i < timeEdges.size(); i++ )
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
        //Find peak yeild in time so you can use it to set the time window for comparing observed and expected charges
        double maxQuantilePerTime=0;
        for ( unsigned int i = 0; i < timeEdges.size(); i++ )
        {
            if ( sumQuantilesPerTime[i] > maxQuantilePerTime ){ maxQuantilePerTime = sumQuantilesPerTime[i]; }
        }
        //If for some reason you encounter a nan value, indicating something bad, record 0 and keep moving
        if ( std::isnan( quantileSum ) ) { unnormedExQsList.push_back( 0 ); }
        else { unnormedExQsList.push_back( quantileSum ); }
        //distance metric
        if ( quantileSum > 0 )
        {
            contribdats[omkey] = WeightedAvg( contribdatsTmp, sumQuantilesPerSegment );
            //Computing the weighted variance, take the square root for the standard deviation
            contribdatsStdDevs[omkey] = pow( WeightedVar( contribdatsTmp, sumQuantilesPerSegment ), .5 );
            log_debug("Filling contribdats");
            log_debug("dat=%f, std=%f",contribdats[omkey],contribdatsStdDevs[omkey]);
            caddatsStdDevs[omkey] = 0;
            cherdatsStdDevs[omkey] = 0;
        }
        else
        {
            contribdats[omkey] = std::numeric_limits<double>::quiet_NaN();
            contribdatsStdDevs[omkey] = std::numeric_limits<double>::quiet_NaN();
            caddatsStdDevs[omkey] = std::numeric_limits<double>::quiet_NaN();
            cherdatsStdDevs[omkey] = std::numeric_limits<double>::quiet_NaN();
        }
        //Find the observed PEs which match with the expected PE (look for yeild to drop below 1% of peak, thus the .01 below)
        std::vector<unsigned int> gtThreshIdxs;        
        for ( unsigned int i = 0; i < timeEdges.size(); i++ )
        {
            log_debug("sum_quantiles_pertime=%f maxQuantilePerTime=%f",sumQuantilesPerTime[i],maxQuantilePerTime);
            if ( sumQuantilesPerTime[i] > .01 * maxQuantilePerTime ){ gtThreshIdxs.push_back(i); }
        }
        log_debug("index=%d",index);
        double minTime=0;
        double maxTime=0;
        if ( gtThreshIdxs.size() > 0 )
        {
            minTime=timeEdges[gtThreshIdxs[0]];
            maxTime=timeEdges[gtThreshIdxs[gtThreshIdxs.size()-1]+1];
        }
        log_debug("minTime=%f maxTime=%f",minTime,maxTime);
        //If for some reason you encounter a nan value, indicating something bad, record 0 and keep moving
        if ( std::isnan( quantileSum ) ) { expectedBgsList.push_back( 0 ); }
        //The calculation below inserts 100 Hz of noise into every time bin
        else { expectedBgsList.push_back( ( timeEdges[1] - timeEdges[0] ) * 100 * pow( 10, -9 ) ); }
        double coincObsQsTmp=0;
        log_debug("q size=%lu",qs.size());
        for ( unsigned int  i = 0; i< qs.size(); i++ )
        {
            log_debug("t=%f q=%f",ts[i],qs[i]);
            if ( ts[i] > minTime && ts[i] < maxTime ){ coincObsQsTmp+=qs[i]; }
        }
        log_debug("coinc obs q= %f",coincObsQsTmp);
        //If for some reason you encounter a nan value, indicating something bad, record 0 and keep moving
        if ( std::isnan( quantileSum ) ) { coincObsQsList.push_back( 0 ); }
        else { coincObsQsList.push_back( coincObsQsTmp ); }
        index++;
    }
    //Book calculation pieces to the frame
    std::string namePrefix (pulsesName+"_"+fitName);
    std::string nameSufixes[14] = { "allObsQs",
        "coincObsQsList", "unnormedExQsList","exbgs_list",
        "cads", "keyToIndexMap", "caddats", "caddatsStdDevs", "cherdats", 
        "cherdatsStdDevs", "contribdats", "contribdatsStdDevs" };
    log_debug("Names set");
    std::string sls = boost::lexical_cast<std::string>( segments->size() );
    std::vector<std::string> names;
    for ( unsigned int i = 0; i < 14; i++)
    {
        std::string n = namePrefix + "_" + nameSufixes[i] + "_" +sls;
        if ( frame->Has(n) ){ frame->Delete(n); }
        names.push_back(n);
        log_debug("%s",n.c_str());
    }
    log_debug("Begin booking");
    int nameIdx=0;
    frame -> Put( names[nameIdx], boost::make_shared<I3MapKeyDouble>(allObsQs) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3VectorDouble>(coincObsQsList) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3VectorDouble>(unnormedExQsList) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3VectorDouble>(expectedBgsList) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3MapKeyDouble>(cads) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3MapKeyUInt>(keyToIndexMap) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3MapKeyDouble>(caddats) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3MapKeyDouble>(caddatsStdDevs) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3MapKeyDouble>(cherdats) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3MapKeyDouble>(cherdatsStdDevs) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3MapKeyDouble>(contribdats) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3MapKeyDouble>(contribdatsStdDevs) ); nameIdx++;
}


double WeightedAvg( std::vector<double> & xs, I3VectorDouble & ws )
{
    if ( xs.size() == ws.size() )
    {
        //Taken from https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Mathematical_definition
        //Designed to return the weighted average of a set of numbers which is needed by the WeightedVar
        double num = 0;
        double denom = 0;
        for ( unsigned int i = 0; i < xs.size(); i++ )
        {
            log_debug("w_avg i=%ul ws=%f xs=%f",i,ws[i],xs[i]);
            num+= ws[i] * xs[i];
            denom+= ws[i];
            log_debug("w_avg num=%f denom=%f",num,denom);
        }
        log_debug("w_avg num=%f denom=%f",num,denom);
        return num/denom;
    }
    else { return std::numeric_limits<double>::quiet_NaN(); }
}
double WeightedVar( std::vector<double> & xs, I3VectorDouble & ws )
{
    if ( xs.size() == ws.size() )
    {
        //Taken from https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance
        //Designed to return the weighted variance of a set of numbers (xs) given the weights (ws)
        double muStar = WeightedAvg( xs, ws );
        double denom = 0;
        double num = 0;
        unsigned int i;
        for ( i = 0; i < xs.size(); i++)
        {
            log_debug("w_var i=%ul ws=%f xs=%f",i,ws[i],xs[i]);
            num+= ws[i] * pow( xs[i] - muStar, 2 );
            denom+= ws[i];
            log_debug("w_var num=%f denom=%f muStar=%f",num,denom,muStar);
        }
        log_debug("w_var num=%f denom=%f muStar=%f",num,denom,muStar);
        return num/denom;
    }
    else { return std::numeric_limits<double>::quiet_NaN(); }
}

void FindNormalization( I3FramePtr frame, std::string pulsesName,
        std::string fitName, int ls, double tol,
        int iterations, bool supress_stoch, bool cascade )
{
    log_debug("Setting up");
    std::string sls = boost::lexical_cast<std::string>(ls);
    I3VectorDoubleConstPtr obs = frame->Get<I3VectorDoubleConstPtr>(pulsesName+"_"+fitName+"_coincObsQsList_"+sls);
    I3VectorDoubleConstPtr exp = frame->Get<I3VectorDoubleConstPtr>(pulsesName+"_"+fitName+"_unnormedExQsList_"+sls);
    I3VectorDoubleConstPtr bg = frame->Get<I3VectorDoubleConstPtr>(pulsesName+"_"+fitName+"_exbgs_list_"+sls);
    I3MapKeyDoubleConstPtr cherdats = frame->Get<I3MapKeyDoubleConstPtr>(pulsesName+"_"+fitName+"_cherdats_"+sls);
    I3MapKeyDoubleConstPtr cads = frame->Get<I3MapKeyDoubleConstPtr>(pulsesName+"_"+fitName+"_cads_"+sls);
    I3MapKeyUIntConstPtr keyToIndexMap = frame->Get<I3MapKeyUIntConstPtr>(pulsesName+"_"+fitName+"_keyToIndexMap_"+sls);
    std::vector<double> theseCherdats;
    std::vector<OMKey> rejectKeys;
    std::vector<OMKey> rejectKeysCheck;
    std::vector<OMKey> okKeysCheck;
    std::vector<unsigned int> okKeyIdxs;
    double a = 0;
    double llh = std::numeric_limits<double>::quiet_NaN();
    I3VectorDouble normedExQsList;
    if ( (!obs) || (!exp) || (!bg) || (!cherdats) || (!cads) )
    {
        log_warn("Failed to get something out of the frame in FindNormalization");
    }
    else
    {
        log_debug("Starting");
        double minCherdat=pow(10,10);
        //Find the miniumum cherenkov distance where we observed charge,
        //without this we can't do the calculation
        for ( I3MapKeyUInt::const_iterator k = keyToIndexMap->begin();
                k != keyToIndexMap->end(); k++ )
        {
            const OMKey& omkey = k->first;
            if ( (*cherdats).find(omkey)->second <= minCherdat && (*obs)[(*k).second] > 0 )
            {
                minCherdat= (*cherdats).find(omkey)->second;
            }
        }
        log_debug("minCherdat=%f",minCherdat);
        //Stochastics can artificially increase the yeild, we only want this if they are at the beggining
        //because that implies a hadronic interaction from a neutrino
        if ( supress_stoch )
        {
            double obsSum=0;
            int obs_cnt=0;
            for ( I3MapKeyUInt::const_iterator k = keyToIndexMap->begin();
                    k != keyToIndexMap->end(); k++ )
            {
                const OMKey& omkey = k->first;
                if ( (*cherdats).find(omkey)->second >= minCherdat )
                {
                    obsSum += (*obs)[(*k).second];
                    obs_cnt +=1;
                }
            }
            double obsAvg=obsSum/obs_cnt;
            log_debug ("obsAvg=%f",obsAvg);
            //if the amount of charge observed is too close, too far, or brighter than average we are probably looking at the effect from a stochastic
            //To check for this we do the following:
            //1) If a DOM observes light is far from the track where it is not expected this is likely due to a stochastic loss.
            //   Remove for DOMs where less than .1 PE is expected but a photon is observed
            //2) If a DOM expects to observe more than 1 PE of light from a mimimum ionizing muon, it is likely very close to the track.
            //   Remove DOMs where the expected PEs are larger than 1
            //3) If a DOM observes much more light than the average of all DOM observations, then the excess is likely due to a stochastic loss
            //   Remove DOMs where more than 3 times the average(rounded up) PEs are observed.
            //we should check that we have enough good keys so store everything temporarily first
            for ( I3MapKeyUInt::const_iterator k = keyToIndexMap->begin();
                    k != keyToIndexMap->end(); k++ )
            {
                const OMKey& omkey = k->first;
                log_debug("idx=%u obs = %f",(*k).second,(*obs)[(*k).second]);
                log_debug("idx=%u exp = %f",(*k).second,(*exp)[(*k).second]);
                if ( (*cherdats).find(omkey)->second <= minCherdat ) { rejectKeysCheck.push_back(omkey); }
                else if ( (*exp)[(*k).second] <= .1 && (*obs)[(*k).second] > 0 ){ rejectKeysCheck.push_back(omkey); }
                else if ( (*exp)[(*k).second] >= 1 ){ rejectKeysCheck.push_back(omkey); }
                else if ( (*obs)[(*k).second] > 3 * ceil(obsAvg) ){ rejectKeysCheck.push_back(omkey); }
                else{ okKeysCheck.push_back(omkey); }
            }
            //if there are enough keys that check out then we will mask the rejected keys, otherwise do the normal calc to prevent failure
            if ( okKeysCheck.size() > 10 )
            {
                for ( std::vector<OMKey>::const_iterator rkc = rejectKeysCheck.begin();
                        rkc != rejectKeysCheck.end(); rkc++ )
                {
                    rejectKeys.push_back(*rkc);
                }
            }
        }
        //skip if we don't have anything to calculate with
        //otherwise build the list of keys you want to calculate with
        if ( minCherdat != pow(10,10) )
        {
            for ( I3MapKeyUInt::const_iterator k = keyToIndexMap->begin();
                    k != keyToIndexMap->end(); k++ )
            {
                const OMKey& omkey = k->first;
                bool keyOK=true;
                for ( std::vector<OMKey>::const_iterator rk = rejectKeys.begin();
                        rk != rejectKeys.end(); rk++ )
                {
                    if ( (*rk) == omkey ) { keyOK=false; break; }//|| omkey.GetString() >= 79 
                }
                if (cascade){
                    if ( keyOK && (*obs)[(*k).second] > 0)
                    {
                    if ( !supress_stoch && (*cherdats).find(omkey)->second > minCherdat )
                    	{ okKeyIdxs.push_back((*k).second); }
                    if ( supress_stoch )
                        { okKeyIdxs.push_back((*k).second); }
                	}
                	else if ( supress_stoch && (*cherdats).find(omkey)->second > minCherdat && (*cherdats).find(omkey)->second < minCherdat + 250 )
                    	{ okKeyIdxs.push_back((*k).second); }
                }
                else{
                    if ( keyOK )
                    {
                    if ( !supress_stoch && (*cherdats).find(omkey)->second > minCherdat )
                    	{ okKeyIdxs.push_back((*k).second); }
                    if ( supress_stoch )
                        { okKeyIdxs.push_back((*k).second); }
                	}
                	else if ( supress_stoch && (*cherdats).find(omkey)->second > minCherdat && (*cherdats).find(omkey)->second < minCherdat + 250 )
                    	{ okKeyIdxs.push_back((*k).second); }
	            }
            }    
        }
        //This is a hand rolled one variable minimization, probably should replace it with something more official sometime.
        //Start at the normalization value equal to the ratio between the observed and expected charge
        //Probe the llh values some step size (1/1000000th smaller than the tested value) larger and smaller,
        //move to in the direction of the llh improvement
        //The step size constantly decreases by a factor of .9
        //Stop when your llh change is smaller than the tolerance
        log_debug("Normalizing");
        if ( okKeyIdxs.size() > 0 )
        {
            double obsSum=0;
            double expSum=0;
            log_debug("Compute exp and obs sums");
            for ( std::vector<unsigned int>::const_iterator OK_idx = okKeyIdxs.begin();
                    OK_idx != okKeyIdxs.end(); OK_idx++ )
            {
                log_debug("OK_idx=%d obs = %f",(*OK_idx),(*obs)[(*OK_idx)]);
                log_debug("OK_idx=%d exp = %f",(*OK_idx),(*exp)[(*OK_idx)]);
                if ( !std::isnan( (*obs)[(*OK_idx)] ) && !std::isnan( (*exp)[(*OK_idx)] ) )
                {
                    obsSum += (*obs)[(*OK_idx)];
                    expSum += (*exp)[(*OK_idx)];
                }
            }
            log_debug("obs sum=%f expSum=%f",obsSum,expSum);
            if ( obsSum > 0 )
            {
                double aInit = obsSum / expSum;
                double stepSize = aInit / 1000000;
                double llhInit = GetPoissonLlh( aInit, *obs, *exp, *bg );
                log_debug("aInit=%f llhInit=%f",aInit,llhInit);
                a = aInit + stepSize;
                llh = GetPoissonLlh( a, *obs, *exp, *bg );
                log_debug("a=%f llh=%f",a,llh);
                double change = -1 * ( llhInit - llh ) / ( aInit - a );
                double aNext = GetNextA( a, change, stepSize );
                double llhNext = GetPoissonLlh( aNext, *obs, *exp, *bg );
                int iters = 0;
                while ( fabs( llh - llhNext ) > tol )
                {
                    iters += 1;
                    stepSize *= .9;
                    if ( iters > iterations ) { break; }
                    a = aNext;
                    aNext = GetNextA( a, change, stepSize);
                    llh = llhNext;
                    log_debug("a=%f llh=%f",a,llh);
                    llhNext = GetPoissonLlh( aNext, *obs, *exp, *bg );
                }
                llh = llhNext;
                a = aNext;
            }
        }
        //pow(10,10) indicates a failure either to have DOMs to calculate on or 
        //a mismatch in the length of the obs and exp vectors of a DOM, both are
        //a fatal failure and the veto output nothing useful now, and indicate the fit
        //is bad. This is accomplished by setting a to 0. Nothing can come out of the
        //calculations down the line with a=0.
        if ( llh==pow(10,10) )
        {
            a=0;
            log_debug("The StartingTrackVeto cannot compute the normalization, and is failing for this event");
        }
        //Book calculation pieces to the frame
        log_debug("Let's make the normed expectations with a=%f",a);
        for ( I3VectorDouble::const_iterator v = exp->begin();
                v != exp->end(); v++ )
        {
            log_debug("Value was %f now is %f",(*v), (*v) * a);
            normedExQsList.push_back( (*v) * a );
        }
    }
    log_debug("Beginning the booking");
    std::string namePrefix (pulsesName+"_"+fitName);
    std::string nameSufixes[4] = { "normedExQsList", "norm_val", "norm_llh", "norm_llh_ov_nobs" };
    std::vector<std::string> names;
    for ( unsigned int i = 0; i < 4; i++)
    {
        std::string n = namePrefix + "_" + nameSufixes[i] + "_" +sls;
        if ( frame->Has(n) ){ frame->Delete(n); }
        names.push_back(n);
    }
    int nameIdx=0;
    frame -> Put( names[nameIdx], boost::make_shared<I3VectorDouble>(normedExQsList) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3Double>(a) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3Double>(llh) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3Double>(llh/(*obs).size()) ); nameIdx++;
    log_debug("Done with the booking");
}

//give the new normalization value with the constraint that it can't be smaller than 0 or 10 times larger than it is now
double GetNextA(double a, double change, double stepSize )
{
    double aNext = a + ( a * change * stepSize );
    if ( aNext < 0 ){ return a/10; }
    else { return std::min( aNext, a*10 ); }
}

//calculate the poisson probability (llh) of the normalized expected hits in comparison to the observed hits
double GetPoissonLlh( double a, const I3VectorDouble & obs, const I3VectorDouble & exp, 
        const I3VectorDouble & bg )
{
    double llh=0;
    if ( obs.size() == exp.size() )
    {
        for( unsigned int i = 0; i < obs.size(); i++ )
        {
            double lambda = ( exp[i] * a ) + bg[i];
            if ( lambda > 0 )
            {
                boost::math::poisson_distribution<> p( lambda );
                double pTmp=pdf( p, round( obs[i] ) );
                log_debug("exp=%f a=%f bg=%f obs=%f pTmp=%f",exp[i], a , bg[i], obs[i], pTmp);
                llh += log( pTmp );
            }
        }
    }
    else
    {
        llh=pow(10,10);
    }
    return llh;
}

void FindDistances( I3FramePtr frame, std::string pulsesName,
        std::string fitName, int ls )
{
    std::string sls = boost::lexical_cast<std::string>(ls);
    I3MapKeyDoubleConstPtr caddats = frame->Get<I3MapKeyDoubleConstPtr>(pulsesName+"_"+fitName+"_caddats_"+sls);
    I3MapKeyDoubleConstPtr cherdats = frame->Get<I3MapKeyDoubleConstPtr>(pulsesName+"_"+fitName+"_cherdats_"+sls);
    I3MapKeyDoubleConstPtr contribdats = frame->Get<I3MapKeyDoubleConstPtr>(pulsesName+"_"+fitName+"_contribdats_"+sls);
    I3MapKeyDoubleConstPtr caddatsStdDevs = frame->Get<I3MapKeyDoubleConstPtr>(pulsesName+"_"+fitName+"_caddatsStdDevs_"+sls);
    I3MapKeyDoubleConstPtr cherdatsStdDevs = frame->Get<I3MapKeyDoubleConstPtr>(pulsesName+"_"+fitName+"_cherdatsStdDevs_"+sls);
    I3MapKeyDoubleConstPtr contribdatsStdDevs = frame->Get<I3MapKeyDoubleConstPtr>(pulsesName+"_"+fitName+"_contribdatsStdDevs_"+sls);
    I3MapKeyUIntConstPtr keyToIndexMap = frame->Get<I3MapKeyUIntConstPtr>(pulsesName+"_"+fitName+"_keyToIndexMap_"+sls);
    I3VectorDoubleConstPtr normedExp = frame->Get<I3VectorDoubleConstPtr>(pulsesName+"_"+fitName+"_normedExQsList_"+sls);
    I3VectorDoubleConstPtr bg = frame->Get<I3VectorDoubleConstPtr>(pulsesName+"_"+fitName+"_exbgs_list_"+sls);
    I3VectorDoubleConstPtr obs = frame->Get<I3VectorDoubleConstPtr>(pulsesName+"_"+fitName+"_coincObsQsList_"+sls);
    double mincontribdat = pow( 10, 10 );
    I3VectorOMKey mincontribdat_k (1);
    double mincontribdatStdDev = 0;
    double maxcontribdat = -1*pow( 10, 10 );
    I3VectorOMKey maxcontribdatK (1);
    double maxcontribdatStdDev = 0;
    double mincaddat = pow( 10, 10 );
    I3VectorOMKey mincaddatK (1);
    double mincaddatStdDev = 0;
    double maxcaddat = -1*pow( 10, 10 );
    I3VectorOMKey maxcaddatK (1);
    double maxcaddatStdDev = 0;
    double mincherdat = pow( 10, 10 );
    I3VectorOMKey mincherdatK (1);
    double mincherdatStdDev = 0;
    double maxcherdat = -1*pow( 10, 10 );
    I3VectorOMKey maxcherdatK (1);
    double maxcherdatStdDev = 0;
    if ( (!obs) || (!normedExp) || (!bg) || (!cherdats) || (!caddats) || (!contribdats) || (!caddatsStdDevs) || (!cherdatsStdDevs) || (!contribdatsStdDevs) || (!keyToIndexMap) )
    {
        log_warn("Failed to get something out of the frame in FindDistances");
    }
    else
    {
    	//Using the distance matrix maps, find the miniumum and maximum distance along the track of the observed hits
    	//Also find the keys these values originate from, and make an estimate of their uncertainty 
    	for ( I3MapKeyUInt::const_iterator k = keyToIndexMap->begin();
            	k != keyToIndexMap->end(); k++ )
    	{
        	unsigned int index = (*k).second;
            const OMKey& omkey = k->first;
        	if ( (*normedExp)[index] > .01 and (*obs)[index] > 0 )
        	{
            	if ( (*caddats).find(omkey)->second < mincaddat )
            	{
                	mincaddat = (*caddats).find(omkey)->second;
                	mincaddatK[0] = omkey;
                	mincaddatStdDev = (*caddatsStdDevs).find(omkey)->second;
            	}
            	if ( (*caddats).find(omkey)->second > maxcaddat )
            	{
                	maxcaddat = (*caddats).find(omkey)->second;
                	maxcaddatK[0] = omkey;
                	maxcaddatStdDev = (*caddatsStdDevs).find(omkey)->second;
            	}
            	if ( (*cherdats).find(omkey)->second < mincherdat )
            	{
                	mincherdat = (*cherdats).find(omkey)->second;
                	mincherdatK[0] = omkey ;
                	mincherdatStdDev = (*cherdatsStdDevs).find(omkey)->second;
            	}
            	if ( (*cherdats).find(omkey)->second > maxcherdat )
            	{
                	maxcherdat = (*cherdats).find(omkey)->second;
                	maxcherdatK[0] = omkey;
                	maxcherdatStdDev = (*cherdatsStdDevs).find(omkey)->second;
            	}
            	if ( (*contribdats).find(omkey)->second < mincontribdat )
            	{
                	mincontribdat = (*contribdats).find(omkey)->second;
                	mincontribdat_k[0] = omkey;
                	mincontribdatStdDev = (*contribdatsStdDevs).find(omkey)->second;
            	}
            	if ( (*contribdats).find(omkey)->second > maxcontribdat )
            	{
                	maxcontribdat = (*contribdats).find(omkey)->second;
                	maxcontribdatK[0] = omkey;
                	maxcontribdatStdDev = (*contribdatsStdDevs).find(omkey)->second;
            	}
            	log_debug("min cad dat=%f max cad dat=%f",mincaddat,maxcaddat);
            	log_debug("min cad dat std=%f max cad dat std=%f",mincaddatStdDev,maxcaddatStdDev);
            	log_debug("min cher dat=%f max cher dat=%f",mincherdat,maxcherdat);
            	log_debug("min cher dat std=%f max cher dat std=%f",mincherdatStdDev,maxcherdatStdDev);
            	log_debug("min contrib dat=%f max contrib dat=%f",mincontribdat,maxcontribdat);
            	log_debug("min contrib dat std=%f max contrib dat std=%f",mincontribdatStdDev,maxcontribdatStdDev);
        	}
    	}
    }
    log_debug("Final distance results");
    log_debug("min cad dat=%f max cad dat=%f",mincaddat,maxcaddat);
    log_debug("min cad dat std=%f max cad dat std=%f",mincaddatStdDev,maxcaddatStdDev);
    log_debug("min cher dat=%f max cher dat=%f",mincherdat,maxcherdat);
    log_debug("min cher dat std=%f max cher dat std=%f",mincherdatStdDev,maxcherdatStdDev);
    log_debug("min contrib dat=%f max contrib dat=%f",mincontribdat,maxcontribdat);
    log_debug("min contrib dat std=%f max contrib dat std=%f",mincontribdatStdDev,maxcontribdatStdDev);
    //Find the distance between the maximum and minimum values
    double cadD = maxcaddat - mincaddat;
    double cherD = maxcherdat - mincherdat;
    double contribD = maxcontribdat - mincontribdat;
    //Book calculation pieces to the frame
    std::string namePrefix (pulsesName+"_"+fitName);
    std::string nameSufixes[21] = { "mincaddat", "mincaddatK", "maxcaddat",
        "maxcaddatK", "cadD", "mincaddatStdDev", "maxcaddatStdDev",
        "mincherdat", "mincherdatK", "maxcherdat",
        "maxcherdatK", "cherD", "mincherdatStdDev", "maxcherdatStdDev",
        "mincontribdat", "mincontribdat_k", "maxcontribdat",
        "maxcontribdatK", "contribD", "mincontribdatStdDev", "maxcontribdatStdDev" };
    std::vector<std::string> names;
    for ( unsigned int i = 0; i < 21; i++)
    {
        std::string n = namePrefix + "_" + nameSufixes[i] + "_" +sls;
        if ( frame->Has(n) ){ frame->Delete(n); }
        names.push_back(n);
    }
    int nameIdx=0;
    frame -> Put( names[nameIdx], boost::make_shared<I3Double>( mincaddat ) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3VectorOMKey>( mincaddatK ) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3Double>( maxcaddat ) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3VectorOMKey>( maxcaddatK ) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3Double>( cadD ) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3Double>( mincaddatStdDev ) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3Double>( maxcaddatStdDev ) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3Double>( mincherdat ) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3VectorOMKey>( mincherdatK ) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3Double>( maxcherdat ) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3VectorOMKey>( maxcherdatK ) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3Double>( cherD ) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3Double>( mincherdatStdDev ) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3Double>( maxcherdatStdDev ) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3Double>( mincontribdat ) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3VectorOMKey>( mincontribdat_k ) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3Double>( maxcontribdat ) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3VectorOMKey>( maxcontribdatK ) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3Double>( contribD ) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3Double>( mincontribdatStdDev ) ); nameIdx++;
    frame -> Put( names[nameIdx], boost::make_shared<I3Double>( maxcontribdatStdDev ) ); nameIdx++;
}

double AssessProbability( I3FramePtr frame, std::string pulsesName,
        std::string fitName, int ls, std::string dat_type, double sds )
{
    log_debug("Getting things from the frame for a prob calc");
    log_debug("Using dat_type=%s",dat_type.c_str());
    std::string sls = boost::lexical_cast<std::string>(ls);
    I3MapKeyUIntConstPtr keyToIndexMap = frame->Get<I3MapKeyUIntConstPtr>(pulsesName+"_"+fitName+"_keyToIndexMap_"+sls);
    I3VectorDoubleConstPtr normedExp = frame->Get<I3VectorDoubleConstPtr>(pulsesName+"_"+fitName+"_normedExQsList_"+sls);
    I3VectorDoubleConstPtr bg = frame->Get<I3VectorDoubleConstPtr>(pulsesName+"_"+fitName+"_exbgs_list_"+sls);
    I3VectorDoubleConstPtr obs = frame->Get<I3VectorDoubleConstPtr>(pulsesName+"_"+fitName+"_coincObsQsList_"+sls);
    I3MapKeyDoubleConstPtr dats = frame->Get<I3MapKeyDoubleConstPtr>(pulsesName+"_"+fitName+"_"+dat_type+"s_"+sls);
    I3DoubleConstPtr mindat = frame->Get<I3DoubleConstPtr>(pulsesName+"_"+fitName+"_min"+dat_type+"_"+sls);
    I3MapKeyDoubleConstPtr datStdDev = frame->Get<I3MapKeyDoubleConstPtr>(pulsesName+"_"+fitName+"_"+dat_type+"sStdDevs_"+sls);
    I3DoubleConstPtr mindatStdDev = frame->Get<I3DoubleConstPtr>(pulsesName+"_"+fitName+"_min"+dat_type+"StdDev_"+sls);
    log_debug( "mindat=%f mindatStdDev=%f",mindat->value, mindatStdDev->value);
    log_debug("OK got it");
    double cumPMiss=1;
    if ( (!obs) || (!normedExp) || (!bg) || (!dats) || (!mindat) || (!datStdDev) || (!mindatStdDev) || (!keyToIndexMap) )
    {
        log_warn("Failed to get something out of the frame in AssessProbability");
    }
    else
    {
    	//Find the probability of not observing any hits on DOMs before the earliest hit.
    	//If a standard deviation was estimated don't include unhit DOMS within a distance  
    	//(sds) times the earliest hits standard deviation of the earliest hit
    	for ( I3MapKeyUInt::const_iterator k = keyToIndexMap->begin();
            	k != keyToIndexMap->end(); k++ )
    	{
            const OMKey& omkey = k->first;
        	unsigned int index = (*k).second;
        	log_debug("index= %u",index);
        	log_debug("omkey= %d,%d",omkey.GetString(),omkey.GetOM());
        	log_debug("dat=%f",(*dats).find(omkey)->second);
        	log_debug("dat std=%f",(*datStdDev).find(omkey)->second);
        	double lhs = (*dats).find(omkey)->second + ( (*datStdDev).find(omkey)->second * sds );
        	log_debug("%f < ",lhs);
        	double rhs = mindat->value - ( mindatStdDev->value * sds );
        	log_debug("%f ?",rhs);
        	if ( lhs < rhs )
        	{
            	double lambda=( (*normedExp)[index] ) + (*bg)[index];
            	if ( (*obs)[index] == 0 and (*normedExp)[index] != 0  and !std::isnan(lambda) )
            	{
                	boost::math::poisson_distribution<> p( lambda );
                	double pMiss = pdf( p, round( (*obs)[index] ) );
                	//cap the miss probability at 10**-3 to mitigate the effect of tracks which artificially pass too close to a DOM
                	if ( pMiss < pow( 10, -3 ) ){ pMiss = pow( 10, -3 ); }
                	cumPMiss *= pMiss;
                	log_debug("pmiss=%f cum_pmis=%f",pMiss,cumPMiss);
            	}
        	}
    	}
    }
    //Book calculation pieces to the frame
    log_debug("Setting up booking");
    std::string namePrefix (pulsesName+"_"+fitName);
    std::string nameSufixes[1] = { "prob_obs_0s" };
    std::vector<std::string> names;
    for ( unsigned int i = 0; i < 1; i++)
    {
        std::string n = namePrefix + "_" + nameSufixes[i] + "_" + dat_type + "_" +sls;
        if ( frame->Has(n) ){ frame->Delete(n); }
        names.push_back(n);
    }
    log_debug("Booking");
    int nameIdx=0;
    frame -> Put( names[nameIdx], boost::make_shared<I3Double>( cumPMiss ) ); nameIdx++;
    log_debug("Done Booking");
    return cumPMiss;
}
} // namespace StartingTrackVetoUtils
