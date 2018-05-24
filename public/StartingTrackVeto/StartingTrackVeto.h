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


#ifndef STARTINGTRACKVETO_H_INCLUDED
#define STARTINGTRACKVETO_H_INCLUDED

#include <icetray/OMKey.h>
#include <dataclasses/geometry/I3Geometry.h>
#include <dataclasses/status/I3DetectorStatus.h>
#include <dataclasses/physics/I3RecoPulse.h>
#include <dataclasses/I3Double.h>
#include <dataclasses/I3Vector.h>
#include <photonics-service/I3PhotonicsService.h>

#include <boost/foreach.hpp>
#include <algorithm>
#include <vector>
#include <StartingTrackVeto/StartingTrackVetoUtils.h>


/**
 * @class \StartingTrackVeto
 * @brief Estimates an events startingness using photon tables,
 */

class StartingTrackVeto : public I3ConditionalModule {
public:
    StartingTrackVeto(const I3Context &);

    void Configure();
    void Physics(I3FramePtr frame);
    void Geometry(I3FramePtr frame);
    void DetectorStatus(I3FramePtr frame);

private:
    std::string pulsesName_;
    std::string fitName_;
    I3PhotonicsServicePtr photonicsService_;
    std::string particleSegmentsName_;
    double timeEdgeMin_;
    double timeEdgeMax_;
    double timeEdgeNSteps_;
    double minCADDist_;
    double percent_;
    double stdDevSp_;
    bool supressStoch_;
    bool cascade_;
    bool norm_;
    double missProbThresh_;
    std::string datType_;
    std::string geoName_;
    std::string badDOMsName_;
    I3GeometryConstPtr geo_;
    I3VectorOMKey badDOMs_;

    SET_LOGGER ("StartingTrackVeto");
};
#endif // STARTINGTRACKVETO_H_INCLUDED
