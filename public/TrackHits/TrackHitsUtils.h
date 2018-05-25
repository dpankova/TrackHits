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

#ifndef TRACKHITSUTILS_H_INCLUDED
#define TRACKHITSUTILS_H_INCLUDED

#include <dataclasses/geometry/I3Geometry.h>
#include <dataclasses/status/I3DetectorStatus.h>
#include <dataclasses/I3Vector.h>
#include <dataclasses/physics/I3Particle.h>
#include <photonics-service/I3PhotonicsService.h>
#include <photonics-service/I3PhotonicsService.h>
#include <math.h>

namespace TrackHitsUtils 
{
    /* This is a namespace for functions and variables needed to find TrackHits */
    void ObtainEventInfo(I3FramePtr frame, 
			 I3VectorI3ParticleConstPtr segments,
			 I3PhotonicsServicePtr photonService,
			 std::vector<double> & timeEdges, 
			 std::string pulsesName,
			 std::string fitName, 
			 double mincadDist,
			 double percent,
			 I3GeometryConstPtr geo, 
			 I3VectorOMKeyPtr badDOMs );
    
};

#endif

