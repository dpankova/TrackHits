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

#ifndef STARTINGTRACKVETOUTILS_H_INCLUDED
#define STARTINGTRACKVETOUTILS_H_INCLUDED

#include <dataclasses/geometry/I3Geometry.h>
#include <dataclasses/status/I3DetectorStatus.h>
#include <dataclasses/I3Vector.h>
#include <dataclasses/physics/I3Particle.h>
#include <photonics-service/I3PhotonicsService.h>
#include <photonics-service/I3PhotonicsService.h>
#include <math.h>

namespace StartingTrackVetoUtils {
    /* This is a namespace for functions and variables needed to calculate a Track Veto
     *
     */
    void ObtainEventInfo( I3FramePtr frame, I3VectorI3ParticleConstPtr segments,
            I3PhotonicsServicePtr photonService,
            std::vector<double> & timeEdges, std::string pulsesName,
            std::string fitName, double mincadDist,
            I3GeometryConstPtr geo, I3VectorOMKeyPtr badDOMs );
    double WeightedAvg( std::vector<double> & xs, I3VectorDouble & ws );
    double WeightedVar( std::vector<double> & xs, I3VectorDouble & ws );
    void FindNormalization( I3FramePtr frame, std::string pulsesName,
            std::string fitName, int sls, double tol, int iterations,
            bool supressStoch, bool cascade );
    double GetNextA( double a, double change, double stepSize );
    double GetPoissonLlh( double a, const I3VectorDouble & obs, const I3VectorDouble & exp,
            const I3VectorDouble & bg );
    void FindDistances( I3FramePtr frame, std::string pulsesName,
            std::string fitName, int sls );
    double AssessProbability( I3FramePtr frame, std::string pulsesName,
            std::string fitName, int sls, std::string datType, double sds );
};

#endif

