#include "PhysicsTools/NanoAOD/interface/SimpleFlatTableProducer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"

typedef SimpleFlatTableProducer<TrackingParticle> SimpleTrackingParticleFlatTableProducer;
typedef SimpleFlatTableProducer<CaloParticle> SimpleCaloParticleFlatTableProducer;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SimpleTrackingParticleFlatTableProducer);
DEFINE_FWK_MODULE(SimpleCaloParticleFlatTableProducer);