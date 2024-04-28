#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include <iostream>

PCaloHit::PCaloHit(float eEM, float eHad, float t, int i, uint16_t d, math::XYZPoint position) : myTime(t), myItra(i), myDepth(d),myPosition(position) {
  std::cout << "3rd Constructor" << std::endl;

  myEnergy = eEM + eHad;
  myEMFraction = (myEnergy <= 0.f ? 1.f : eEM / myEnergy);
}

PCaloHit::PCaloHit(unsigned int id, float eEM, float eHad, float t, int i, uint16_t d, math::XYZPoint position)
//PCaloHit::PCaloHit(unsigned int id, float eEM, float eHad, float t, int i, uint16_t d)
    : myTime(t), myItra(i), detId(id), myDepth(d),myPosition(position) {
  std::cout << "4th Constructor" << std::endl;

  myEnergy = eEM + eHad;
  myEMFraction = (myEnergy <= 0.f ? 1.f : eEM / myEnergy);
}

std::ostream& operator<<(std::ostream& o, const PCaloHit& hit) {
  o << "0x" << std::hex << hit.id() << std::dec << ": Energy (EM) " << hit.energyEM() << " GeV "
    << ": Energy (Had) " << hit.energyHad() << " GeV "
    << " Tof " << hit.time() << " ns "
    << " Geant track #" << hit.geantTrackId() << " Encoded depth " << hit.depth();

  return o;
}
