#include <cassert>
#include "duns/Particle.hpp"

using namespace duns;

void Particle::integrate(real duration)
{
    //Don't integrate solid objects
    if (inverseMass <= 0.0f) return;
    
    assert(duration > 0.0);
    
    position.addScaledVector(velocity, duration);
    
    Vector3 resultingAcc = acceleration;
    
    velocity.addScaledVector(resultingAcc, duration);
    
    velocity *= real_pow(damping, duration);
    
    clearAccumulator();
}
