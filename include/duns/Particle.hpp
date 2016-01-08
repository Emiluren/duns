#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "duns/Precision.hpp"
#include "duns/Vector3.hpp"

namespace duns
{
    class Particle
    {
    protected:
        Vector3 position;
        Vector3 velocity;
        Vector3 acceleration;
        
        //Damping to remove numerical instability
        real damping;
        //1 / mass to allow infinite mass
        real inverseMass;
        
    public:
        void integrate(real duration);
        void setMass(const real mass) { inverseMass = 1.0 / mass; };
        void setInverseMass(const real invMass) { inverseMass = invMass; };
        real getKineticEnergy;
    }
}

#endif //PARTICLE_HPP
