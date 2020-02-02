use std::boxed::Box;

use slotmap::DenseSlotMap;

use crate::particle::ParticleMap;
use crate::particle_contact::{
    ParticleContact,
    ParticleContactGenerator,
    ParticleContactResolver,
};
use crate::particle_force_generator::ParticleForceRegistry;
use crate::vec3::vec3;

pub struct ParticleWorld {
    pub particles: ParticleMap,
    pub contact_generators: Vec<Box<dyn ParticleContactGenerator>>,
    pub registry: ParticleForceRegistry,
    pub resolver: ParticleContactResolver,
    pub contacts: Vec<ParticleContact>,
    pub max_contacts: usize,
}

impl ParticleWorld {
    pub fn new(max_contacts: usize, iterations: Option<usize>) -> Self {
        ParticleWorld {
            particles: DenseSlotMap::with_key(),
            contact_generators: Vec::new(),
            registry: Vec::new(),
            resolver: ParticleContactResolver::new(
                iterations.unwrap_or(max_contacts * 2)
            ),
            contacts: Vec::with_capacity(max_contacts),
            max_contacts
        }
    }

    pub fn start_frame(&mut self) {
        for particle in self.particles.values_mut() {
            particle.accumulated_force = vec3(0., 0., 0.);
        }
    }

    pub fn generate_contacts(&mut self) -> usize {
        self.contacts.clear();
        for gen in &self.contact_generators {
            if let Some(contact) = gen.add_contact(&self.particles) {
                self.contacts.push(contact);
                if self.contacts.len() == self.contacts.capacity() {
                    break;
                }
            }
        }
        self.contacts.len()
    }

    pub fn integrate(&mut self, duration: f32) {
        for particle in self.particles.values_mut() {
            particle.integrate(duration);
        }
    }

    pub fn run_physics(&mut self, duration: f32) {
        for (particle_key, gen) in &self.registry {
            gen.update_force(*particle_key, duration, &mut self.particles);
        }
        self.integrate(duration);
        let num_contacts = self.generate_contacts();

        if num_contacts > 0 {
            self.resolver.resolve_contacts(
                &self.contacts,
                duration,
                &mut self.particles
            );
        }
    }
}
