use std::cell::RefCell;

use crate::particle::Particle;
use crate::vec3::{Vec3, vec3};

pub trait ParticleForceGenerator {
    fn update_force(&self, particle: &mut Particle, duration: f32);
}

pub type ParticleForceRegistry<'a> = Vec<(
    &'a RefCell<Particle>,
    &'a RefCell<dyn ParticleForceGenerator>
)>;

pub fn update_forces(particle_force_registry: &ParticleForceRegistry, duration: f32) {
    for (particle_cell, force_generator_cell) in particle_force_registry {
        let mut particle = particle_cell.borrow_mut();
        let force_generator = force_generator_cell.borrow();

        force_generator.update_force(&mut particle, duration);
    }
}

pub struct ParticleGravity {
    pub gravity: Vec3,
}

impl ParticleGravity {
    pub fn new(gravity: Vec3) -> Self {
        ParticleGravity { gravity }
    }
}

impl ParticleForceGenerator for ParticleGravity {
    fn update_force(&self, particle: &mut Particle, _duration: f32) {
        if particle.inverse_mass == 0. {
            return;
        }

        particle.accumulated_force += self.gravity * particle.mass()
    }
}

pub struct ParticleDrag {
    pub k1: f32,
    pub k2: f32,
}

impl ParticleDrag {
    pub fn new(k1: f32, k2: f32) -> Self {
        ParticleDrag { k1, k2 }
    }
}

impl ParticleForceGenerator for ParticleDrag {
    fn update_force(&self, particle: &mut Particle, _duration: f32) {
        let speed = particle.velocity.norm();
        let drag_coeff =
            self.k1 * speed +
            self.k2 * speed * speed;

        let force = particle.velocity / speed * -drag_coeff;
        particle.accumulated_force += force;
    }
}

pub struct ParticleSpring<'a> {
    pub other: &'a RefCell<Particle>,
    pub spring_constant: f32,
    pub rest_length: f32,
}

impl<'a> ParticleSpring<'a> {
    pub fn new(other: &'a RefCell<Particle>, spring_constant: f32, rest_length: f32) -> Self {
        ParticleSpring { other, spring_constant, rest_length }
    }
}

fn spring_force(p1: Vec3, p2: Vec3, rest_length: f32, spring_constant: f32) -> Vec3 {
    let to_other = p1 - p2;
    let dist = to_other.norm();
    let force = to_other / dist *
        -(dist - rest_length).abs() *
        spring_constant;
    force
}

impl ParticleForceGenerator for ParticleSpring<'_> {
    fn update_force(&self, particle: &mut Particle, _duration: f32) {
        particle.accumulated_force += spring_force(
            particle.position,
            self.other.borrow().position,
            self.rest_length,
            self.spring_constant
        );
    }
}

pub struct ParticleAnchoredSpring {
    anchor: Vec3,
    spring_constant: f32,
    rest_length: f32,
}

impl ParticleAnchoredSpring {
    pub fn new(anchor: Vec3, spring_constant: f32, rest_length: f32) -> Self {
        ParticleAnchoredSpring { anchor, spring_constant, rest_length }
    }
}

impl ParticleForceGenerator for ParticleAnchoredSpring {
    fn update_force(&self, particle: &mut Particle, _duration: f32) {
        particle.accumulated_force += spring_force(
            particle.position,
            self.anchor,
            self.rest_length,
            self.spring_constant
        );
    }
}

pub struct ParticleBungee<'a> {
    pub other: &'a RefCell<Particle>,
    pub spring_constant: f32,
    pub rest_length: f32,
}

impl<'a> ParticleBungee<'a> {
    pub fn new(other: &'a RefCell<Particle>, spring_constant: f32, rest_length: f32) -> Self {
        ParticleBungee { other, spring_constant, rest_length }
    }
}

impl ParticleForceGenerator for ParticleBungee<'_> {
    fn update_force(&self, particle: &mut Particle, _duration: f32) {
        let other = self.other.borrow();
        let to_other = particle.position - other.position;
        let dist = to_other.norm();

        if dist <= self.rest_length {
            return;
        }

        let force = to_other / dist *
            (self.rest_length - dist) *
            self.spring_constant;
        particle.accumulated_force += force;
    }
}

pub struct ParticelBuoyancy {
    pub max_depth: f32,
    pub volume: f32,
    pub water_height: f32,
    pub liquid_density: f32, // In kg / m^3
}

impl ParticelBuoyancy {
    pub fn new(max_depth: f32, volume: f32, water_height: f32) -> Self {
        ParticelBuoyancy {
            max_depth,
            volume,
            water_height,
            liquid_density: 1000.
        }
    }
}

impl ParticleForceGenerator for ParticelBuoyancy {
    fn update_force(&self, particle: &mut Particle, _duration: f32) {
        let depth = particle.position.y;
        let out_of_water = depth >= self.water_height + self.max_depth;
        if out_of_water {
            return;
        }

        let at_max_depth = depth <= self.water_height - self.max_depth;
        if at_max_depth {
            particle.accumulated_force += vec3(0., self.liquid_density * self.volume, 0.);
        } else {
            particle.accumulated_force += vec3(
                0.,
                self.liquid_density * self.volume *
                    (depth - self.max_depth - self.water_height) / 2. * self.max_depth,
                0.
            )
        }
    }
}

pub struct ParticleFakeSpring {
    pub anchor: Vec3,
    pub spring_constant: f32,
    pub damping: f32,
}

impl ParticleFakeSpring {
    pub fn new(anchor: Vec3, spring_constant: f32, damping: f32) -> Self {
        ParticleFakeSpring { anchor, spring_constant, damping }
    }
}

impl ParticleForceGenerator for ParticleFakeSpring {
    fn update_force(&self, particle: &mut Particle, duration: f32) {
        if particle.inverse_mass == 0. {
            return;
        }

        let to_anchor = particle.position - self.anchor;
        let gamma = 0.5 * (4. * self.spring_constant - self.damping.powi(2)).sqrt();
        if gamma == 0. {
            return;
        }

        let c = to_anchor * (self.damping / (2. * gamma)) +
            particle.velocity * (1. / gamma);

        let target = (to_anchor * (gamma * duration).cos() +
                      c * (gamma * duration).sin()) *
            (-0.5 * duration * self.damping).exp();

        let accel = (target - to_anchor) * (1. / duration*duration) -
            particle.velocity * duration;
        particle.accumulated_force += accel * particle.mass();
    }
}
