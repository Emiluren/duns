use crate::vec3::{vec3, Vec3};

#[derive(Clone)]
pub struct Particle {
    pub position: Vec3,
    pub velocity: Vec3,
    pub acceleration: Vec3,
    pub accumulated_force: Vec3,

    pub damping: f32,
    pub inverse_mass: f32,
}

impl Particle {
    pub fn new() -> Self {
        Particle {
            position: vec3(0., 0., 0.),
            velocity: vec3(0., 0., 0.),
            acceleration: vec3(0., 0., 0.),
            accumulated_force: vec3(0., 0., 0.),

            damping: 0.99,
            inverse_mass: 0.,
        }
    }

    pub fn set_mass(&mut self, mass: f32) {
        self.inverse_mass = 1. / mass;
    }

    pub fn integrate(&mut self, duration: f32) {
        self.position += self.velocity * duration;

        let acceleration = self.acceleration +
            self.accumulated_force * self.inverse_mass;

        self.velocity += acceleration * duration;

        self.velocity *= self.damping.powf(duration);
    }
}
