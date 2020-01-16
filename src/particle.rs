use crate::vec3::Vec3;

pub struct Particle {
    pub position: Vec3,
    pub velocity: Vec3,
    pub acceleration: Vec3,
    pub accumulated_force: Vec3,

    pub damping: f32,
    pub inverse_mass: f32,
}

impl Particle {
    pub fn integrate(&mut self, duration: f32) {
        self.position += self.velocity * duration;

        let acceleration = self.acceleration +
            self.accumulated_force * self.inverse_mass;

        self.velocity += acceleration * duration;

        self.velocity *= self.damping.powf(duration);
    }
}
