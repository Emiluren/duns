use ultraviolet::vec::Vec3;

use crate::particle::{ParticleKey, ParticleMap};

#[derive(Debug)]
pub struct ParticleContact {
    pub p1: ParticleKey,
    pub p2_opt: Option<ParticleKey>,
    pub restitution: f32,
    pub contact_normal: Vec3,
    pub penetration: f32,
}

impl ParticleContact {
    pub fn calculate_separating_speed(&self, particles: &ParticleMap) -> f32 {
        let p2_velocity = self.p2_opt
            .map(|p2| particles[p2].velocity)
            .unwrap_or(Vec3::new(0., 0., 0.));
        let relative_velocity = particles[self.p1].velocity - p2_velocity;

        relative_velocity.dot(self.contact_normal)
    }

    pub fn resolve(&self, duration: f32, particles: &mut ParticleMap) {
        self.resolve_velocity(duration, particles);
        self.resolve_interpenetration(duration, particles);
    }

    pub fn resolve_velocity(&self, duration: f32, particles: &mut ParticleMap) {
        let separating_speed = self.calculate_separating_speed(particles);

        if separating_speed > 0. {
            return;
        }

        let acc_caused_velocity = particles[self.p1].acceleration +
            self.p2_opt.map(|p2| particles[p2].acceleration)
            .unwrap_or(Vec3::new(0., 0., 0.));
        let acc_caused_sep_speed = acc_caused_velocity.dot(self.contact_normal) *
            duration;

        let mut new_sep_speed = -separating_speed * self.restitution;
        if acc_caused_sep_speed < 0. {
            new_sep_speed =
                (new_sep_speed + self.restitution * acc_caused_sep_speed).max(0.);
        }

        let delta_speed = new_sep_speed - separating_speed;

        let total_inverse_mass = particles[self.p1].inverse_mass +
            self.p2_opt.map(|p2| particles[p2].inverse_mass).unwrap_or(0.);

        if total_inverse_mass <= 0. {
            return;
        }

        let impulse = delta_speed / total_inverse_mass;
        let impulse_per_imass = self.contact_normal * impulse;

        let p1_inv_mass = particles[self.p1].inverse_mass;
        particles[self.p1].velocity += impulse_per_imass * p1_inv_mass;

        if let Some(part2) = self.p2_opt {
            let p2_inv_mass = particles[part2].inverse_mass;
            particles[part2].velocity += impulse_per_imass * p2_inv_mass;
        }
    }

    pub fn resolve_interpenetration(&self, _duration: f32, particles: &mut ParticleMap) {
        if self.penetration <= 0. {
            return;
        }

        let total_inverse_mass = particles[self.p1].inverse_mass +
            self.p2_opt.map(|p| particles[p].inverse_mass).unwrap_or(0.);

        if total_inverse_mass <= 0. {
            return;
        }

        let move_per_imass =
            self.contact_normal * self.penetration / total_inverse_mass;

        let p1_inv_mass = particles[self.p1].inverse_mass;
        particles[self.p1].position += move_per_imass * p1_inv_mass;

        if let Some(part2) = self.p2_opt {
            let p2_inv_mass = particles[part2].inverse_mass;
            particles[part2].position += move_per_imass * p2_inv_mass;
        }
    }
}

pub struct ParticleContactResolver {
    pub iterations: usize,
    pub iterations_used: usize,
}

impl ParticleContactResolver {
    pub fn new(iterations: usize) -> Self {
        ParticleContactResolver { iterations, iterations_used: 0 }
    }

    pub fn resolve_contacts(
        &mut self,
        contacts: &[ParticleContact],
        duration: f32,
        particles: &mut ParticleMap
    ) {
        self.iterations_used = 0;

        while self.iterations_used < self.iterations {
            let max_closing_opt = contacts.iter()
                .filter(|c| c.calculate_separating_speed(particles) < 0. &&
                        c.penetration > 0.)
                .min_by(|c, c2| c.calculate_separating_speed(particles).partial_cmp(
                    &c2.calculate_separating_speed(particles)
                ).unwrap());

            if let Some(max_closing) = max_closing_opt {
                max_closing.resolve(duration, particles);
            } else {
                break;
            }

            self.iterations_used += 1;
        }
    }
}

pub trait ParticleContactGenerator {
    fn add_contact(
        &self,
        particles: &ParticleMap
    ) -> Option<ParticleContact>;
}

pub struct ParticleLink {
    pub p1: ParticleKey,
    pub p2: ParticleKey,
    pub max_length: f32,
    pub restitution: f32,
}

impl ParticleLink {
    fn current_length(&self, particles: &ParticleMap) -> f32 {
        (particles[self.p1].position - particles[self.p2].position).mag()
    }
}

impl ParticleContactGenerator for ParticleLink {
    fn add_contact(
        &self,
        particles: &ParticleMap
    ) -> Option<ParticleContact> {
        let length = self.current_length(particles);

        if length < self.max_length {
            return None;
        }

        let normal =
            (particles[self.p2].position - particles[self.p1].position).normalized();

        Some(ParticleContact {
            p1: self.p1,
            p2_opt: Some(self.p2),
            restitution: self.restitution,
            contact_normal: normal,
            penetration: length - self.max_length,
        })
    }
}

pub struct ParticleRod {
    pub p1: ParticleKey,
    pub p2: ParticleKey,
    pub length: f32,
}

impl ParticleRod {
    fn current_length(&self, particles: &ParticleMap) -> f32 {
        (particles[self.p1].position - particles[self.p2].position).mag()
    }
}

impl ParticleContactGenerator for ParticleRod {
    fn add_contact(
        &self,
        particles: &ParticleMap
    ) -> Option<ParticleContact> {
        let current_length = self.current_length(particles);

        if current_length == self.length {
            return None;
        }

        let normal =
            (particles[self.p2].position - particles[self.p1].position).normalized();

        let (normal, penetration) = if current_length > self.length {
            (normal, current_length - self.length)
        } else {
            (-normal, self.length - current_length)
        };

        Some(ParticleContact {
            p1: self.p1,
            p2_opt: Some(self.p2),
            restitution: 0.,
            contact_normal: normal,
            penetration: penetration,
        })
    }
}

pub struct SpherePlane {
    pub sphere_key: ParticleKey,
    pub radius: f32,
}

impl ParticleContactGenerator for SpherePlane {
    fn add_contact(
        &self,
        particles: &ParticleMap
    ) -> Option<ParticleContact> {
        let pos = particles[self.sphere_key].position;
        if pos.y > self.radius {
            return None;
        }

        Some(ParticleContact {
            p1: self.sphere_key,
            p2_opt: None,
            restitution: 1.,
            contact_normal: Vec3::new(0., 1., 0.),
            penetration: self.radius - pos.y
        })
    }
}
