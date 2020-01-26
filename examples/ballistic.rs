use duns::particle::Particle;
use duns::vec3::vec3;

use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use sdl2::pixels::Color;
use sdl2::rect::Rect;
use std::time::{Duration, Instant};
use strum::IntoEnumIterator;
use strum_macros::EnumIter;

pub fn modulo(x: i32, div: i32) -> i32 {
    (x % div + div) % div
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, EnumIter, PartialOrd, Ord)]
enum ShotType {
    Pistol, Artillery, Fireball, Laser
}

#[derive(Clone)]
struct Shot {
    particle: Particle,
    shot_type: ShotType,
    start_time: Instant,
}

fn fire_shot(shot_type: ShotType) -> Shot {
    let shooting_position = vec3(0., 5., 0.);

    let mut particle = Particle::new();
    particle.position = shooting_position;

    match shot_type {
        ShotType::Pistol => {
            particle.set_mass(2.);
            particle.velocity = vec3(35., 0., 0.);
            particle.acceleration = vec3(0., -1., 0.);
        },
        ShotType::Artillery => {
            particle.set_mass(200.);
            particle.velocity = vec3(40., 30., 0.); // 50 m/s
            particle.acceleration = vec3(0., -20., 0.);
        },
        ShotType::Fireball => {
            particle.set_mass(1.);
            particle.velocity = vec3(10., 0., 0.);
            particle.acceleration = vec3(0., 0.6, 0.);
            particle.damping = 0.9;
        },
        ShotType::Laser => {
            particle.set_mass(0.1);
            particle.velocity = vec3(100., 0., 0.);
            particle.acceleration = vec3(0., 0., 0.);
        },
    }

    Shot {
        particle,
        shot_type,
        start_time: Instant::now(),
    }
}

fn main() {
    let sdl_context = sdl2::init().unwrap();
    let video_subsystem = sdl_context.video().unwrap();

    let window = video_subsystem.window("duns ballistic demo", 800, 600)
        .resizable()
        .build()
        .unwrap();

    let mut canvas = window.into_canvas().build().unwrap();
 
    canvas.set_draw_color(Color::RGB(0, 50, 80));
    canvas.clear();
    canvas.present();

    let shot_types: Vec<_> = ShotType::iter().collect();

    let mut shot_index = 0;
    let mut ammo = vec![
        None; 100
    ];

    let dt = 1. / 60.;
    let mut now = Instant::now();
    let mut accumulated_time = 0.;

    let mut event_pump = sdl_context.event_pump().unwrap();
    'running: loop {
        for event in event_pump.poll_iter() {
            match event {
                Event::Quit {..} => {
                    break 'running
                }
                Event::KeyDown { keycode: Some(kc), .. } => {
                    match kc {
                        Keycode::Space => {
                            for shot in &mut ammo {
                                if shot.is_none() {
                                    *shot = Some(fire_shot(shot_types[shot_index]));
                                    break;
                                }
                            }
                        }
                        Keycode::Right => {
                            shot_index = if shot_index == shot_types.len() - 1 {
                                0
                            } else {
                                shot_index + 1
                            };
                            println!(
                                "Changing shot type to {:?}",
                                shot_types[shot_index]
                            );
                        }
                        Keycode::Left => {
                            shot_index = if shot_index == 0 {
                                shot_types.len() - 1
                            } else {
                                shot_index - 1
                            };
                            println!(
                                "Changing shot type to {:?}",
                                shot_types[shot_index]
                            );
                        }
                        _ => {}
                    }
                }
                _ => {}
            }
        }

        let (screen_w, screen_h) = canvas.output_size().unwrap();

        let new_now = Instant::now();
        accumulated_time += new_now.duration_since(now).as_secs_f32();
        now = new_now;

        while accumulated_time >= dt {
            accumulated_time -= dt;

            for opt_shot in &mut ammo {
                let mut remove_shot = false;
                if let Some(shot) = opt_shot {
                    shot.particle.integrate(dt);
                    let pos = shot.particle.position;
                    if pos.y < 0. ||
                        now.duration_since(shot.start_time) > Duration::from_secs(5) ||
                        (pos.x * 10.) as i32 > screen_w as i32 {
                            remove_shot = true;
                    }
                }
                if remove_shot {
                    *opt_shot = None;
                }
            }
        }

        canvas.set_draw_color(Color::RGB(0, 50, 80));
        canvas.clear();

        canvas.set_draw_color(Color::RGB(255, 255, 80));
        for opt_shot in &ammo {
            if let Some(shot) = opt_shot {
                let p = shot.particle.position;
                canvas.fill_rect(Rect::new(
                    (p.x * 10.) as i32,
                    screen_h as i32 - (p.y * 10.) as i32,
                    20,
                    20
                )).unwrap();
            }
        }
        canvas.present();
    }
}
