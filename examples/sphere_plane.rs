use duns::particle::Particle;
use duns::particle::contact;
use duns::particle::world::ParticleWorld;

use sdl2::event::Event;
use sdl2::pixels::Color;
use sdl2::render::Canvas;
use sdl2::video::Window;
use ultraviolet::vec::Vec3;

use std::time::Instant;

fn draw_circle((px, py): (i32, i32), radius: f32, canvas: &mut Canvas<Window>) {
    let circle_coords = |i: usize| -> (i32, i32) {
        let angle = i as f32 / 16. * std::f32::consts::PI * 2.;
        (
            px + (angle.cos() * radius) as i32,
            py + (angle.sin() * radius) as i32
        )
    };

    canvas.set_draw_color(Color::RGB(255, 255, 255));
    for i in 0..16 {
        canvas.draw_line(
            circle_coords(i), circle_coords(i + 1)
        ).unwrap();
    }
}

fn main() {
    let sdl_context = sdl2::init().unwrap();
    let video_subsystem = sdl_context.video().unwrap();

    let window = video_subsystem.window("duns sphere_plane demo", 800, 600)
        .resizable()
        .build()
        .unwrap();

    let mut canvas = window.into_canvas().build().unwrap();

    let bg_color = Color::RGB(0, 50, 80);
    canvas.set_draw_color(bg_color);
    canvas.clear();
    canvas.present();

    let mut world = ParticleWorld::new(100, None);

    let sphere_key = world.particles.insert(Particle::new());
    let mut sphere = world.particles.get_mut(sphere_key).unwrap();
    sphere.position = Vec3::new(0., 5., 0.);
    sphere.acceleration = Vec3::new(0., -10., 0.);
    sphere.set_mass(1.0);

    let sphere_radius = 0.5;
    world.contact_generators.push(Box::new(contact::SpherePlane{
        sphere_key, radius: sphere_radius
    }));

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
                _ => {}
            }
        }

        let (screen_w, screen_h) = canvas.output_size().unwrap();

        let new_now = Instant::now();
        accumulated_time += new_now.duration_since(now).as_secs_f32();
        now = new_now;

        while accumulated_time >= dt {
            accumulated_time -= dt;

            world.start_frame();
            world.run_physics(dt);
        }

        canvas.set_draw_color(bg_color);
        canvas.clear();

        let world_pos = world.particles[sphere_key].position;
        let screen_pos = (
            screen_w as i32 / 2 + (world_pos.x * 20.) as i32,
            screen_h as i32 - (world_pos.y * 20.) as i32
        );
        draw_circle(screen_pos, sphere_radius * 20., &mut canvas);

        canvas.present();
    }
}
