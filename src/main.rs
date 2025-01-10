use ggez::{Context, GameResult, graphics::{self, Color}, event};
use mint::Point2;

mod simulation;
use simulation::{Simulation, WINDOW_WIDTH, WINDOW_HEIGHT};

const BASE_RADIUS: f32 = 3.0;

struct MainState {
    simulation: Simulation,
}

impl MainState {
    fn new() -> GameResult<MainState> {
        Ok(MainState {
            simulation: Simulation::new(),
        })
    }

    fn get_particle_color(&self, charge: f32) -> Color {
        // Find maximum absolute charge for color scaling
        let max_charge = self.simulation.particles.iter()
            .map(|p| p.charge.abs())
            .fold(0.0f32, f32::max)
            .max(1e-6); // Prevent division by zero

        let intensity = charge / max_charge;
        if intensity > 0.0 {
            Color::new(intensity, 0.0, 0.0, 1.0) // Red for positive
        } else {
            Color::new(0.0, 0.0, -intensity, 1.0) // Blue for negative
        }
    }

    fn get_particle_radius(&self, mass: f32) -> f32 {
        BASE_RADIUS * mass.sqrt()
    }
}

impl event::EventHandler for MainState {
    fn update(&mut self, _ctx: &mut Context) -> GameResult {
        let dt = 1.0 / 60.0; // Fixed time step
        self.simulation.update(dt);
        Ok(())
    }

    fn draw(&mut self, ctx: &mut Context) -> GameResult {
        let mut canvas = graphics::Canvas::from_frame(ctx, Color::WHITE);

        for particle in &self.simulation.particles {
            let point = Point2 {
                x: particle.position.x,
                y: particle.position.y,
            };
            
            let circle = graphics::Mesh::new_circle(
                ctx,
                graphics::DrawMode::fill(),
                point,
                self.get_particle_radius(particle.mass),
                0.1,
                self.get_particle_color(particle.charge),
            )?;
            canvas.draw(&circle, graphics::DrawParam::default());
        }

        canvas.finish(ctx)?;
        Ok(())
    }
}

fn main() -> GameResult {
    let cb = ggez::ContextBuilder::new("charged_particles", "cursor")
        .window_setup(ggez::conf::WindowSetup::default().title("Charged Particles Simulation"))
        .window_mode(ggez::conf::WindowMode::default().dimensions(WINDOW_WIDTH, WINDOW_HEIGHT));
    
    let (ctx, event_loop) = cb.build()?;
    let state = MainState::new()?;
    event::run(ctx, event_loop, state)
}
