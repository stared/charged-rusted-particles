use ggez::{Context, GameResult, graphics::{self, Color}, event};
use mint::Point2;
use nalgebra::Vector2;

mod simulation;
use simulation::{Simulation, WINDOW_WIDTH, WINDOW_HEIGHT};

const BASE_RADIUS: f32 = 3.0;
const TRAIL_LENGTH: usize = 20;
const TRAIL_DECAY: f32 = 0.9;

#[derive(Clone)]
struct ParticleTrail {
    positions: Vec<Vector2<f32>>,
}

impl ParticleTrail {
    fn new() -> Self {
        ParticleTrail {
            positions: Vec::with_capacity(TRAIL_LENGTH),
        }
    }

    fn update(&mut self, new_pos: Vector2<f32>) {
        self.positions.insert(0, new_pos);
        if self.positions.len() > TRAIL_LENGTH {
            self.positions.pop();
        }
    }
}

struct MainState {
    simulation: Simulation,
    trails: Vec<ParticleTrail>,
}

impl MainState {
    fn new() -> GameResult<MainState> {
        let simulation = Simulation::new();
        let trails = vec![ParticleTrail::new(); simulation.particles.len()];
        
        Ok(MainState {
            simulation,
            trails,
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

        // Update trails
        for (i, particle) in self.simulation.particles.iter().enumerate() {
            self.trails[i].update(particle.position);
        }
        
        Ok(())
    }

    fn draw(&mut self, ctx: &mut Context) -> GameResult {
        let mut canvas = graphics::Canvas::from_frame(ctx, Color::BLACK);

        // Draw trails
        for (trail, particle) in self.trails.iter().zip(self.simulation.particles.iter()) {
            let base_color = self.get_particle_color(particle.charge);
            
            // Draw trail segments
            for i in 0..trail.positions.len().saturating_sub(1) {
                let pos1 = &trail.positions[i];
                let pos2 = &trail.positions[i + 1];
                let opacity = TRAIL_DECAY.powi(i as i32);
                let color = Color::new(base_color.r, base_color.g, base_color.b, opacity);
                
                let line = graphics::Mesh::new_line(
                    ctx,
                    &[
                        Point2 { x: pos1.x, y: pos1.y },
                        Point2 { x: pos2.x, y: pos2.y },
                    ],
                    2.0,
                    color,
                )?;
                canvas.draw(&line, graphics::DrawParam::default());
            }
        }

        // Draw particles
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
