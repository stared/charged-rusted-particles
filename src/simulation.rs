use nalgebra::Vector2;
use rand::Rng;

pub const WINDOW_WIDTH: f32 = 800.0;
pub const WINDOW_HEIGHT: f32 = 600.0;
pub static mut COULOMB_CONSTANT: f32 = 1_000_000.0;
pub static mut DAMPING: f32 = 0.00;
pub static mut PARTICLE_COUNT: usize = 30;
pub static mut VELOCITY_RANGE: f32 = 2.0;
pub static mut USE_VERLET: bool = true;

pub struct Particle {
    pub position: Vector2<f32>,
    pub velocity: Vector2<f32>,
    pub mass: f32,
    pub charge: f32,
}

impl Particle {
    /// Former `update` method, renamed to `step_euler` for reference.
    /// This uses a simple (forward) Euler integration approach.
    pub fn step_euler(&mut self, dt: f32, force: Vector2<f32>) {
        // F = m*a, so a = F / m
        self.velocity += force * (dt / self.mass);
        unsafe {
            self.velocity *= 1.0 - DAMPING;
        }
        self.position += self.velocity * dt;

        Self::handle_wall_collisions(&mut self.position, &mut self.velocity, self.mass);
    }

    fn handle_wall_collisions(position: &mut Vector2<f32>, velocity: &mut Vector2<f32>, mass: f32) {
        let radius = 3.0 * mass.sqrt();

        // Bounce off walls
        if position.x < radius {
            position.x = radius;
            velocity.x = velocity.x.abs();
        } else if position.x > WINDOW_WIDTH - radius {
            position.x = WINDOW_WIDTH - radius;
            velocity.x = -velocity.x.abs();
        }

        if position.y < radius {
            position.y = radius;
            velocity.y = velocity.y.abs();
        } else if position.y > WINDOW_HEIGHT - radius {
            position.y = WINDOW_HEIGHT - radius;
            velocity.y = -velocity.y.abs();
        }
    }
}

pub struct Simulation {
    pub particles: Vec<Particle>,
}

impl Simulation {
    pub fn new() -> Self {
        let mut rng = rand::thread_rng();
        let mut particles = Vec::new();

        unsafe {
            for _ in 0..PARTICLE_COUNT {
                let mass: f32 = rng.gen_range(0.5..4.0);
                let radius = 3.0 * mass.sqrt();
                let x = rng.gen_range(radius..WINDOW_WIDTH - radius);
                let y = rng.gen_range(radius..WINDOW_HEIGHT - radius);
                let vx = rng.gen_range(-VELOCITY_RANGE..VELOCITY_RANGE);
                let vy = rng.gen_range(-VELOCITY_RANGE..VELOCITY_RANGE);
                let charge = rng.gen_range(-2.0..2.0);

                particles.push(Particle {
                    position: Vector2::new(x, y),
                    velocity: Vector2::new(vx, vy),
                    mass,
                    charge,
                });
            }
        }

        Simulation { particles }
    }

    /// Example of how you might integrate using Euler,
    /// leaving the new Velocity Verlet method for you to adapt:
    pub fn update(&mut self, dt: f32) {
        unsafe {
            if USE_VERLET {
                // First calculate initial forces
                let old_forces = self.calculate_forces();
                
                // Update positions and half-step velocities
                for (particle, force) in self.particles.iter_mut().zip(old_forces.iter()) {
                    // First half of verlet step
                    let half_dt = 0.5 * dt;
                    particle.velocity += *force * (half_dt / particle.mass);
                    particle.position += particle.velocity * dt;
                    Particle::handle_wall_collisions(&mut particle.position, &mut particle.velocity, particle.mass);
                }

                // Calculate new forces at new positions
                let new_forces = self.calculate_forces();

                // Complete velocity updates
                for (particle, new_force) in self.particles.iter_mut().zip(new_forces.iter()) {
                    // Second half of verlet step
                    let half_dt = 0.5 * dt;
                    particle.velocity += *new_force * (half_dt / particle.mass);
                    particle.velocity *= 1.0 - DAMPING;
                }
            } else {
                let forces = self.calculate_forces();
                for (particle, force) in self.particles.iter_mut().zip(forces.iter()) {
                    particle.step_euler(dt, *force);
                }
            }
        }
    }

    fn calculate_forces(&self) -> Vec<Vector2<f32>> {
        let mut forces = vec![Vector2::zeros(); self.particles.len()];

        for i in 0..self.particles.len() {
            for j in (i + 1)..self.particles.len() {
                let diff = self.particles[j].position - self.particles[i].position;
                let distance = diff.magnitude();
                let min_distance =
                    3.0 * (self.particles[i].mass.sqrt() + self.particles[j].mass.sqrt());

                if distance < min_distance {
                    continue; // Skip to avoid extreme forces
                }

                unsafe {
                    let force_magnitude = COULOMB_CONSTANT
                        * self.particles[i].charge
                        * self.particles[j].charge
                        / (distance * distance);

                    let force = diff.normalize() * force_magnitude;
                    forces[i] -= force;
                    forces[j] += force;
                }
            }
        }

        forces
    }
}  