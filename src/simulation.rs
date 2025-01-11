use nalgebra::Vector2;
use rand::Rng;

pub const WINDOW_WIDTH: f32 = 800.0;
pub const WINDOW_HEIGHT: f32 = 600.0;
pub const COULOMB_CONSTANT: f32 = 1000000.0;
const DAMPING: f32 = 0.00;
pub const PARTICLE_COUNT: usize = 30;
pub const VELOCITY_RANGE: f32 = 100.0;

pub struct Particle {
    pub position: Vector2<f32>,
    pub velocity: Vector2<f32>,
    pub mass: f32,
    pub charge: f32,
}

impl Particle {
    pub fn new(x: f32, y: f32, vx: f32, vy: f32, mass: f32, charge: f32) -> Self {
        Particle {
            position: Vector2::new(x, y),
            velocity: Vector2::new(vx, vy),
            mass,
            charge,
        }
    }

    fn handle_wall_collisions(
        position: &mut Vector2<f32>,
        velocity: &mut Vector2<f32>,
        mass: f32,
    ) {
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

    pub fn update(&mut self, dt: f32, force: Vector2<f32>) {
        // F = ma, so a = F/m
        self.velocity += force * (dt / self.mass);
        self.velocity *= 1.0 - DAMPING;
        self.position += self.velocity * dt;

        Self::handle_wall_collisions(
            &mut self.position,
            &mut self.velocity,
            self.mass,
        );
    }
}

pub struct Simulation {
    pub particles: Vec<Particle>,
}

impl Simulation {
    pub fn new() -> Self {
        let mut rng = rand::thread_rng();
        let mut particles = Vec::new();

        // Create some random particles
        for _ in 0..PARTICLE_COUNT {
            let mass: f32 = rng.gen_range(0.5..4.0);
            let radius = 3.0 * mass.sqrt();
            let x = rng.gen_range(radius..WINDOW_WIDTH - radius);
            let y = rng.gen_range(radius..WINDOW_HEIGHT - radius);
            let vx = rng.gen_range(-VELOCITY_RANGE..VELOCITY_RANGE);
            let vy = rng.gen_range(-VELOCITY_RANGE..VELOCITY_RANGE);
            let charge = rng.gen_range(-2.0..2.0);
            particles.push(Particle::new(x, y, vx, vy, mass, charge));
        }

        Simulation { particles }
    }

    pub fn update(&mut self, dt: f32) {
        let forces = self.calculate_forces();
        for (particle, force) in self.particles.iter_mut().zip(forces.iter()) {
            particle.update(dt, *force);
        }
    }

    fn calculate_forces(&self) -> Vec<Vector2<f32>> {
        let mut forces = vec![Vector2::new(0.0, 0.0); self.particles.len()];
        
        // Calculate forces between all pairs of particles
        for i in 0..self.particles.len() {
            for j in (i + 1)..self.particles.len() {
                let diff = self.particles[j].position - self.particles[i].position;
                let distance = diff.magnitude();
                let min_distance = 3.0 * (self.particles[i].mass.sqrt() + self.particles[j].mass.sqrt());
                
                if distance < min_distance {
                    continue; // Prevent division by zero and extreme forces
                }

                let force_magnitude = COULOMB_CONSTANT * self.particles[i].charge * self.particles[j].charge 
                    / (distance * distance);
                let force = diff.normalize() * force_magnitude;

                forces[i] -= force;
                forces[j] += force;
            }
        }
        
        forces
    }
}  