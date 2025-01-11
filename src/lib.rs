use wasm_bindgen::prelude::*;
use web_sys::{CanvasRenderingContext2d, HtmlCanvasElement};

mod simulation;
use simulation::*;

#[wasm_bindgen]
pub struct SimulationWrapper {
    simulation: Simulation,
    context: CanvasRenderingContext2d,
}

#[wasm_bindgen]
impl SimulationWrapper {
    #[wasm_bindgen(constructor)]
    pub fn new(canvas: HtmlCanvasElement) -> Result<SimulationWrapper, JsValue> {
        let context = canvas
            .get_context("2d")?
            .unwrap()
            .dyn_into::<CanvasRenderingContext2d>()?;

        // Set initial background
        context.set_fill_style_str("#000000");
        context.fill_rect(0.0, 0.0, WINDOW_WIDTH as f64, WINDOW_HEIGHT as f64);

        Ok(SimulationWrapper {
            simulation: Simulation::new(),
            context,
        })
    }

    pub fn update(&mut self, dt: f32) {
        self.simulation.update(dt);
    }

    pub fn draw(&self, decay_factor: f32) {
        // Apply trail effect by drawing a transparent black rectangle
        // The opacity is inverse of decay_factor to create a linear fade
        self.context.set_global_alpha(1.0 - decay_factor as f64);
        self.context.set_fill_style_str("#000000");
        self.context.fill_rect(0.0, 0.0, WINDOW_WIDTH as f64, WINDOW_HEIGHT as f64);
        
        // Reset alpha for particle drawing
        self.context.set_global_alpha(1.0);
        
        for particle in &self.simulation.particles {
            self.context.begin_path();
            let radius = 3.0 * particle.mass.sqrt();
            
            // Color based on charge with glow effect
            let opacity = (particle.charge / 2.0).abs();
            let color = if particle.charge > 0.0 {
                "rgb(255, 50, 50)"
            } else {
                "rgb(50, 50, 255)"
            };
            
            self.context.set_global_alpha(opacity as f64);
            self.context.set_fill_style_str(color);
            
            self.context.arc(
                particle.position.x as f64,
                particle.position.y as f64,
                radius as f64,
                0.0,
                2.0 * std::f64::consts::PI,
            ).unwrap();
            
            self.context.fill();
        }
    }

    pub fn set_coulomb_constant(&mut self, value: f32) {
        unsafe {
            COULOMB_CONSTANT = value;
        }
    }

    pub fn set_damping(&mut self, value: f32) {
        unsafe {
            DAMPING = value;
        }
    }

    pub fn set_velocity_range(&mut self, value: f32) {
        unsafe {
            VELOCITY_RANGE = value;
        }
    }

    pub fn reset_with_particle_count(&mut self, count: usize) {
        unsafe {
            PARTICLE_COUNT = count;
        }
        self.simulation = Simulation::new();
        
        // Reset background
        self.context.set_fill_style_str("#000000");
        self.context.fill_rect(0.0, 0.0, WINDOW_WIDTH as f64, WINDOW_HEIGHT as f64);
    }

    pub fn set_decay_factor(&mut self, _value: f32) {
        // No need to store the decay factor as it's passed directly to draw
    }
} 