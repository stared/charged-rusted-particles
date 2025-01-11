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

        Ok(SimulationWrapper {
            simulation: Simulation::new(),
            context,
        })
    }

    pub fn update(&mut self, dt: f32) {
        self.simulation.update(dt);
    }

    pub fn draw(&self) {
        // Clear canvas
        self.context.clear_rect(0.0, 0.0, WINDOW_WIDTH as f64, WINDOW_HEIGHT as f64);
        
        for particle in &self.simulation.particles {
            self.context.begin_path();
            let radius = 3.0 * particle.mass.sqrt();
            
            // Color based on charge
            let color = if particle.charge > 0.0 {
                format!("rgba(255, 0, 0, {})", (particle.charge / 2.0).abs())
            } else {
                format!("rgba(0, 0, 255, {})", (particle.charge / 2.0).abs())
            };
            
            self.context.set_fill_style(&JsValue::from_str(&color));
            
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
    }
} 