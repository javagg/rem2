//! rem-bem — Laplace/Helmholtz Boundary Element Method solver
//!
//! Solves exterior Laplace and Helmholtz problems using the boundary integral
//! equation approach with P0 (constant) and P1 (linear) basis functions.
//!
//! # Formulation
//!
//! For the exterior Laplace problem (electrostatics):
//! ```text
//! ½ φ(r) + ∫_S ∂G/∂n'(r,r') φ(r') dS' = ∫_S G(r,r') σ(r') dS'
//! ```
//! where G = 1/(4πR) is the Laplace Green function,
//! σ = ∂φ/∂n is the normal flux (surface charge / ε₀).
//!
//! # Architecture
//! ```text
//! SurfaceMesh (from rem-mom)
//!     ↓
//! assemble_laplace_bem (V + K matrices)
//!     ↓
//! solve (LU via faer)
//!     ↓
//! postprocess (capacitance, potential)
//! ```

pub mod kernel;
pub mod assemble;
pub mod solve;
pub mod postprocess;
