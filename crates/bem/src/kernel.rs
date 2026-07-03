//! BEM kernels — currently only free-space Laplace.

use std::f64::consts::PI;

pub trait BemKernel: Send + Sync {
    fn g(&self, r: &[f64; 3], rp: &[f64; 3]) -> f64;
    fn dg_dn(&self, r: &[f64; 3], rp: &[f64; 3], n_prime: &[f64; 3]) -> f64;
}

pub struct LaplaceKernel;

impl BemKernel for LaplaceKernel {
    fn g(&self, r: &[f64; 3], rp: &[f64; 3]) -> f64 {
        let d2 = (r[0]-rp[0]).powi(2)+(r[1]-rp[1]).powi(2)+(r[2]-rp[2]).powi(2);
        if d2 < 1e-28 { return 0.0; }
        1.0 / (4.0 * PI * d2.sqrt())
    }
    fn dg_dn(&self, r: &[f64; 3], rp: &[f64; 3], n: &[f64; 3]) -> f64 {
        let rx=r[0]-rp[0]; let ry=r[1]-rp[1]; let rz=r[2]-rp[2];
        let d2 = rx*rx+ry*ry+rz*rz;
        if d2 < 1e-28 { return 0.0; }
        -(rx*n[0]+ry*n[1]+rz*n[2]) / (4.0*PI*d2*d2.sqrt())
    }
}

pub fn laplace_G(r: &[f64;3], rp: &[f64;3]) -> f64 { LaplaceKernel.g(r, rp) }
pub fn laplace_dG_dn(r: &[f64;3], rp: &[f64;3], n: &[f64;3]) -> f64 { LaplaceKernel.dg_dn(r, rp, n) }

#[cfg(test)]
mod tests {
    use super::*;
    #[test] fn free_1o4pi() { let r=[1.0,0.0,0.0]; let rp=[0.0;3];
        assert!((laplace_G(&r,&rp)-1.0/(4.0*PI)).abs()<1e-14); }
    #[test] fn symmetric() { let r=[1.0,2.0,3.0]; let rp=[0.5,0.1,0.7];
        assert!((laplace_G(&r,&rp)-laplace_G(&rp,&r)).abs()<1e-14); }
}
