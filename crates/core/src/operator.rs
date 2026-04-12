//! Unified linear operator interface supporting various matrix types.
//!
//! This module provides a trait-based abstraction for matrix-vector operations,
//! enabling uniform access to dense matrices, sparse matrices, and matrix-free operators.
//!
//! # Design
//!
//! - **LinearOperator<T>**: Trait for y ← A*x with type parameter T (f64, Complex64)
//! - **Adapters**: Implementations for DMatrix, CsrMatrix, etc.
//! - **Solvers**: Generic GMRES/BiCGSTAB using LinearOperator
//!
//! # Example
//!
//! ```ignore
//! use nalgebra::{DMatrix, DVector};
//! use num_complex::Complex64;
//! use rem_core::LinearOperator;
//!
//! let z: DMatrix<Complex64> = /* ... */;
//! let b = DVector::zeros(z.ncols());
//!
//! // z implements LinearOperator<Complex64> automatically
//! let (m, n) = z.size();
//! let mut y = DVector::zeros(m);
//! z.matvec(&b, &mut y)?;
//! ```

use nalgebra::{ComplexField, DMatrix, DVector};
use num_complex::Complex64;

/// Trait for matrix-vector operations with type-based genericity.
///
/// Implementations must satisfy:
/// - `matvec(x, y)` computes y ← A*x (may accumulate if y is nonzero)
/// - `matvec_adjoint(x, y)` computes y ← A^H*x (for iterative solvers)
/// - `diagonal()` returns diag(A) if available (for preconditioning)
///
/// # Safety
///
/// All methods must respect vector dimensions (nrows, ncols).
/// Implementations should panic or return Err on dimension mismatch.
pub trait LinearOperator<T: ComplexField>: Send + Sync {
    /// Return (nrows, ncols) of the operator A.
    fn size(&self) -> (usize, usize);

    /// Return (ncols, nrows) - transposed dimensions for adjoint operations.
    fn size_adjoint(&self) -> (usize, usize) {
        let (m, n) = self.size();
        (n, m)
    }

    /// Compute y ← A*x.
    ///
    /// # Errors
    ///
    /// Returns Err if dimensions don't match or operation fails.
    fn matvec(&self, x: &DVector<T>, y: &mut DVector<T>) -> Result<(), String>;

    /// Compute y ← A^H*x (adjoint/conjugate transpose for Complex64).
    ///
    /// Default implementation returns Err. Solvers should provide fallback behavior.
    fn matvec_adjoint(&self, _x: &DVector<T>, _y: &mut DVector<T>) -> Result<(), String> {
        Err("matvec_adjoint not implemented".to_string())
    }

    /// Extract main diagonal diag(A).
    ///
    /// Default returns None. Useful for preconditioning.
    fn diagonal(&self) -> Option<DVector<T>> {
        None
    }

    /// Return matrix density estimate (0.0 = sparse, 1.0 = dense).
    /// Used by adaptive solver selection.
    fn density(&self) -> f64 {
        let (m, n) = self.size();
        if m * n == 0 {
            0.0
        } else {
            1.0  // Conservative default
        }
    }
}

// ---------------------------------------------------------------------------
// Adapters for common matrix types
// ---------------------------------------------------------------------------

/// Adapter for dense real matrices.
impl LinearOperator<f64> for DMatrix<f64> {
    fn size(&self) -> (usize, usize) {
        (self.nrows(), self.ncols())
    }

    fn matvec(&self, x: &DVector<f64>, y: &mut DVector<f64>) -> Result<(), String> {
        let (m, n) = self.size();
        if x.len() != n || y.len() != m {
            return Err(format!(
                "matvec dimension mismatch: matrix {}×{}, x len {}, y len {}",
                m, n, x.len(), y.len()
            ));
        }
        *y = self * x;
        Ok(())
    }

    fn matvec_adjoint(&self, x: &DVector<f64>, y: &mut DVector<f64>) -> Result<(), String> {
        let (m, n) = self.size();
        if x.len() != m || y.len() != n {
            return Err(format!(
                "matvec_adjoint dimension mismatch: matrix {}×{}, x len {}, y len {}",
                m, n, x.len(), y.len()
            ));
        }
        *y = self.transpose() * x;
        Ok(())
    }

    fn diagonal(&self) -> Option<DVector<f64>> {
        Some(self.diagonal())
    }

    fn density(&self) -> f64 {
        1.0  // Dense
    }
}

/// Adapter for dense complex matrices.
impl LinearOperator<Complex64> for DMatrix<Complex64> {
    fn size(&self) -> (usize, usize) {
        (self.nrows(), self.ncols())
    }

    fn matvec(&self, x: &DVector<Complex64>, y: &mut DVector<Complex64>) -> Result<(), String> {
        let (m, n) = self.size();
        if x.len() != n || y.len() != m {
            return Err(format!(
                "matvec dimension mismatch: matrix {}×{}, x len {}, y len {}",
                m, n, x.len(), y.len()
            ));
        }
        *y = self * x;
        Ok(())
    }

    fn matvec_adjoint(
        &self,
        x: &DVector<Complex64>,
        y: &mut DVector<Complex64>,
    ) -> Result<(), String> {
        let (m, n) = self.size();
        if x.len() != m || y.len() != n {
            return Err(format!(
                "matvec_adjoint dimension mismatch: matrix {}×{}, x len {}, y len {}",
                m, n, x.len(), y.len()
            ));
        }
        *y = self.adjoint() * x;
        Ok(())
    }

    fn diagonal(&self) -> Option<DVector<Complex64>> {
        Some(self.diagonal())
    }

    fn density(&self) -> f64 {
        1.0  // Dense
    }
}

// ---------------------------------------------------------------------------
// Solver trait
// ---------------------------------------------------------------------------

/// Generic linear solver interface.
///
/// Implementations: GMRES, BiCGSTAB, CG, AMG, direct solver, etc.
pub trait LinearSolver<T: ComplexField> {
    /// Solve A*x = b for x.
    ///
    /// # Errors
    ///
    /// Returns Err if solver fails to converge or encounters invalid input.
    fn solve(
        &mut self,
        op: &dyn LinearOperator<T>,
        b: &DVector<T>,
    ) -> Result<DVector<T>, String>;

    /// Return iteration count from last solve (if applicable).
    fn iterations(&self) -> Option<usize> {
        None
    }

    /// Return residual norm from last solve (if applicable).
    fn residual(&self) -> Option<f64> {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dmatrix_real_matvec() {
        let a = DMatrix::from_fn(3, 2, |i, j| (i + j) as f64);
        let x = DVector::from_vec(vec![1.0, 2.0]);
        let mut y = DVector::zeros(3);

        a.matvec(&x, &mut y).unwrap();
        // a = [[0, 1], [1, 2], [2, 3]]
        // y = a * x = [0*1 + 1*2, 1*1 + 2*2, 2*1 + 3*2] = [2, 5, 8]
        assert_eq!(y[0], 2.0);
        assert_eq!(y[1], 5.0);
        assert_eq!(y[2], 8.0);
    }

    #[test]
    fn test_dmatrix_complex_matvec() {
        let a = DMatrix::from_fn(2, 2, |i, j| Complex64::new((i + j) as f64, 0.0));
        let x = DVector::from_vec(vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 1.0),
        ]);
        let mut y = DVector::zeros(2);

        a.matvec(&x, &mut y).unwrap();
        // a = [[0, 1], [1, 2]]
        // y = a * x = [0*1 + 1*(0+i), 1*1 + 2*(0+i)] = [i, 1+2i]
        assert_eq!(y[0], Complex64::new(0.0, 1.0));
        assert_eq!(y[1], Complex64::new(1.0, 2.0));
    }

    #[test]
    fn test_dmatrix_adjoint() {
        let a = DMatrix::from_fn(2, 3, |i, j| Complex64::new(i as f64, j as f64));
        let x = DVector::from_vec(vec![
            Complex64::new(1.0, 0.0),
            Complex64::new(0.0, 1.0),
        ]);
        let mut y = DVector::zeros(3);

        a.matvec_adjoint(&x, &mut y).unwrap();
        // y = A^H * x
        assert_eq!(y.len(), 3);
    }

    #[test]
    fn test_size_adjoint() {
        let a = DMatrix::<f64>::zeros(5, 3);
        assert_eq!(a.size(), (5, 3));
        assert_eq!(a.size_adjoint(), (3, 5));
    }

    #[test]
    fn test_dimension_mismatch() {
        let a = DMatrix::<f64>::zeros(3, 2);
        let x = DVector::zeros(3);  // Wrong dimension
        let mut y = DVector::zeros(3);

        let result = a.matvec(&x, &mut y);
        assert!(result.is_err());
    }
}
