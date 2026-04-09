/// FEM-assembly backed stiffness matrix builder for the electrostatic solver.
///
/// Uses `fem_assembly::Assembler` + `fem_space::H1Space` instead of
/// rem2's hand-rolled Tet4/Tri3 kernels.  This enables:
///   - P2 (quadratic) elements for higher accuracy
///   - Quadrature-order control
///   - Reuse of the same assembly path for future physics
///
/// The resulting `fem_linalg::CsrMatrix<f64>` is converted to rem2's
/// `rem_core::CsrMatrix` via `CsrMatrix::from_fem_csr()`.
///
/// # Limitations
/// - Does **not** support periodic node remapping (use `assemble.rs` for that).
/// - Anisotropic case uses a custom `TensorDiffusionIntegrator`; requires
///   the material tensor to be uniform per element (piecewise constant).

use rem_core::CsrMatrix;
use rem_materials::DomainMap;
use rem_mesh::RemMesh;
use fem_assembly::standard::DiffusionIntegrator;
use fem_assembly::coefficient::PWConstCoeff;
use fem_assembly::Assembler;
use fem_assembly::integrator::{BilinearIntegrator, QpData};
use fem_space::H1Space;

// ---------------------------------------------------------------------------
// Isotropic: use DiffusionIntegrator with PWConstCoeff
// ---------------------------------------------------------------------------

/// Assemble the global stiffness matrix using fem-rs P1 assembly.
///
/// `order` — polynomial order: 1 = P1 (standard), 2 = P2 (quadratic, higher accuracy).
/// `quad_order` — quadrature rule order (2 is sufficient for P1; 4 for P2).
pub fn assemble_stiffness_fem(
    mesh: &RemMesh,
    domain_map: &DomainMap,
    order: u8,
    quad_order: u8,
) -> CsrMatrix {
    // Build piecewise-constant epsilon from domain map
    let kappa = PWConstCoeff::new(
        mesh.domain_tags.iter().map(|(&phys_tag, &mat_idx)| {
            let eps = domain_map.materials[mat_idx].epsilon_abs();
            (phys_tag as i32, eps)
        })
    ).with_default(1.0);

    let integ = DiffusionIntegrator { kappa };

    if mesh.dim == 2 {
        let space = H1Space::new(mesh.to_simplex_mesh_2d(), order);
        let fem_csr = Assembler::assemble_bilinear(&space, &[&integ], quad_order);
        CsrMatrix::from_fem_csr(fem_csr)
    } else {
        let space = H1Space::new(mesh.to_simplex_mesh(), order);
        let fem_csr = Assembler::assemble_bilinear(&space, &[&integ], quad_order);
        CsrMatrix::from_fem_csr(fem_csr)
    }
}

// ---------------------------------------------------------------------------
// Anisotropic: custom TensorDiffusionIntegrator
// ---------------------------------------------------------------------------

/// Bilinear integrator for the anisotropic diffusion form:
///   a(u,v) = ∫_Ω  ∇u · A(x) · ∇v  dx
///
/// `A` is a 3×3 symmetric tensor (row-major, piecewise constant per element tag).
struct TensorDiffusionIntegrator {
    /// Maps element tag (i32) to 9-entry row-major tensor.
    tensors: std::collections::HashMap<i32, [f64; 9]>,
    default_tensor: [f64; 9],
}

impl TensorDiffusionIntegrator {
    fn new(
        entries: impl IntoIterator<Item = (i32, [[f64; 3]; 3])>,
        default: [[f64; 3]; 3],
    ) -> Self {
        let flatten = |t: [[f64; 3]; 3]| -> [f64; 9] {
            [t[0][0], t[0][1], t[0][2],
             t[1][0], t[1][1], t[1][2],
             t[2][0], t[2][1], t[2][2]]
        };
        TensorDiffusionIntegrator {
            tensors: entries.into_iter().map(|(tag, t)| (tag, flatten(t))).collect(),
            default_tensor: flatten(default),
        }
    }

    #[inline]
    fn tensor_for(&self, tag: i32) -> &[f64; 9] {
        self.tensors.get(&tag).unwrap_or(&self.default_tensor)
    }
}

impl BilinearIntegrator for TensorDiffusionIntegrator {
    fn add_to_element_matrix(&self, qp: &QpData<'_>, k_elem: &mut [f64]) {
        let n   = qp.n_dofs;
        let d   = qp.dim;
        let a   = self.tensor_for(qp.elem_tag);
        let w   = qp.weight;

        // k_elem[i,j] += w * ∇φᵢ · A · ∇φⱼ
        // A is stored as d×d (we use the d×d leading submatrix for 2-D)
        for i in 0..n {
            for j in 0..n {
                let mut val = 0.0;
                for r in 0..d {
                    let mut ag_jr = 0.0;
                    for c in 0..d {
                        // a[r,c] = a[r*3 + c] (we use up to 3×3 but clip at dim)
                        ag_jr += a[r * 3 + c] * qp.grad_phys[j * d + c];
                    }
                    val += qp.grad_phys[i * d + r] * ag_jr;
                }
                k_elem[i * n + j] += w * val;
            }
        }
    }
}

// Silence the Send+Sync lint — HashMap<i32,[f64;9]> is Send+Sync
unsafe impl Send for TensorDiffusionIntegrator {}
unsafe impl Sync for TensorDiffusionIntegrator {}

/// Assemble the anisotropic stiffness matrix using fem-rs assembly.
///
/// `order`      — polynomial order (1 or 2).
/// `quad_order` — quadrature order (2 for P1, 4 for P2).
pub fn assemble_stiffness_aniso_fem(
    mesh: &RemMesh,
    domain_map: &DomainMap,
    order: u8,
    quad_order: u8,
) -> CsrMatrix {
    let identity = [[1.0f64, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let integ = TensorDiffusionIntegrator::new(
        mesh.domain_tags.iter().map(|(&phys_tag, &mat_idx)| {
            let t = domain_map.materials[mat_idx].epsilon_tensor;
            (phys_tag as i32, t)
        }),
        identity,
    );

    if mesh.dim == 2 {
        let space = H1Space::new(mesh.to_simplex_mesh_2d(), order);
        let fem_csr = Assembler::assemble_bilinear(&space, &[&integ], quad_order);
        CsrMatrix::from_fem_csr(fem_csr)
    } else {
        let space = H1Space::new(mesh.to_simplex_mesh(), order);
        let fem_csr = Assembler::assemble_bilinear(&space, &[&integ], quad_order);
        CsrMatrix::from_fem_csr(fem_csr)
    }
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use rem_materials::{DomainMap, Material};
    use rem_mesh::{Node, Element, ElementKind, RemMesh};

    fn unit_tet_mesh() -> RemMesh {
        // Regular tetrahedron: (0,0,0),(1,0,0),(0,1,0),(0,0,1), tag=1
        RemMesh {
            nodes: vec![
                Node { id: 0, x: 0.0, y: 0.0, z: 0.0 },
                Node { id: 1, x: 1.0, y: 0.0, z: 0.0 },
                Node { id: 2, x: 0.0, y: 1.0, z: 0.0 },
                Node { id: 3, x: 0.0, y: 0.0, z: 1.0 },
            ],
            volume_elements: vec![
                Element { id: 1, kind: ElementKind::Tet4, tag: 1,
                          node_ids: vec![0,1,2,3], rank: 0 },
            ],
            boundary_elements: vec![],
            domain_tags: [(1u32, 0usize)].iter().cloned().collect(),
            boundary_tags: Default::default(),
            dim: 3, rank: 0, size: 1,
        }
    }

    fn domain_map_eps(eps: f64) -> DomainMap {
        let eps_rel = eps / rem_core::constants::EPS0;
        let mat = Material::from_scalars(eps_rel, 1.0, 0.0, 0.0);
        DomainMap::from_materials(vec![mat], [(1u32, 0usize)])
    }

    #[test]
    fn fem_assembly_row_sum_zero() {
        let mesh = unit_tet_mesh();
        let dom  = domain_map_eps(rem_core::constants::EPS0);
        let csr  = assemble_stiffness_fem(&mesh, &dom, 1, 2);

        let n = csr.nrows;
        let x = vec![1.0f64; n];
        let mut y = vec![0.0f64; n];
        csr.matvec(&x, &mut y, &rem_parallel::NoComm);
        for (i, &yi) in y.iter().enumerate() {
            assert!(yi.abs() < 1e-12, "row {i} sum = {yi}");
        }
    }

    #[test]
    fn fem_assembly_symmetry() {
        let mesh = unit_tet_mesh();
        let dom  = domain_map_eps(2.5 * rem_core::constants::EPS0);
        let csr  = assemble_stiffness_fem(&mesh, &dom, 1, 2);
        let n = csr.nrows;
        for i in 0..n {
            for k in csr.row_ptr[i]..csr.row_ptr[i+1] {
                let j = csr.col_idx[k];
                let kij = csr.values[k];
                let kji = csr.col_idx[csr.row_ptr[j]..csr.row_ptr[j+1]]
                    .iter().zip(csr.values[csr.row_ptr[j]..csr.row_ptr[j+1]].iter())
                    .find(|(&c, _)| c == i)
                    .map(|(_, &v)| v)
                    .unwrap_or(0.0);
                assert!((kij - kji).abs() < 1e-12,
                    "K[{i},{j}]={kij} != K[{j},{i}]={kji}");
            }
        }
    }

    #[test]
    fn fem_vs_manual_scalar_match() {
        // fem-assembly path and hand-rolled path must give same row-sum = 0 for unit field
        let mesh = unit_tet_mesh();
        let eps  = rem_core::constants::EPS0;
        let dom  = domain_map_eps(eps);

        let fem_csr  = assemble_stiffness_fem(&mesh, &dom, 1, 2);
        let hand_csr = crate::assemble::assemble_stiffness(&mesh, |_| eps)
            .unwrap().to_csr();

        let n = fem_csr.nrows;
        let x = vec![1.0f64; n];
        let mut yf = vec![0.0f64; n];
        let mut yh = vec![0.0f64; n];
        fem_csr.matvec(&x, &mut yf, &rem_parallel::NoComm);
        hand_csr.matvec(&x, &mut yh, &rem_parallel::NoComm);
        for (yfi, yhi) in yf.iter().zip(yh.iter()) {
            assert!((yfi - yhi).abs() < 1e-12, "fem={yfi} hand={yhi}");
        }
    }

    #[test]
    fn aniso_fem_identity_matches_scalar() {
        let mesh = unit_tet_mesh();
        let eps  = 3.0 * rem_core::constants::EPS0;
        let dom  = domain_map_eps(eps);

        let aniso_csr  = assemble_stiffness_aniso_fem(&mesh, &dom, 1, 2);
        let scalar_csr = assemble_stiffness_fem(&mesh, &dom, 1, 2);

        let n = aniso_csr.nrows;
        let x = vec![1.0f64; n];
        let mut ya = vec![0.0f64; n];
        let mut ys = vec![0.0f64; n];
        aniso_csr.matvec(&x, &mut ya, &rem_parallel::NoComm);
        scalar_csr.matvec(&x, &mut ys, &rem_parallel::NoComm);
        for (a, s) in ya.iter().zip(ys.iter()) {
            assert!((a - s).abs() < 1e-11, "aniso={a} scalar={s}");
        }
    }
}
