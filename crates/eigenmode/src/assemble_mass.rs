//! P1 finite element mass matrix assembly.
//!
//! M_ij = Σ_e ε_e ∫_Ωe φ_i φ_j dΩ
//!
//! For Tri3: using exact integration of P1 products over a triangle.
//! For Tet4: using exact integration over a tetrahedron.

use rem_core::{TripletMatrix, RemError, RemResult};
use rem_mesh::{RemMesh, ElementKind};

/// Assemble the global mass matrix with a scalar coefficient per element.
///
/// `coeff_fn` maps physical group tag → permittivity ε (or similar coefficient).
pub fn assemble_mass(
    mesh: &RemMesh,
    coeff_fn: impl Fn(u32) -> f64,
) -> RemResult<TripletMatrix> {
    let n = mesh.n_nodes();
    let cap = mesh.n_volume_elements() * 16;
    let mut triplet = TripletMatrix::with_capacity(n, n, cap);

    for elem in &mesh.volume_elements {
        if mesh.size > 1 && elem.rank != mesh.rank {
            continue;
        }
        let eps = coeff_fn(elem.tag);
        match elem.kind {
            ElementKind::Tri3 => mass_tri3(mesh, elem, eps, &mut triplet)?,
            ElementKind::Tet4 => mass_tet4(mesh, elem, eps, &mut triplet)?,
            ElementKind::Tet10 => mass_tet4_corners(mesh, elem, eps, &mut triplet)?,
            other => {
                log::warn!("Mass matrix: element {:?} not supported — skipping", other);
            }
        }
    }

    Ok(triplet)
}

/// Assemble the global mass matrix with a full anisotropic tensor coefficient.
///
/// The mass integral `M_ij = ∫ ε̄ φ_i φ_j dΩ` uses a scalar effective permittivity
/// derived from the tensor as `ε̄ = trace(A) / 3`.  For isotropic rotation this equals
/// ε₀εᵣ exactly; for genuinely anisotropic media it gives the isotropic equivalent.
///
/// `tensor_fn` maps physical group tag → absolute permittivity tensor [F/m], 3×3 row-major.
pub fn assemble_mass_aniso(
    mesh: &RemMesh,
    tensor_fn: impl Fn(u32) -> [[f64; 3]; 3],
) -> RemResult<TripletMatrix> {
    assemble_mass(mesh, |tag| {
        let a = tensor_fn(tag);
        (a[0][0] + a[1][1] + a[2][2]) / 3.0
    })
}

fn mass_tri3(
    mesh: &RemMesh,
    elem: &rem_mesh::Element,
    eps: f64,
    triplet: &mut TripletMatrix,
) -> RemResult<()> {
    debug_assert_eq!(elem.node_ids.len(), 3);
    let [n0, n1, n2] = [elem.node_ids[0], elem.node_ids[1], elem.node_ids[2]];
    let (x0, y0) = (mesh.nodes[n0].x, mesh.nodes[n0].y);
    let (x1, y1) = (mesh.nodes[n1].x, mesh.nodes[n1].y);
    let (x2, y2) = (mesh.nodes[n2].x, mesh.nodes[n2].y);

    let det_j = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
    let area = 0.5 * det_j.abs();
    if area < 1e-300 {
        return Err(RemError::Mesh(format!(
            "Degenerate Tri3 element {} in mass assembly", elem.id
        )));
    }

    // ∫_Tri φ_i φ_j dA = area * (1 + δ_{ij}) / 12   (exact for linear triangle)
    let nodes = [n0, n1, n2];
    let diag = eps * area / 6.0;
    let off  = eps * area / 12.0;

    for i in 0..3 {
        for j in 0..3 {
            let val = if i == j { diag } else { off };
            triplet.add(nodes[i], nodes[j], val);
        }
    }
    Ok(())
}

fn mass_tet4(
    mesh: &RemMesh,
    elem: &rem_mesh::Element,
    eps: f64,
    triplet: &mut TripletMatrix,
) -> RemResult<()> {
    mass_tet4_first4(mesh, &elem.node_ids, elem.id, eps, triplet)
}

/// Use first 4 (corner) nodes of a Tet10 element for P1 approximation.
fn mass_tet4_corners(
    mesh: &RemMesh,
    elem: &rem_mesh::Element,
    eps: f64,
    triplet: &mut TripletMatrix,
) -> RemResult<()> {
    mass_tet4_first4(mesh, &elem.node_ids[..4], elem.id, eps, triplet)
}

fn mass_tet4_first4(
    mesh: &RemMesh,
    node_ids: &[usize],
    elem_id: usize,
    eps: f64,
    triplet: &mut TripletMatrix,
) -> RemResult<()> {
    let [n0, n1, n2, n3] = [node_ids[0], node_ids[1], node_ids[2], node_ids[3]];
    let nodes = [n0, n1, n2, n3];

    let x = [mesh.nodes[n0].x, mesh.nodes[n1].x, mesh.nodes[n2].x, mesh.nodes[n3].x];
    let y = [mesh.nodes[n0].y, mesh.nodes[n1].y, mesh.nodes[n2].y, mesh.nodes[n3].y];
    let z = [mesh.nodes[n0].z, mesh.nodes[n1].z, mesh.nodes[n2].z, mesh.nodes[n3].z];

    // det(J) = 6 * vol
    let j = [
        [x[1]-x[0], x[2]-x[0], x[3]-x[0]],
        [y[1]-y[0], y[2]-y[0], y[3]-y[0]],
        [z[1]-z[0], z[2]-z[0], z[3]-z[0]],
    ];
    let det = j[0][0]*(j[1][1]*j[2][2]-j[1][2]*j[2][1])
            - j[0][1]*(j[1][0]*j[2][2]-j[1][2]*j[2][0])
            + j[0][2]*(j[1][0]*j[2][1]-j[1][1]*j[2][0]);
    let vol = det.abs() / 6.0;
    if vol < 1e-300 {
        return Err(RemError::Mesh(format!(
            "Degenerate Tet4 element {} in mass assembly", elem_id
        )));
    }

    // ∫_Tet φ_i φ_j dV = vol * (1 + δ_{ij}) / 20   (exact for linear tetrahedron)
    let diag = eps * vol / 10.0;
    let off  = eps * vol / 20.0;

    for i in 0..4 {
        for j in 0..4 {
            let val = if i == j { diag } else { off };
            triplet.add(nodes[i], nodes[j], val);
        }
    }
    Ok(())
}
