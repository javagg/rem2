// MoM impedance matrix fill — EFIE / MFIE kernels.
//
// kernel_type == 0: EFIE  — scalar Green's function G
// kernel_type == 1: MFIE  — gradient Green's function ∇G (vector)
//
// When green_table_n > 0: G(ρ) from log-spaced lookup table (layered media).
// For MFIE with table: ∇G_horiz from table finite-diff, ∇G_z ≈ 0.

struct Params {
    n: u32,
    quad_n: u32,
    omega_mu0: f32,
    k_divisor: f32,
    k: f32,
    green_table_n: u32,
    green_rho_min: f32,
    green_rho_max: f32,
    kernel_type: u32,
};

struct FaceData {
    nodes: vec3<u32>,
    centroid: vec4<f32>,
    normal: vec4<f32>,
    area: f32,
};

struct BasisData {
    edge_idx: u32,
    plus_face: u32,
    minus_face: u32,
    free_node_plus: u32,
    free_node_minus: u32,
    length: f32,
};

@group(0) @binding(0) var<storage, read> nodes:        array<vec4<f32>>;
@group(0) @binding(1) var<storage, read> faces:        array<FaceData>;
@group(0) @binding(2) var<storage, read> bases:        array<BasisData>;
@group(0) @binding(3) var<storage, read> quad_data:    array<vec4<f32>>;  // bary.xyz + weight.w
@group(0) @binding(4) var<storage, read>   params:       Params;
@group(0) @binding(5) var<storage, read_write> out_re:  array<f32>;
@group(0) @binding(6) var<storage, read_write> out_im:  array<f32>;
@group(0) @binding(7) var<storage, read> green_table:   array<vec4<f32>>;

const INV_4PI: f32 = 0.07957747154594767;

fn dist3(ax: f32, ay: f32, az: f32, bx: f32, by: f32, bz: f32) -> f32 {
    let dx = ax - bx; let dy = ay - by; let dz = az - bz;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

fn dist2(ax: f32, ay: f32, bx: f32, by: f32) -> f32 {
    let dx = ax - bx; let dy = ay - by;
    return sqrt(dx * dx + dy * dy);
}

fn table_lookup(rho: f32) -> vec2<f32> {
    if (rho < 1e-7) { return vec2<f32>(0.0); }
    let log_min = params.green_rho_min;
    let log_max = params.green_rho_max;
    let t = (log(rho) - log_min) / (log_max - log_min);
    let idx_f = t * f32(params.green_table_n - 1u);
    let i0 = min(u32(idx_f), params.green_table_n - 1u);
    let i1 = min(i0 + 1u, params.green_table_n - 1u);
    let f = idx_f - f32(i0);
    return mix(green_table[i0].xy, green_table[i1].xy, vec2<f32>(f, f));
}

fn green_g(rm: vec4<f32>, rn: vec4<f32>) -> vec2<f32> {
    if (params.green_table_n > 0u) {
        return table_lookup(dist2(rm.x, rm.y, rn.x, rn.y));
    }
    let R = dist3(rm.x, rm.y, rm.z, rn.x, rn.y, rn.z);
    if (R < 1e-7) { return vec2<f32>(0.0); }
    let phase = params.k * R;
    return vec2<f32>(cos(phase), -sin(phase)) * (INV_4PI / R);
}

// ∇G component along (dr_c): returns (re, im)
fn green_dg_dc(rm: vec4<f32>, rn: vec4<f32>, dr: f32, R: f32, rho: f32) -> vec2<f32> {
    if (params.green_table_n > 0u) {
        // ∂G/∂c = dG_drho * dr_c / rho   (horizontal only; vertical ≈ 0)
        if (rho < 1e-7 || abs(dr) < 1e-14) { return vec2<f32>(0.0); }
        let h = rho * 0.001;
        let gp = table_lookup(rho + h);
        let gm = table_lookup(rho - h);
        return (gp - gm) / (2.0 * h) * (dr / rho);
    }
    // Analytic: ∇G_c = -G · (1/R + jk)/R · dr_c
    if (R < 1e-7) { return vec2<f32>(0.0); }
    let g = green_g(rm, rn);
    let invR = 1.0 / R;
    let factor_re = -(g.x * invR - g.y * params.k) * invR;
    let factor_im = -(g.y * invR + g.x * params.k) * invR;
    return vec2<f32>(factor_re * dr, factor_im * dr);
}

fn global_point(bary: vec4<f32>, fi: u32) -> vec4<f32> {
    let f = faces[fi];
    let n0 = nodes[f.nodes.x];
    let n1 = nodes[f.nodes.y];
    let n2 = nodes[f.nodes.z];
    return vec4<f32>(
        bary.x * n0.x + bary.y * n1.x + bary.z * n2.x,
        bary.x * n0.y + bary.y * n1.y + bary.z * n2.y,
        bary.x * n0.z + bary.y * n1.z + bary.z * n2.z, 0.0);
}

fn mfie_overlap(m: u32, bm: BasisData, qn: u32) -> f32 {
    var overlap = 0.0f;
    for (var pi = 0u; pi < 2u; pi++) {
        let mfi = select(bm.minus_face, bm.plus_face, pi == 0u);
        let m_plus = pi == 0u;
        let f = faces[mfi];
        let free = select(nodes[bm.free_node_minus], nodes[bm.free_node_plus], m_plus);
        let sc = select(-bm.length / (2.0 * f.area), bm.length / (2.0 * f.area), m_plus);
        for (var qi = 0u; qi < qn; qi++) {
            let pt = quad_data[qi];
            let r = global_point(pt, mfi);
            let fx = sc * (r.x - free.x);
            let fy = sc * (r.y - free.y);
            let fz = sc * (r.z - free.z);
            overlap += (fx*fx + fy*fy + fz*fz) * pt.w * 2.0 * f.area;
        }
    }
    return overlap;
}

@compute @workgroup_size(8, 8)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    let m = id.x;
    let n = id.y;
    if (m >= params.n || n >= params.n) { return; }

    let bm = bases[m];
    let bn = bases[n];
    let qn = params.quad_n;
    let is_mfie = params.kernel_type == 1u;

    var efie_re = 0.0f;
    var efie_im = 0.0f;
    var mfie_re = 0.0f;
    var mfie_im = 0.0f;

    for (var p = 0u; p < 2u; p++) {
        let mfi    = select(bm.minus_face, bm.plus_face,  p == 0u);
        let m_free = select(bm.free_node_minus, bm.free_node_plus, p == 0u);
        let m_sign = select(-1.0, 1.0, p == 0u);

        let fm = faces[mfi];
        let free_m = nodes[m_free];
        let area_m = fm.area;
        let scale_m = m_sign * bm.length / (2.0 * area_m);
        let div_m   = m_sign * bm.length / area_m;
        let nm_x = fm.normal.x;
        let nm_y = fm.normal.y;
        let nm_z = fm.normal.z;

        for (var q = 0u; q < 2u; q++) {
            let nfi    = select(bn.minus_face, bn.plus_face,  q == 0u);
            let n_free = select(bn.free_node_minus, bn.free_node_plus, q == 0u);
            let n_sign = select(-1.0, 1.0, q == 0u);

            let f_n = faces[nfi];
            let free_n = nodes[n_free];
            let area_n = f_n.area;
            let scale_n = n_sign * bn.length / (2.0 * area_n);
            let div_n   = n_sign * bn.length / area_n;

            let area_scale = 4.0 * area_m * area_n;
            let scalar_part = params.k_divisor * div_m * div_n;

                for (var qi = 0u; qi < qn; qi++) {
                    let bm_pt = quad_data[qi];
                    let wm = bm_pt.w;
                    let rm = global_point(bm_pt, mfi);

                    let fm_x = scale_m * (rm.x - free_m.x);
                    let fm_y = scale_m * (rm.y - free_m.y);
                    let fm_z = scale_m * (rm.z - free_m.z);

                    for (var qj = 0u; qj < qn; qj++) {
                        let bn_pt = quad_data[qj];
                        let wn = bn_pt.w;
                    let rn = global_point(bn_pt, nfi);

                    let fn_x = scale_n * (rn.x - free_n.x);
                    let fn_y = scale_n * (rn.y - free_n.y);
                    let fn_z = scale_n * (rn.z - free_n.z);

                    let dx = rm.x - rn.x;
                    let dy = rm.y - rn.y;
                    let dz = rm.z - rn.z;
                    let R   = sqrt(dx*dx + dy*dy + dz*dz);
                    let rho = sqrt(dx*dx + dy*dy);

                    // EFIE: G · (f_m·f_n - (1/k²)·div_m·div_n)
                    let dot_ff = fm_x * fn_x + fm_y * fn_y + fm_z * fn_z;
                    let kern   = dot_ff - scalar_part;
                    let g = green_g(rm, rn);
                    let wgt = wm * wn * area_scale;
                    efie_re += (g.x * kern) * wgt;
                    efie_im += (g.y * kern) * wgt;

                    // MFIE: n̂_m · [∇G × f_n] × f_m
                    let dgx = green_dg_dc(rm, rn, dx, R, rho);
                    let dgy = green_dg_dc(rm, rn, dy, R, rho);
                    let dgz = green_dg_dc(rm, rn, dz, R, rho);

                    // ∇G × f_n  (complex vector × real vector)
                    let cgx_re = dgy.x * fn_z - dgz.x * fn_y;
                    let cgx_im = dgy.y * fn_z - dgz.y * fn_y;
                    let cgy_re = dgz.x * fn_x - dgx.x * fn_z;
                    let cgy_im = dgz.y * fn_x - dgx.y * fn_z;
                    let cgz_re = dgx.x * fn_y - dgy.x * fn_x;
                    let cgz_im = dgx.y * fn_y - dgy.y * fn_x;

                    // n̂ × (∇G × f_n)  (real vector × complex vector)
                    let ncg_x_re = nm_y * cgz_re - nm_z * cgy_re;
                    let ncg_x_im = nm_y * cgz_im - nm_z * cgy_im;
                    let ncg_y_re = nm_z * cgx_re - nm_x * cgz_re;
                    let ncg_y_im = nm_z * cgx_im - nm_x * cgz_im;
                    let ncg_z_re = nm_x * cgy_re - nm_y * cgx_re;
                    let ncg_z_im = nm_x * cgy_im - nm_y * cgx_im;

                    // f_m · (n̂ × (∇G × f_n))
                    let dc_re = fm_x * ncg_x_re + fm_y * ncg_y_re + fm_z * ncg_z_re;
                    let dc_im = fm_x * ncg_x_im + fm_y * ncg_y_im + fm_z * ncg_z_im;

                    mfie_re += dc_re * wgt;
                    mfie_im += dc_im * wgt;
                }
            }
        }
    }

    let idx = m + n * params.n;
    if (is_mfie) {
        let ident = select(0.0, mfie_overlap(m, bm, qn), m == n);
        out_re[idx] = 0.5 * ident + mfie_re;
        out_im[idx] = mfie_im;
    } else {
        out_re[idx] = params.omega_mu0 * efie_im;
        out_im[idx] = -params.omega_mu0 * efie_re;
    }
}
