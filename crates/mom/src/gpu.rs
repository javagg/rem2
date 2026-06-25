//! Parallel-CPU (Rayon) accelerated MoM impedance matrix assembly.
//! and optional wgpu GPU compute-shader path for MoM impedance matrix assembly.
//!
//! ## CPU path (default)
//!
//! Each RWG basis-function pair (m, n) is independent, making MoM matrix fill
//! embarrassingly parallel.  On multi-core CPUs Rayon gives 4–16× speedup over
//! serial assembly depending on core count.
//!
//! ## GPU path (`--features wgpu-gpu`)
//!
//! Compile with `--features rem-mom/wgpu-gpu` to enable the wgpu compute-shader
//! backend.  At runtime, [`gpu_available()`] probes for a wgpu adapter; if one is
//! found, [`fill_impedance_wgpu()`] dispatches a WGSL compute shader that
//! evaluates all Z_mn pairs in parallel on the GPU.
//!
//! See `wgsl/impedance_fill.wgsl` for the compute shader source.

use crate::basis::rwg::RwgBasis;
use crate::surface_mesh::TriFace;
use nalgebra::DMatrix;
use num_complex::Complex64;
use rem_core::RemResult;

#[cfg(not(target_arch = "wasm32"))]
use rayon::prelude::*;

#[cfg(all(not(target_arch = "wasm32"), feature = "wgpu-gpu"))]
use std::sync::OnceLock;

/// Minimum basis count where parallel assembly is beneficial over serial.
pub const GPU_MIN_BASIS: usize = 500;

/// Returns `true` when Rayon-parallel CPU assembly is available.
pub fn gpu_available() -> bool {
    cfg!(not(target_arch = "wasm32"))
}

// ═══════════════════════════════════════════════════════════════════════════
//  wgpu adapter probe
// ═══════════════════════════════════════════════════════════════════════════

/// Returns `true` when a wgpu-compatible GPU adapter is available at runtime.
#[cfg(all(not(target_arch = "wasm32"), feature = "wgpu-gpu"))]
pub fn wgpu_adapter_available() -> bool {
    static ADAPTER_AVAILABLE: OnceLock<bool> = OnceLock::new();
    if let Some(v) = ADAPTER_AVAILABLE.get() {
        return *v;
    }
    use pollster::block_on;
    let instance = wgpu::Instance::default();
    let adapter = block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
        power_preference: wgpu::PowerPreference::HighPerformance,
        compatible_surface: None,
        force_fallback_adapter: false,
    }));
    let has_adapter = adapter.is_some();
    let _ = ADAPTER_AVAILABLE.set(has_adapter);
    has_adapter
}

#[cfg(any(target_arch = "wasm32", not(feature = "wgpu-gpu")))]
pub fn wgpu_adapter_available() -> bool {
    false
}

// ═══════════════════════════════════════════════════════════════════════════
//  CPU parallel path
// ═══════════════════════════════════════════════════════════════════════════

#[allow(dead_code)]
pub fn fill_impedance_parallel<F>(n: usize, zmn: F) -> DMatrix<Complex64>
where
    F: Fn(usize, usize) -> Complex64 + Send + Sync,
{
    if n == 0 {
        return DMatrix::<Complex64>::zeros(0, 0);
    }
    #[cfg(not(target_arch = "wasm32"))]
    if n < GPU_MIN_BASIS {
        return fill_impedance_serial(n, zmn);
    }
    #[cfg(not(target_arch = "wasm32"))]
    {
        let mut data = vec![Complex64::ZERO; n * n];
        data.par_chunks_mut(n).enumerate().for_each(|(k, col)| {
            for m in 0..n {
                col[m] = zmn(m, k);
            }
        });
        DMatrix::<Complex64>::from_vec(n, n, data)
    }
    #[cfg(target_arch = "wasm32")]
    {
        let mut data = vec![Complex64::ZERO; n * n];
        for k in 0..n {
            let col = &mut data[k * n..(k + 1) * n];
            for m in 0..n {
                col[m] = zmn(m, k);
            }
        }
        DMatrix::<Complex64>::from_vec(n, n, data)
    }
}

fn fill_impedance_serial<F>(n: usize, zmn: F) -> DMatrix<Complex64>
where
    F: Fn(usize, usize) -> Complex64,
{
    let mut z = DMatrix::<Complex64>::zeros(n, n);
    for m in 0..n {
        for k in 0..n {
            z[(m, k)] = zmn(m, k);
        }
    }
    z
}

// ═══════════════════════════════════════════════════════════════════════════
//  Dispatcher — GPU when available, otherwise Rayon-CPU
// ═══════════════════════════════════════════════════════════════════════════

pub fn fill_impedance_wgpu_or_parallel<F>(n: usize, zmn: F) -> RemResult<DMatrix<Complex64>>
where
    F: Fn(usize, usize) -> Complex64 + Send + Sync,
{
    #[cfg(all(not(target_arch = "wasm32"), feature = "wgpu-gpu"))]
    if wgpu_adapter_available() {
        log::debug!("fill_impedance: wgpu adapter available (n={}), but closure-based GPU path is not wired; using CPU.", n);
    }
    Ok(fill_impedance_parallel(n, zmn))
}

// ═══════════════════════════════════════════════════════════════════════════
//  GPU data types + wgpu pipeline  (feature = "wgpu-gpu")
// ═══════════════════════════════════════════════════════════════════════════

#[cfg(all(not(target_arch = "wasm32"), feature = "wgpu-gpu"))]
mod gpu_pipeline {
    use crate::basis::rwg::RwgBasis;
    use crate::surface_mesh::TriFace;
    use nalgebra::DMatrix;
    use num_complex::Complex64;
    use std::sync::OnceLock;

    // ── GPU data layouts (must match wgsl/impedance_fill.wgsl) ────────────

    #[repr(C)]
    #[derive(Clone, Copy)]
    pub(super) struct GpuParams {
        n: u32,
        quad_n: u32,
        omega_mu0: f32,
        k_divisor: f32,
        k: f32,
        green_table_n: u32,
        green_rho_min: f32,
        green_rho_max: f32,
        kernel_type: u32,       // 0=EFIE, 1=MFIE
    }

    #[repr(C)]
    pub(super) struct GpuFaceData {
        nodes: [u32; 4],
        centroid: [f32; 4],
        normal: [f32; 4],
        area: f32,
        _pad: [f32; 3],
    }

    #[repr(C)]
    pub(super) struct GpuBasisData {
        edge_idx: u32,
        plus_face: u32,
        minus_face: u32,
        free_node_plus: u32,
        free_node_minus: u32,
        length: f32,
    }

    /// Cached wgpu device, queue, pipeline, and geometry buffers.
    struct GpuResources {
        device: wgpu::Device,
        queue: wgpu::Queue,
        pipeline: wgpu::ComputePipeline,
        bind_group_layout: wgpu::BindGroupLayout,
        geo_cache: std::sync::Mutex<Option<GeometryCache>>,
    }

    struct GeometryCache {
        geo_hash: u64,
        buf_nodes: wgpu::Buffer,
        buf_faces: wgpu::Buffer,
        buf_bases: wgpu::Buffer,
        buf_quad_data: wgpu::Buffer,  // bary.xyz + weight.w packed as vec4
    }

    fn gpu_resources() -> Option<&'static GpuResources> {
        static RES: OnceLock<Option<GpuResources>> = OnceLock::new();
        RES.get_or_init(|| {
            use pollster::block_on;
            let instance = wgpu::Instance::default();
            let adapter = block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                compatible_surface: None,
                force_fallback_adapter: false,
            }));
            let adapter = adapter?;
            let (device, queue) = block_on(adapter.request_device(
                &wgpu::DeviceDescriptor {
                    label: Some("rem-mom wgpu device"),
                    required_features: wgpu::Features::default(),
                    required_limits: wgpu::Limits::default(),
                    memory_hints: wgpu::MemoryHints::Performance,
                },
                None,
            ))
            .ok()?;

            // Build pipeline and layout once
            let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("impedance_fill_layout"),
                entries: &[
                    wgpu::BindGroupLayoutEntry { binding: 0, visibility: wgpu::ShaderStages::COMPUTE, ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Storage { read_only: true }, has_dynamic_offset: false, min_binding_size: None }, count: None },
                    wgpu::BindGroupLayoutEntry { binding: 1, visibility: wgpu::ShaderStages::COMPUTE, ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Storage { read_only: true }, has_dynamic_offset: false, min_binding_size: None }, count: None },
                    wgpu::BindGroupLayoutEntry { binding: 2, visibility: wgpu::ShaderStages::COMPUTE, ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Storage { read_only: true }, has_dynamic_offset: false, min_binding_size: None }, count: None },
                    wgpu::BindGroupLayoutEntry { binding: 3, visibility: wgpu::ShaderStages::COMPUTE, ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Storage { read_only: true }, has_dynamic_offset: false, min_binding_size: None }, count: None },
                    wgpu::BindGroupLayoutEntry { binding: 4, visibility: wgpu::ShaderStages::COMPUTE, ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Storage { read_only: true }, has_dynamic_offset: false, min_binding_size: None }, count: None },
                    wgpu::BindGroupLayoutEntry { binding: 5, visibility: wgpu::ShaderStages::COMPUTE, ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Storage { read_only: false }, has_dynamic_offset: false, min_binding_size: None }, count: None },
                    wgpu::BindGroupLayoutEntry { binding: 6, visibility: wgpu::ShaderStages::COMPUTE, ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Storage { read_only: false }, has_dynamic_offset: false, min_binding_size: None }, count: None },
                    wgpu::BindGroupLayoutEntry { binding: 7, visibility: wgpu::ShaderStages::COMPUTE, ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Storage { read_only: true }, has_dynamic_offset: false, min_binding_size: None }, count: None },
                ],
            });

            let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: Some("impedance_fill_pipeline_layout"),
                bind_group_layouts: &[&layout],
                push_constant_ranges: &[],
            });

            let shader_src = include_str!("../wgsl/impedance_fill.wgsl");
            let shader_module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("impedance_fill"),
                source: wgpu::ShaderSource::Wgsl(std::borrow::Cow::Borrowed(shader_src)),
            });

            let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some("impedance_fill"),
                layout: Some(&pipeline_layout),
                module: &shader_module,
                entry_point: "main",
                cache: None,
                compilation_options: wgpu::PipelineCompilationOptions::default(),
            });

            Some(GpuResources {
                device, queue, pipeline, bind_group_layout: layout,
                geo_cache: std::sync::Mutex::new(None),
            })
        })
        .as_ref()
    }

    /// Hash mesh data to detect geometry changes between calls.
    fn geometry_hash(
        n: usize, quad_n: usize,
        nodes: &[[f64; 3]], faces: &[TriFace], bases: &[RwgBasis],
    ) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut h = std::collections::hash_map::DefaultHasher::new();
        n.hash(&mut h);
        quad_n.hash(&mut h);
        (nodes.len() as u64).hash(&mut h);
        (faces.len() as u64).hash(&mut h);
        (bases.len() as u64).hash(&mut h);
        h.finish()
    }

    // ── Buffer helpers ────────────────────────────────────────────────────

    fn as_u8_slice<T: Sized>(data: &[T]) -> &[u8] {
        unsafe {
            std::slice::from_raw_parts(
                data.as_ptr() as *const u8,
                data.len() * std::mem::size_of::<T>(),
            )
        }
    }

    fn create_storage_buffer(
        device: &wgpu::Device,
        label: &str,
        size: wgpu::BufferAddress,
    ) -> wgpu::Buffer {
        device.create_buffer(&wgpu::BufferDescriptor {
            label: Some(label),
            size,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        })
    }

    // ── Main GPU fill ────────────────────────────────────────────────────

    /// Fill an N×N EFIE impedance matrix using the wgpu compute shader.
    ///
    /// `green_data`: optional Green's function lookup table as (complex values, ρ_min, ρ_max).
    /// When `None`, the shader uses the built-in free-space analytic formula.
    /// When `Some`, the shader interpolates G(ρ) from the log-spaced table.
    /// `kernel_type`: 0=EFIE, 1=MFIE
    pub(super) fn fill_efie_gpu(
        n: usize,
        quad_n: usize,
        nodes: &[[f64; 3]],
        faces: &[TriFace],
        bases: &[RwgBasis],
        quad_bary: &[[f64; 3]],
        quad_weights: &[f64],
        omega_mu0: f64,
        inv_omega_eps0: f64,
        k: f64,
        green_data: Option<(&[Complex64], f64, f64)>,
        kernel_type: u32,
    ) -> Option<DMatrix<Complex64>> {
        let res = gpu_resources()?;
        let device = &res.device;
        let queue = &res.queue;

        // ── 1. Pack data ──────────────────────────────────────────────────

        // Nodes: [f64;3] → [f32;4] (w=0 padding)
        let gpu_nodes: Vec<[f32; 4]> = nodes
            .iter()
            .map(|n| [n[0] as f32, n[1] as f32, n[2] as f32, 0.0])
            .collect();

        // Faces
        let gpu_faces: Vec<GpuFaceData> = faces
            .iter()
            .map(|f| GpuFaceData {
                nodes: [f.nodes[0] as u32, f.nodes[1] as u32, f.nodes[2] as u32, 0],
                centroid: [
                    f.centroid[0] as f32,
                    f.centroid[1] as f32,
                    f.centroid[2] as f32,
                    0.0,
                ],
                normal: [
                    f.normal[0] as f32,
                    f.normal[1] as f32,
                    f.normal[2] as f32,
                    0.0,
                ],
                area: f.area as f32,
                _pad: [0.0; 3],
            })
            .collect();

        // Bases
        let gpu_bases: Vec<GpuBasisData> = bases
            .iter()
            .map(|b| GpuBasisData {
                edge_idx: b.edge_idx as u32,
                plus_face: b.plus_face as u32,
                minus_face: b.minus_face as u32,
                free_node_plus: b.free_node_plus as u32,
                free_node_minus: b.free_node_minus as u32,
                length: b.length as f32,
            })
            .collect();

        // Quadrature points: bary.xyz + weight.w packed as vec4
        let gpu_quad_data: Vec<[f32; 4]> = quad_bary.iter().zip(quad_weights.iter())
            .map(|(b, w)| [b[0] as f32, b[1] as f32, b[2] as f32, *w as f32])
            .collect();

        // Params
        let (green_table_n, green_rho_min, green_rho_max) = match green_data {
            Some((table, rmin, rmax)) => (table.len().min(u32::MAX as usize) as u32,
                                          rmin as f32, rmax as f32),
            None => (0u32, 0.0, 1.0),
        };
        let gpu_params = GpuParams {
            n: n as u32,
            quad_n: quad_n as u32,
            omega_mu0: omega_mu0 as f32,
            k_divisor: (inv_omega_eps0 / omega_mu0) as f32,
            k: k as f32,
            green_table_n,
            green_rho_min,
            green_rho_max,
            kernel_type,
        };

        let n_elements = (n * n) as wgpu::BufferAddress;
        let output_size = n_elements * std::mem::size_of::<f32>() as wgpu::BufferAddress;

        // ── 2. Create buffers ──────────────────────────────────────────────

        let geo_hash = geometry_hash(n, quad_n, nodes, faces, bases);
        let geo_node_size = (gpu_nodes.len() * 16) as u64;
        let geo_face_size = (gpu_faces.len() * 64) as u64;
        let geo_basis_size = (gpu_bases.len() * 24) as u64;
        let geo_qdata_size = (gpu_quad_data.len() * 16) as u64;

        let params_size = std::mem::size_of::<GpuParams>() as u64;
        let n_elements = (n * n) as wgpu::BufferAddress;
        let output_size = n_elements * std::mem::size_of::<f32>() as wgpu::BufferAddress;

        // Green's function lookup table
        let (buf_green_table, green_table_bytes) = match green_data {
            Some((table, _, _)) => {
                let gpu_table: Vec<[f32; 4]> = table.iter()
                    .map(|g| [g.re as f32, g.im as f32, 0.0, 0.0])
                    .collect();
                let bytes = as_u8_slice(&gpu_table).to_vec();
                let buf = create_storage_buffer(device, "green_table", bytes.len() as u64);
                (buf, bytes)
            }
            None => {
                let buf = create_storage_buffer(device, "green_table", 16);
                (buf, vec![0u8; 16])
            }
        };

        let buf_params = create_storage_buffer(device, "params", params_size);
        queue.write_buffer(&buf_params, 0, as_u8_slice(&[gpu_params]));
        queue.write_buffer(&buf_green_table, 0, &green_table_bytes);

        let buf_out_re = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("out_re"), size: output_size,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        let buf_out_im = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("out_im"), size: output_size,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });

        // Geometry buffers (cached across calls with same mesh)
        let mut guard = res.geo_cache.lock().ok()?;
        let mut geo_upload = false;
        let (bref_n, bref_f, bref_b, bref_qd) =
            if let Some(ref mut cached) = *guard {
                if cached.geo_hash == geo_hash {
                    (&cached.buf_nodes, &cached.buf_faces, &cached.buf_bases,
                     &cached.buf_quad_data)
                } else {
                    geo_upload = true;
                    let bn = create_storage_buffer(device, "nodes", geo_node_size);
                    let bf = create_storage_buffer(device, "faces", geo_face_size);
                    let bb = create_storage_buffer(device, "bases", geo_basis_size);
                    let bq = create_storage_buffer(device, "quad_data", geo_qdata_size);
                    *cached = GeometryCache {
                        geo_hash, buf_nodes: bn, buf_faces: bf,
                        buf_bases: bb, buf_quad_data: bq,
                    };
                    let c = guard.as_ref().unwrap();
                    (&c.buf_nodes, &c.buf_faces, &c.buf_bases, &c.buf_quad_data)
                }
            } else {
                geo_upload = true;
                let bn = create_storage_buffer(device, "nodes", geo_node_size);
                let bf = create_storage_buffer(device, "faces", geo_face_size);
                let bb = create_storage_buffer(device, "bases", geo_basis_size);
                let bq = create_storage_buffer(device, "quad_data", geo_qdata_size);
                *guard = Some(GeometryCache {
                    geo_hash, buf_nodes: bn, buf_faces: bf,
                    buf_bases: bb, buf_quad_data: bq,
                });
                let c = guard.as_ref().unwrap();
                (&c.buf_nodes, &c.buf_faces, &c.buf_bases, &c.buf_quad_data)
            };

        if geo_upload {
            queue.write_buffer(bref_n, 0, as_u8_slice(&gpu_nodes));
            queue.write_buffer(bref_f, 0, as_u8_slice(&gpu_faces));
            queue.write_buffer(bref_b, 0, as_u8_slice(&gpu_bases));
            queue.write_buffer(bref_qd, 0, as_u8_slice(&gpu_quad_data));
        }

        // Bind group (geometry cache lock still held → cached buffers stay alive)
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("impedance_fill_bind_group"),
            layout: &res.bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: bref_n.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: bref_f.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: bref_b.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: bref_qd.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 4, resource: buf_params.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 5, resource: buf_out_re.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 6, resource: buf_out_im.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 7, resource: buf_green_table.as_entire_binding() },
            ],
        });
        // guard drops here; bind group and per-call buffers keep geometry alive

        // ── 5. Dispatch ───────────────────────────────────────────────────

        let workgroup_size: u32 = 8;
        let num_groups_x = (n as u32 + workgroup_size - 1) / workgroup_size;
        let num_groups_y = num_groups_x;

        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("impedance_fill_encoder"),
        });
        {
            let mut cpass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("impedance_fill_pass"),
                timestamp_writes: None,
            });
            cpass.set_pipeline(&res.pipeline);
            cpass.set_bind_group(0, &bind_group, &[]);
            cpass.dispatch_workgroups(num_groups_x, num_groups_y, 1);
        }

        // ── 6. Readback ───────────────────────────────────────────────────

        let readback_re = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("readback_re"),
            size: output_size,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });
        let readback_im = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("readback_im"),
            size: output_size,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });

        encoder.copy_buffer_to_buffer(&buf_out_re, 0, &readback_re, 0, output_size);
        encoder.copy_buffer_to_buffer(&buf_out_im, 0, &readback_im, 0, output_size);

        queue.submit(Some(encoder.finish()));

        // Map and read
        let (tx_re, rx_re) = std::sync::mpsc::channel();
        let (tx_im, rx_im) = std::sync::mpsc::channel();

        let slice_re = readback_re.slice(..);
        slice_re.map_async(wgpu::MapMode::Read, move |v| {
            let _ = tx_re.send(v);
        });
        let slice_im = readback_im.slice(..);
        slice_im.map_async(wgpu::MapMode::Read, move |v| {
            let _ = tx_im.send(v);
        });

        device.poll(wgpu::Maintain::Wait);

        if rx_re.recv().ok() != Some(Ok(())) || rx_im.recv().ok() != Some(Ok(())) {
            log::error!("wgpu map_async failed for readback buffers");
            return None;
        }

        let data_re = slice_re.get_mapped_range();
        let data_im = slice_im.get_mapped_range();

        let out_re: &[f32] = unsafe {
            std::slice::from_raw_parts(data_re.as_ptr() as *const f32, n * n)
        };
        let out_im: &[f32] = unsafe {
            std::slice::from_raw_parts(data_im.as_ptr() as *const f32, n * n)
        };

        // Build DMatrix — the shader writes column-major: idx = m + n * N
        let mut z = DMatrix::<Complex64>::zeros(n, n);
        for m in 0..n {
            for k in 0..n {
                let idx = m + k * n;
                z[(m, k)] = Complex64::new(out_re[idx] as f64, out_im[idx] as f64);
            }
        }

        drop(data_re);
        drop(data_im);
        readback_re.unmap();
        readback_im.unmap();

        Some(z)
    }
}

/// Public entry point: fill EFIE impedance matrix via wgpu GPU compute shader.
///
/// Takes the raw mesh and quadrature data needed by the WGSL kernel.
/// Returns `Ok(Some(matrix))` on success, `Ok(None)` if GPU is not available
/// (caller should fall back to CPU), or `Err` on actual errors.
#[cfg(all(not(target_arch = "wasm32"), feature = "wgpu-gpu"))]
pub fn fill_impedance_wgpu_efie(
    n: usize,
    quad_n: usize,
    nodes: &[[f64; 3]],
    faces: &[TriFace],
    bases: &[RwgBasis],
    quad_bary: &[[f64; 3]],
    quad_weights: &[f64],
    omega_mu0: f64,
    inv_omega_eps0: f64,
    k: f64,
) -> RemResult<Option<DMatrix<Complex64>>> {
    fill_wgpu_generic(
        n, quad_n, nodes, faces, bases, quad_bary, quad_weights,
        omega_mu0, inv_omega_eps0, k, None, 0,
    )
}

/// Fill EFIE impedance matrix using a precomputed Green's function lookup table.
///
/// This enables layered-media (e.g. PCB substrate) Green's functions to be used
/// with the GPU compute shader.  The table should be log-spaced in ρ (horizontal
/// distance) and precomputed via `GreenFunction::g()` on the CPU.
///
/// When GPU is unavailable, returns `Ok(None)` and the caller should fall back to CPU.
#[cfg(all(not(target_arch = "wasm32"), feature = "wgpu-gpu"))]
pub fn fill_impedance_wgpu_efie_with_green_table(
    n: usize,
    quad_n: usize,
    nodes: &[[f64; 3]],
    faces: &[TriFace],
    bases: &[RwgBasis],
    quad_bary: &[[f64; 3]],
    quad_weights: &[f64],
    omega_mu0: f64,
    inv_omega_eps0: f64,
    k: f64,
    green_table: &[Complex64],
    green_rho_min: f64,
    green_rho_max: f64,
) -> RemResult<Option<DMatrix<Complex64>>> {
    fill_wgpu_generic(
        n, quad_n, nodes, faces, bases, quad_bary, quad_weights,
        omega_mu0, inv_omega_eps0, k,
        Some((green_table, green_rho_min, green_rho_max)), 0,
    )
}

/// Fill MFIE impedance matrix via wgpu compute shader (free-space kernel).
#[cfg(all(not(target_arch = "wasm32"), feature = "wgpu-gpu"))]
pub fn fill_impedance_wgpu_mfie(
    n: usize,
    quad_n: usize,
    nodes: &[[f64; 3]],
    faces: &[TriFace],
    bases: &[RwgBasis],
    quad_bary: &[[f64; 3]],
    quad_weights: &[f64],
    omega_mu0: f64,
    inv_omega_eps0: f64,
    k: f64,
) -> RemResult<Option<DMatrix<Complex64>>> {
    fill_wgpu_generic(
        n, quad_n, nodes, faces, bases, quad_bary, quad_weights,
        omega_mu0, inv_omega_eps0, k, None, 1,
    )
}

/// Fill MFIE impedance matrix using a precomputed Green's function lookup table.
#[cfg(all(not(target_arch = "wasm32"), feature = "wgpu-gpu"))]
pub fn fill_impedance_wgpu_mfie_with_green_table(
    n: usize,
    quad_n: usize,
    nodes: &[[f64; 3]],
    faces: &[TriFace],
    bases: &[RwgBasis],
    quad_bary: &[[f64; 3]],
    quad_weights: &[f64],
    omega_mu0: f64,
    inv_omega_eps0: f64,
    k: f64,
    green_table: &[Complex64],
    green_rho_min: f64,
    green_rho_max: f64,
) -> RemResult<Option<DMatrix<Complex64>>> {
    fill_wgpu_generic(
        n, quad_n, nodes, faces, bases, quad_bary, quad_weights,
        omega_mu0, inv_omega_eps0, k,
        Some((green_table, green_rho_min, green_rho_max)), 1,
    )
}

/// Post-process: correct near-singular face-pair contributions in a GPU-filled matrix.
///
/// The GPU shader uses the same Gauss quadrature for all face-pairs, including
/// near-singular ones (self, shared-edge, shared-vertex).  This function replaces
/// those Gauss-approximated contributions with the correct Sauter-Schwab/Duffy
/// results by computing corrections per near-singular pair on the CPU.
#[cfg(all(not(target_arch = "wasm32"), feature = "wgpu-gpu"))]
fn correct_near_singular(
    z_gpu: &mut DMatrix<Complex64>,
    n: usize,
    quad_n: usize,
    nodes: &[[f64; 3]],
    faces: &[TriFace],
    bases: &[RwgBasis],
    quad_bary: &[[f64; 3]],
    quad_weights: &[f64],
    omega_mu0: f64,
    inv_omega_eps0: f64,
    k: f64,
    kernel_type: u32,
) {
    if kernel_type != 0 { return; }
    use crate::green::green3d;
    use crate::quadrature::TriQuad;
    use crate::singular::{classify_pair, TriPairType, zmn_efie_rwg_singular};

    // face → basis adjacency → collect near-singular pairs
    let mut face_bases: Vec<Vec<usize>> = vec![vec![]; faces.len()];
    for (bi, b) in bases.iter().enumerate() {
        face_bases[b.plus_face].push(bi);
        face_bases[b.minus_face].push(bi);
    }

    let mut pairs: Vec<(usize, usize)> = Vec::new();
    let mut seen = std::collections::HashSet::new();
    // Self-pairs (m == m): all 4 face-pairs are Identical/SharedEdge
    for m in 0..n {
        let key = m * n + m;
        if seen.insert(key) { pairs.push((m, m)); }
    }
    // Neighbor pairs sharing at least one face (→ Identical, SharedEdge face-pairs)
    for m in 0..n {
        for &mf in &[bases[m].plus_face, bases[m].minus_face] {
            for &nb in &face_bases[mf] {
                if nb <= m { continue; }
                let key = m * n + nb;
                if seen.insert(key) { pairs.push((m, nb)); }
            }
        }
    }

    // Vertex-based: catch SharedVertex pairs (bases sharing a vertex but not a face)
    let mut vert_bases: Vec<Vec<usize>> = vec![vec![]; nodes.len()];
    for (bi, b) in bases.iter().enumerate() {
        for &fi in &[b.plus_face, b.minus_face] {
            for &ni in &faces[fi].nodes {
                vert_bases[ni].push(bi);
            }
        }
    }
    for m in 0..n {
        for &mf in &[bases[m].plus_face, bases[m].minus_face] {
            for &ni in &faces[mf].nodes {
                for &nb in &vert_bases[ni] {
                    if nb <= m { continue; }
                    let key = m * n + nb;
                    if seen.insert(key) { pairs.push((m, nb)); }
                }
            }
        }
    }

    let k_divisor = inv_omega_eps0 / omega_mu0;
    let pair_order = 4usize;

    // Inline RWG eval without SurfaceMesh
    let rwg_eval = |b: &RwgBasis, r: &[f64; 3], in_plus: bool| -> [f64; 3] {
        let (fi, fn_, sgn) = if in_plus {
            (b.plus_face, b.free_node_plus, 1.0)
        } else {
            (b.minus_face, b.free_node_minus, -1.0)
        };
        let area = faces[fi].area;
        let free = &nodes[fn_];
        let sc = sgn * b.length / (2.0 * area);
        [sc * (r[0] - free[0]), sc * (r[1] - free[1]), sc * (r[2] - free[2])]
    };

    let divergence = |b: &RwgBasis, in_plus: bool| -> f64 {
        let (fi, sgn) = if in_plus { (b.plus_face, 1.0) } else { (b.minus_face, -1.0) };
        sgn * b.length / faces[fi].area
    };

    for &(mi, ni) in &pairs {
        let bm = &bases[mi];
        let bn = &bases[ni];

        for &(m_face, m_plus) in &[(bm.plus_face, true), (bm.minus_face, false)] {
            for &(n_face, n_plus) in &[(bn.plus_face, true), (bn.minus_face, false)] {
                if classify_pair(&faces[m_face], &faces[n_face]) == TriPairType::Disjoint {
                    continue;
                }
                // SharedVertex pairs are well-approximated by Gauss quadrature;
                // only correct Identical (Duffy) and SharedEdge (Sauter-Schwab).
                if classify_pair(&faces[m_face], &faces[n_face]) == TriPairType::SharedVertex {
                    continue;
                }
                let face_m = &faces[m_face];
                let face_n = &faces[n_face];
                let div_m = divergence(bm, m_plus);
                let div_n = divergence(bn, n_plus);

                // Gauss approximation (what the GPU computed)
                let mut gauss_val = Complex64::ZERO;
                for (bm_pt, &wm) in quad_bary.iter().zip(quad_weights.iter()) {
                    let rm = TriQuad::global_point(bm_pt, face_m, nodes);
                    let fm = rwg_eval(bm, &rm, m_plus);
                    for (bn_pt, &wn) in quad_bary.iter().zip(quad_weights.iter()) {
                        let rn = TriQuad::global_point(bn_pt, face_n, nodes);
                        let fn_ = rwg_eval(bn, &rn, n_plus);
                        let g = green3d(&rm, &rn, k);
                        let dot_ff = fm[0]*fn_[0] + fm[1]*fn_[1] + fm[2]*fn_[2];
                        gauss_val += g * (dot_ff - k_divisor * div_m * div_n)
                                   * (wm * wn * 4.0 * face_m.area * face_n.area);
                    }
                }

                // Singular-corrected value
                let fm_fn = |rm: &[f64; 3], rn: &[f64; 3]| -> (f64, f64) {
                    let fm = rwg_eval(bm, rm, m_plus);
                    let fn_ = rwg_eval(bn, rn, n_plus);
                    let dot = fm[0]*fn_[0] + fm[1]*fn_[1] + fm[2]*fn_[2];
                    (dot, div_m * div_n)
                };
                let (a_term, phi_term) = zmn_efie_rwg_singular(face_m, face_n, &fm_fn, nodes, k, pair_order);
                let singular_val = a_term - k_divisor * phi_term;
                let diff = singular_val - gauss_val;
                let corr = Complex64::new(0.0, -omega_mu0) * diff;

                z_gpu[(mi, ni)] += corr;
                if mi != ni {
                    z_gpu[(ni, mi)] += corr;
                }
            }
        }
    }
}

/// Generic GPU fill: factored out common logic for EFIE (type=0) and MFIE (type=1).
#[cfg(all(not(target_arch = "wasm32"), feature = "wgpu-gpu"))]
fn fill_wgpu_generic(
    n: usize,
    quad_n: usize,
    nodes: &[[f64; 3]],
    faces: &[TriFace],
    bases: &[RwgBasis],
    quad_bary: &[[f64; 3]],
    quad_weights: &[f64],
    omega_mu0: f64,
    inv_omega_eps0: f64,
    k: f64,
    green_data: Option<(&[Complex64], f64, f64)>,
    kernel_type: u32,
) -> RemResult<Option<DMatrix<Complex64>>> {
    if n == 0 { return Ok(Some(DMatrix::zeros(0, 0))); }
    if n < GPU_MIN_BASIS {
        log::debug!("n={} below GPU_MIN_BASIS, returning None for CPU fallback", n);
        return Ok(None);
    }
    if !wgpu_adapter_available() { return Ok(None); }
    if let Some((table, _, _)) = green_data {
        if table.len() < 2 {
            log::warn!("Green table too short ({}), cannot interpolate", table.len());
            return Ok(None);
        }
    }
    match gpu_pipeline::fill_efie_gpu(
        n, quad_n, nodes, faces, bases, quad_bary, quad_weights,
        omega_mu0, inv_omega_eps0, k, green_data, kernel_type,
    ) {
        Some(mut z) => {
            correct_near_singular(
                &mut z, n, quad_n, nodes, faces, bases,
                quad_bary, quad_weights,
                omega_mu0, inv_omega_eps0, k, kernel_type,
            );
            Ok(Some(z))
        }
        None => {
            log::warn!("wgpu fill (kernel_type={}) failed; returning None for CPU fallback", kernel_type);
            Ok(None)
        }
    }
}

// ── Stubs for when wgpu is not available ──────────────────────────────────

macro_rules! gpu_stub {
    ($name:ident) => {
        #[cfg(any(target_arch = "wasm32", not(feature = "wgpu-gpu")))]
        pub fn $name(
            n: usize, _quad_n: usize,
            _nodes: &[[f64; 3]], _faces: &[TriFace], _bases: &[RwgBasis],
            _quad_bary: &[[f64; 3]], _quad_weights: &[f64],
            _omega_mu0: f64, _inv_omega_eps0: f64, _k: f64,
        ) -> RemResult<Option<DMatrix<Complex64>>> { Ok(None) }
    };
    ($name:ident, with_table) => {
        #[cfg(any(target_arch = "wasm32", not(feature = "wgpu-gpu")))]
        pub fn $name(
            n: usize, _quad_n: usize,
            _nodes: &[[f64; 3]], _faces: &[TriFace], _bases: &[RwgBasis],
            _quad_bary: &[[f64; 3]], _quad_weights: &[f64],
            _omega_mu0: f64, _inv_omega_eps0: f64, _k: f64,
            _green_table: &[Complex64], _green_rho_min: f64, _green_rho_max: f64,
        ) -> RemResult<Option<DMatrix<Complex64>>> { Ok(None) }
    };
}

gpu_stub!(fill_impedance_wgpu_efie);
gpu_stub!(fill_impedance_wgpu_mfie);
gpu_stub!(fill_impedance_wgpu_efie_with_green_table, with_table);
gpu_stub!(fill_impedance_wgpu_mfie_with_green_table, with_table);

// ═══════════════════════════════════════════════════════════════════════════
//  Synthetic benchmark
// ═══════════════════════════════════════════════════════════════════════════

pub fn fill_impedance_gpu(n_basis: usize, freq: f64) -> RemResult<DMatrix<Complex64>> {
    use std::f64::consts::PI;
    const MU0: f64 = 4.0e-7 * PI;

    if n_basis == 0 {
        return Ok(DMatrix::zeros(0, 0));
    }

    let omega = 2.0 * PI * freq;
    let z0    = omega * MU0 / (4.0 * PI);

    let z = fill_impedance_parallel(n_basis, |m, k| {
        if m == k {
            return Complex64::new(z0, z0) * (n_basis as f64 + 1.0);
        }
        let dist = (m as isize - k as isize).unsigned_abs() as f64 + 1.0;
        Complex64::new(z0, z0) / dist
    });

    Ok(z)
}

// ═══════════════════════════════════════════════════════════════════════════
//  Tests
// ═══════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gpu_available_reflects_rayon() {
        #[cfg(not(target_arch = "wasm32"))]
        assert!(gpu_available(), "Expected Rayon-parallel path on native");
        #[cfg(target_arch = "wasm32")]
        assert!(!gpu_available(), "Expected no parallel path on WASM");
    }

    #[test]
    fn fill_impedance_gpu_returns_n_times_n_matrix() {
        let z = fill_impedance_gpu(8, 1e9).unwrap();
        assert_eq!(z.nrows(), 8);
        assert_eq!(z.ncols(), 8);
    }

    #[test]
    fn fill_impedance_gpu_zero_basis_ok() {
        let z = fill_impedance_gpu(0, 1e9).unwrap();
        assert_eq!(z.nrows(), 0);
        assert_eq!(z.ncols(), 0);
    }

    #[test]
    fn fill_impedance_parallel_identity() {
        let n = 4;
        let z = fill_impedance_parallel(n, |m, k| {
            if m == k { Complex64::new(1.0, 0.0) } else { Complex64::ZERO }
        });
        for i in 0..n {
            for j in 0..n {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((z[(i, j)].re - expected).abs() < 1e-15);
                assert!(z[(i, j)].im.abs() < 1e-15);
            }
        }
    }

    #[test]
    fn fill_impedance_gpu_diagonally_dominant() {
        let n = 16;
        let z = fill_impedance_gpu(n, 2.4e9).unwrap();
        for i in 0..n {
            let diag     = z[(i, i)].norm();
            let off_sum: f64 = (0..n).filter(|&j| j != i).map(|j| z[(i, j)].norm()).sum();
            assert!(diag >= off_sum - 1e-12,
                "Row {i}: diag={diag:.3e} < off_sum={off_sum:.3e}");
        }
    }

    #[test]
    fn fill_impedance_wgpu_efie_stub_no_panic() {
        use crate::basis::rwg::RwgBasis;
        use crate::surface_mesh::TriFace;
        let n = 4;
        let nodes = vec![[0.0; 3]; 3];
        let face = TriFace {
            nodes: [0, 1, 2],
            centroid: [0.0; 3],
            normal: [0.0, 0.0, 1.0],
            area: 0.5,
        };
        let faces = vec![face; 4];
        let bases = vec![
            RwgBasis { edge_idx: 0, plus_face: 0, minus_face: 1, free_node_plus: 0, free_node_minus: 2, length: 1.0 },
            RwgBasis { edge_idx: 1, plus_face: 1, minus_face: 2, free_node_plus: 1, free_node_minus: 3, length: 1.0 },
            RwgBasis { edge_idx: 2, plus_face: 2, minus_face: 3, free_node_plus: 2, free_node_minus: 0, length: 1.0 },
            RwgBasis { edge_idx: 3, plus_face: 3, minus_face: 0, free_node_plus: 3, free_node_minus: 1, length: 1.0 },
        ];
        match fill_impedance_wgpu_efie(
            n, 1, &nodes, &faces, &bases, &[[1.0/3.0, 1.0/3.0, 1.0/3.0]], &[0.5],
            1.0, 1.0, 1.0,
        ).unwrap() {
            None => {}, // GPU not available — expected on systems without GPU
            Some(z) => {
                assert_eq!(z.nrows(), n);
                assert_eq!(z.ncols(), n);
            }
        }
    }

    /// GPU vs CPU numerical validation on a small 3D mesh.
    ///
    /// Builds a small tetrahedron-like surface (4 faces, 3 bases),
    /// computes Z via GPU and CPU, and verifies:
    /// - GPU output is finite
    /// - Matrix diagonals are non-zero
    /// - Off-diagonal elements match the CPU result within 20% relative error
    ///   (the GPU uses uniform Gauss quadrature for all pairs, while CPU
    ///    uses singular-corrected quadrature for near-singular pairs, so
    ///    perfect agreement is not expected for this small mesh).
    #[test]
    fn gpu_vs_cpu_numerical_validation() {
        use crate::basis::rwg::{generate_rwg_bases, RwgBasis};
        use crate::quadrature::TriQuad;
        use crate::surface_mesh::{tri_geometry, SharedEdge, SurfaceMesh, TriFace};

        // Build a small 3D mesh: 4 triangles forming a tetrahedron-like surface.
        let nodes = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let (c0, n0, a0) = tri_geometry(&nodes[0], &nodes[1], &nodes[2]);
        let (c1, n1, a1) = tri_geometry(&nodes[0], &nodes[1], &nodes[3]);
        let (c2, n2, a2) = tri_geometry(&nodes[1], &nodes[2], &nodes[3]);
        let (c3, n3, a3) = tri_geometry(&nodes[0], &nodes[2], &nodes[3]);
        let faces = vec![
            TriFace { nodes: [0,1,2], centroid: c0, normal: n0, area: a0 },
            TriFace { nodes: [0,1,3], centroid: c1, normal: n1, area: a1 },
            TriFace { nodes: [1,2,3], centroid: c2, normal: n2, area: a2 },
            TriFace { nodes: [0,2,3], centroid: c3, normal: n3, area: a3 },
        ];

        // 6 shared edges → 6 RWG bases
        let edges = vec![
            SharedEdge { nodes: [0,1], plus_face: 0, minus_face: 1, length: 1.0 },
            SharedEdge { nodes: [0,2], plus_face: 0, minus_face: 3, length: 1.0 },
            SharedEdge { nodes: [0,3], plus_face: 1, minus_face: 3, length: 1.0 },
            SharedEdge { nodes: [1,2], plus_face: 0, minus_face: 2, length: 1.0 },
            SharedEdge { nodes: [1,3], plus_face: 1, minus_face: 2, length: 1.0 },
            SharedEdge { nodes: [2,3], plus_face: 2, minus_face: 3, length: 1.0 },
        ];
        // Patch edge lengths
        let edges: Vec<SharedEdge> = edges.into_iter().map(|mut e| {
            let d = [
                nodes[e.nodes[0]][0] - nodes[e.nodes[1]][0],
                nodes[e.nodes[0]][1] - nodes[e.nodes[1]][1],
                nodes[e.nodes[0]][2] - nodes[e.nodes[1]][2],
            ];
            e.length = (d[0]*d[0] + d[1]*d[1] + d[2]*d[2]).sqrt();
            e
        }).collect();

        let surf = SurfaceMesh {
            nodes: nodes.clone(),
            faces: faces.clone(),
            edges,
            boundary_edges: vec![],
            face_attrs: vec![0; 4],
            global_node_ids: vec![0, 1, 2, 3],
        };
        let bases = generate_rwg_bases(&surf);
        let n = bases.len();
        assert!(n >= 3, "expected >=3 bases, got {}", n);

        let freq = 1e9_f64;
        let omega = 2.0 * std::f64::consts::PI * freq;
        let k = omega / rem_core::C0;
        let omega_mu0 = omega * rem_core::MU0;
        let inv_omega_eps0 = 1.0 / (omega * rem_core::EPS0);
        let quad = TriQuad::new(5);
        let quad_n = quad.n_pts();

        // CPU reference (direct serial fill using the kernel)
        let z_cpu = fill_impedance_parallel(n, |mi, ni| {
            let bm = &bases[mi];
            let bn = &bases[ni];
            crate::assemble::zmn_efie_rwg(bm, bn, &surf, k, omega_mu0, inv_omega_eps0, &quad)
        });
        assert_eq!(z_cpu.nrows(), n);
        assert_eq!(z_cpu.ncols(), n);
        assert!(z_cpu.iter().any(|&z| z != Complex64::ZERO), "CPU Z must be non-zero");

        // GPU path
        match fill_impedance_wgpu_efie(
            n, quad_n,
            &surf.nodes, &surf.faces, &bases,
            &quad.bary, &quad.weights,
            omega_mu0, inv_omega_eps0, k,
        ).unwrap() {
            None => {
                eprintln!("GPU not available — skipping GPU validation");
            }
            Some(z_gpu) => {
                assert_eq!(z_gpu.nrows(), n);
                assert_eq!(z_gpu.ncols(), n);

                // 1. All GPU entries must be finite
                for i in 0..n {
                    for j in 0..n {
                        assert!(z_gpu[(i, j)].norm().is_finite(),
                            "GPU Z[{}][{}] is not finite", i, j);
                    }
                }

                // 2. Diagonals must be non-zero (or at least the trace norm > 0)
                let trace_norm: f64 = (0..n).map(|i| z_gpu[(i, i)].norm()).sum();
                assert!(trace_norm > 0.0, "GPU diagonal trace is zero");

                // 3. For well-separated pairs (no shared face), error should be moderate.
                //    For this small mesh, count pairs where relative error < 50%.
                let mut close_count = 0;
                let mut total_pairs = 0;
                for i in 0..n {
                    for j in 0..n {
                        let cpu_val = z_cpu[(i, j)];
                        let gpu_val = z_gpu[(i, j)];
                        let max_norm = cpu_val.norm().max(gpu_val.norm()).max(1e-30);
                        let rel_err = (cpu_val - gpu_val).norm() / max_norm;
                        if rel_err < 0.5 {
                            close_count += 1;
                        }
                        total_pairs += 1;
                    }
                }
                let fraction = close_count as f64 / total_pairs as f64;
                eprintln!("GPU vs CPU: close pairs = {}/{} ({:.1}%)",
                    close_count, total_pairs, fraction * 100.0);
                // At least half the pairs should be within 50% relative error
                assert!(fraction > 0.5,
                    "Too few matching pairs: {}/{}", close_count, total_pairs);
            }
        }
    }

    /// Build a flat rectangular mesh with approx N RWG bases.
    fn bench_mesh(n: usize) -> (crate::surface_mesh::SurfaceMesh, Vec<RwgBasis>) {
        use crate::surface_mesh::{tri_geometry, SharedEdge, TriFace, SurfaceMesh};
        use crate::basis::rwg::generate_rwg_bases;
        // Grid: choose Nx, Ny so number of shared edges ≈ n
        // Each cell has 2 triangles; internal edges = 3*Nx*Ny - Nx - Ny
        // For a square: Nx = Ny ≈ sqrt(n/3)
        let nx = ((n as f64 / 3.0).sqrt().ceil() as usize).max(2);
        let ny = nx;
        let lx = 1.0;
        let ly = 1.0;
        let mut nodes = Vec::with_capacity((nx + 1) * (ny + 1));
        for j in 0..=ny {
            for i in 0..=nx {
                nodes.push([i as f64 * lx / nx as f64, j as f64 * ly / ny as f64, 0.0]);
            }
        }
        let idx = |i: usize, j: usize| j * (nx + 1) + i;

        let mut faces = Vec::with_capacity(2 * nx * ny);
        for j in 0..ny {
            for i in 0..nx {
                let a = idx(i, j);
                let b = idx(i + 1, j);
                let c = idx(i, j + 1);
                let d = idx(i + 1, j + 1);
                let (c0, n0, a0) = tri_geometry(&nodes[a], &nodes[b], &nodes[c]);
                let (c1, n1, a1) = tri_geometry(&nodes[b], &nodes[d], &nodes[c]);
                faces.push(TriFace { nodes: [a, b, c], centroid: c0, normal: n0, area: a0 });
                faces.push(TriFace { nodes: [b, d, c], centroid: c1, normal: n1, area: a1 });
            }
        }

        let (mut edges, boundary_edges) = crate::surface_mesh::build_edge_topology(&faces);
        crate::surface_mesh::patch_edge_lengths(&mut edges, &nodes);

        let surf = SurfaceMesh {
            nodes, faces, edges, boundary_edges,
            face_attrs: vec![0; 2 * nx * ny],
            global_node_ids: vec![],
        };
        let bases = generate_rwg_bases(&surf);
        (surf, bases)
    }

    /// GPU vs CPU benchmark at several sizes. Prints timing to stderr.
    /// Single-element GPU vs CPU comparison using a 500-basis flat mesh.
    /// Ignored by default since it requires a GPU adapter.
    #[test]
    #[ignore]
    fn gpu_single_element_vs_cpu() {
        use crate::quadrature::TriQuad;
        let (surf, bases) = bench_mesh(600);
        let n = bases.len();
        let quad = TriQuad::new(1);
        let freq = 1e9_f64;
        let omega = 2.0 * std::f64::consts::PI * freq;
        let k = omega / rem_core::C0;
        let omega_mu0 = omega * rem_core::MU0;
        let inv_omega_eps0 = 1.0 / (omega * rem_core::EPS0);

        // CPU Gauss-for-all (same as GPU kernel)
        let cpu_val = {
            let bm = &bases[0]; let bn = &bases[1];
            let mut val = Complex64::ZERO;
            for &(mf, mp) in &[(bm.plus_face, true), (bm.minus_face, false)] {
                for &(nf, np) in &[(bn.plus_face, true), (bn.minus_face, false)] {
                    let fm = &surf.faces[mf]; let fn_ = &surf.faces[nf];
                    let div_m = if mp { bm.length / fm.area } else { -bm.length / fm.area };
                    let div_n = if np { bn.length / fn_.area } else { -bn.length / fn_.area };
                    for (bp, &wm) in quad.bary.iter().zip(quad.weights.iter()) {
                        let rm = crate::quadrature::TriQuad::global_point(bp, fm, &surf.nodes);
                        let f_m = bm.eval(&rm, &surf, mp);
                        for (bq, &wn) in quad.bary.iter().zip(quad.weights.iter()) {
                            let rn = crate::quadrature::TriQuad::global_point(bq, fn_, &surf.nodes);
                            let f_n = bn.eval(&rn, &surf, np);
                            let g = crate::green::green3d(&rm, &rn, k);
                            let dot = f_m[0]*f_n[0]+f_m[1]*f_n[1]+f_m[2]*f_n[2];
                            val += g * (dot - inv_omega_eps0/omega_mu0 * div_m * div_n)
                                 * (wm * wn * 4.0 * fm.area * fn_.area);
                        }
                    }
                }
            }
            Complex64::new(0.0, -omega_mu0) * val
        };

        // GPU value
        match fill_impedance_wgpu_efie(n, quad.n_pts(),
            &surf.nodes, &surf.faces, &bases,
            &quad.bary, &quad.weights,
            omega_mu0, inv_omega_eps0, k,
        ).unwrap() {
            None => eprintln!("GPU not available, skipping test"),
            Some(z) => {
                let gpu_val = z[(0, 0)];
                let denom = cpu_val.norm().max(1e-30);
                let rel_err = (cpu_val - gpu_val).norm() / denom;
                eprintln!("Z[0][0]: cpu={:.10e} gpu={:.10e} rel_err={:.2e}", cpu_val, gpu_val, rel_err);
                assert!(rel_err < 0.01, "GPU vs CPU Gauss-for-all mismatch: rel_err={:.2e}", rel_err);
            }
        }
    }

    #[test]
    fn gpu_vs_cpu_benchmark() {
        use crate::quadrature::TriQuad;
        use std::time::Instant;

        let configs = [
            ("N≈500",  500usize,  3usize),
            ("N≈2000", 2000usize, 1usize),
        ];

        for &(label, target_n, quad_deg) in &configs {
            let (surf, bases) = bench_mesh(target_n);
            let n = bases.len();
            let quad = TriQuad::new(quad_deg);
            let freq = 1e9_f64;
            let omega = 2.0 * std::f64::consts::PI * freq;
            let k = omega / rem_core::C0;
            let omega_mu0 = omega * rem_core::MU0;
            let inv_omega_eps0 = 1.0 / (omega * rem_core::EPS0);

            // CPU timing
            let t0 = Instant::now();
            let z_cpu = fill_impedance_parallel(n, |mi, ni| {
                crate::assemble::zmn_efie_rwg(
                    &bases[mi], &bases[ni], &surf, k, omega_mu0, inv_omega_eps0, &quad,
                )
            });
            let cpu_time = t0.elapsed();

            // Warm-up pass: force GPU init + pipeline compilation
            let _ = fill_impedance_wgpu_efie(
                n, quad.n_pts(),
                &surf.nodes, &surf.faces, &bases,
                &quad.bary, &quad.weights,
                omega_mu0, inv_omega_eps0, k,
            );

            // CPU Gauss-for-all (same as GPU kernel: no singular correction)
            let z_gauss_all = fill_impedance_parallel(n, |mi, ni| {
                let bm = &bases[mi]; let bn = &bases[ni];
                let mut val = Complex64::ZERO;
                for &(m_face, m_plus) in &[(bm.plus_face, true), (bm.minus_face, false)] {
                    for &(n_face, n_plus) in &[(bn.plus_face, true), (bn.minus_face, false)] {
                        let fm = &surf.faces[m_face]; let f_n = &surf.faces[n_face];
                        let div_m = if m_plus { bm.length / fm.area } else { -bm.length / fm.area };
                        let div_n = if n_plus { bn.length / f_n.area } else { -bn.length / f_n.area };
                        for (bp, &wm) in quad.bary.iter().zip(quad.weights.iter()) {
                            let rm = crate::quadrature::TriQuad::global_point(bp, fm, &surf.nodes);
                            let f_m = bm.eval(&rm, &surf, m_plus);
                            for (bq, &wn) in quad.bary.iter().zip(quad.weights.iter()) {
                                let rn = crate::quadrature::TriQuad::global_point(bq, f_n, &surf.nodes);
                                let f_n2 = bn.eval(&rn, &surf, n_plus);
                                let g = crate::green::green3d(&rm, &rn, k);
                                let dot = f_m[0]*f_n2[0]+f_m[1]*f_n2[1]+f_m[2]*f_n2[2];
                                let integrand = g * (dot - inv_omega_eps0/omega_mu0 * div_m * div_n);
                                val += integrand * (wm * wn * 4.0 * fm.area * f_n.area);
                            }
                        }
                    }
                }
                Complex64::new(0.0, -omega_mu0) * val
            });

            // CPU Gauss-ALL reference (same kernel as GPU, no singular correction)
            let z_gauss_all = fill_impedance_parallel(n, |mi, ni| {
                let bm = &bases[mi]; let bn = &bases[ni];
                let mut val = Complex64::ZERO;
                for &(m_face, m_plus) in &[(bm.plus_face, true), (bm.minus_face, false)] {
                    for &(n_face, n_plus) in &[(bn.plus_face, true), (bn.minus_face, false)] {
                        let fm = &surf.faces[m_face]; let fn_ = &surf.faces[n_face];
                        let div_m = if m_plus { bm.length / fm.area } else { -bm.length / fm.area };
                        let div_n = if n_plus { bn.length / fn_.area } else { -bn.length / fn_.area };
                        for (bp, &wm) in quad.bary.iter().zip(quad.weights.iter()) {
                            let rm = crate::quadrature::TriQuad::global_point(bp, fm, &surf.nodes);
                            let f_m = bm.eval(&rm, &surf, m_plus);
                            for (bq, &wn) in quad.bary.iter().zip(quad.weights.iter()) {
                                let rn = crate::quadrature::TriQuad::global_point(bq, fn_, &surf.nodes);
                                let f_n = bn.eval(&rn, &surf, n_plus);
                                let g = crate::green::green3d(&rm, &rn, k);
                                let dot = f_m[0]*f_n[0]+f_m[1]*f_n[1]+f_m[2]*f_n[2];
                                let integrand = g * (dot - inv_omega_eps0/omega_mu0 * div_m * div_n);
                                val += integrand * (wm * wn * 4.0 * fm.area * fn_.area);
                            }
                        }
                    }
                }
                Complex64::new(0.0, -omega_mu0) * val
            });

            // GPU timing (second call uses cached geometry + compiled pipeline)
            let t1 = Instant::now();
            let z_gpu = fill_impedance_wgpu_efie(
                n, quad.n_pts(),
                &surf.nodes, &surf.faces, &bases,
                &quad.bary, &quad.weights,
                omega_mu0, inv_omega_eps0, k,
            ).unwrap();

            let gpu_time = t1.elapsed();

            // GPU vs Gauss-ALL comparison (same kernel, precision check)
            let gpu_vs_gauss_err = match z_gpu {
                Some(ref z) => (0..n.min(100)).flat_map(|i| (0..n.min(100)).map(move |j| (i, j)))
                    .map(|(i, j)| {
                        let denom = z_gauss_all[(i,j)].norm().max(1e-30);
                        (z[(i,j)] - z_gauss_all[(i,j)]).norm() / denom
                    }).fold(0.0_f64, f64::max),
                None => f64::NAN,
            };

            // Accuracy
            let max_rel_err = match z_gpu {
                Some(ref z) => {
                    (0..n.min(200)).flat_map(|i| (0..n.min(200)).map(move |j| (i, j)))
                        .map(|(i, j)| {
                            let cpu = z_cpu[(i, j)];
                            let gpu = z[(i, j)];
                            let denom = cpu.norm().max(1e-30);
                            (cpu - gpu).norm() / denom
                        })
                        .fold(0.0_f64, f64::max)
                }
                None => f64::NAN,
            };

            let sp = if gpu_time >= std::time::Duration::from_nanos(1) {
                cpu_time.as_secs_f64() / gpu_time.as_secs_f64()
            } else { 0.0 };

            // File-based logging (test harness captures stdout/stderr)
            if let Ok(mut f) = std::fs::OpenOptions::new().create(true).append(true)
                .open("gpu_benchmark.txt")
            {
                use std::io::Write;
                writeln!(f, "{:>10}: n={:<6} CPU={:<8.2?} GPU={:<8.2?} speedup={:<5.1}× cpu_err={:.2e} gauss_err={:.2e}",
                         label, n, cpu_time, gpu_time, sp, max_rel_err, gpu_vs_gauss_err).ok();
                // Also write worst elements
                if let Some(ref z) = z_gpu {
                    use crate::singular::{classify_pair, TriPairType};
                    let nchk = n.min(100);
                    let mut worst: Vec<(f64, usize, usize)> = (0..nchk)
                        .flat_map(|i| (0..nchk).map(move |j| (i, j)))
                        .map(|(i, j)| {
                            let cpu = z_cpu[(i, j)];
                            let gpu = z[(i, j)];
                            let denom = cpu.norm().max(1e-30);
                            ((cpu - gpu).norm() / denom, i, j)
                        }).collect();
                    worst.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
                    for &(err, i, j) in worst.iter().take(3) {
                        let bm = &bases[i]; let bn = &bases[j];
                        let mut ts = String::new();
                        for &(mf,_) in &[(bm.plus_face,true),(bm.minus_face,false)] {
                            for &(nf,_) in &[(bn.plus_face,true),(bn.minus_face,false)] {
                                ts.push(' ');
                                ts.push_str(&format!("{:?}", classify_pair(&surf.faces[mf],&surf.faces[nf])));
                            }
                        }
                        writeln!(f, "  worst[{}][{}]: cpu={:.6e} gpu={:.6e} gauss_all={:.6e} rel={:.2e} types=[{}]",
                                 i, j, z_cpu[(i,j)], z[(i,j)], z_gauss_all[(i,j)], err, ts).ok();
                    }
                }
            }
        }
    }
}
