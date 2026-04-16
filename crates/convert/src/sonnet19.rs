use anyhow::{bail, Context};
use quick_xml::events::Event;
use quick_xml::Reader;
use rmsh_io::{save_msh_v2_to_path, save_step_to_path};
use rmsh_model::{Element, ElementType, Mesh, Node};
use serde_json::{json, Value};
use std::collections::BTreeSet;
use std::path::Path;

const DEFAULT_FREQ_MIN: f64 = 1.0e9;
const DEFAULT_FREQ_MAX: f64 = 10.0e9;
const DEFAULT_FREQ_STEP: f64 = 0.1e9;
const DEFAULT_CELLS_X: usize = 80;
const DEFAULT_CELLS_Y: usize = 40;
const DEFAULT_WIDTH_M: f64 = 10.0e-3;
const DEFAULT_HEIGHT_M: f64 = 5.0e-3;
const MAX_GENERATED_CELLS_PRODUCT: usize = 120_000;
const BASE_PEC_TAG: u32 = 1;
const BOTTOM_TAG: u32 = 2;
const TOP_TAG: u32 = 3;
const LEFT_TAG: u32 = 4;
const RIGHT_TAG: u32 = 5;
const PORT_TAG_BASE: u32 = 1001;

#[derive(Debug, Clone)]
struct SonnetDielectricLayer {
    name: String,
    eps_r: f64,           // isotropic permittivity
    loss_tan: f64,        // loss tangent
    mu_r: f64,            // permeability
    thickness_m: f64,     // layer thickness [m]
}

#[derive(Debug, Clone)]
struct SonnetHints {
    width_m: Option<f64>,
    height_m: Option<f64>,
    cells_x: Option<usize>,
    cells_y: Option<usize>,
    freq_min_hz: Option<f64>,
    freq_max_hz: Option<f64>,
    freq_step_hz: Option<f64>,
    ports: Vec<SonnetPort>,
    dielectric_layers: Vec<SonnetDielectricLayer>,
    y_direction_negative: bool,  // true if YDirection="Negative"
    local_origin_y_m: Option<f64>,
    output_folder: Option<String>,
    matrix_solver: Option<String>,
    sweep_type: Option<String>,
    speed_control: Option<String>,
    precision_mode: Option<String>,
    deembed_on: Option<bool>,
    subs_per_lambda: Option<f64>,
    ref_planes: Vec<SonnetRefPlane>,
    conductor_polygons_m: Vec<Vec<(f64, f64)>>,
}

#[derive(Debug, Clone)]
struct SonnetRefPlane {
    side: String,
    plane_type: String,
    cal_length_m: Option<f64>,
}

#[derive(Debug, Clone, Default)]
struct SonnetPort {
    number: Option<i32>,
    resistance_ohm: Option<f64>,
    center_xy_m: Option<(f64, f64)>,
    vertex_1based: Option<usize>,
    direction_hint: Option<String>,
    parent_polygon_points_m: Vec<(f64, f64)>,
}

#[derive(Debug, Clone)]
struct PendingPort {
    port: SonnetPort,
    attr: u32,
}

#[derive(Debug, Clone, Default)]
struct PolygonCapture {
    points_m: Vec<(f64, f64)>,
}

#[derive(Debug, Clone, Copy)]
struct SweepRaw {
    start: f64,
    stop: f64,
    step: Option<f64>,
    target: Option<f64>,
}

#[derive(Debug, Clone, Copy)]
struct ParserState {
    length_unit_scale: f64,
    freq_unit_scale: f64,
}

impl Default for ParserState {
    fn default() -> Self {
        Self {
            // Sonnet files commonly specify UM and GHZ explicitly; these are safe defaults.
            length_unit_scale: 1.0e-6,
            freq_unit_scale: 1.0e9,
        }
    }
}

impl SonnetHints {
    fn new() -> Self {
        Self {
            width_m: None,
            height_m: None,
            cells_x: None,
            cells_y: None,
            freq_min_hz: None,
            freq_max_hz: None,
            freq_step_hz: None,
            ports: Vec::new(),
            dielectric_layers: Vec::new(),
            y_direction_negative: false,
            local_origin_y_m: None,
            output_folder: None,
            matrix_solver: None,
            sweep_type: None,
            speed_control: None,
            precision_mode: None,
            deembed_on: None,
            subs_per_lambda: None,
            ref_planes: Vec::new(),
            conductor_polygons_m: Vec::new(),
        }
    }
}

pub fn convert_xml_to_rem(
    xml_path: &Path,
    out_config: &Path,
    out_msh: &Path,
    freq_min_override: Option<f64>,
    freq_max_override: Option<f64>,
    freq_step_override: Option<f64>,
    debug_step_out: Option<&Path>,
) -> anyhow::Result<()> {
    let xml = std::fs::read_to_string(xml_path)
        .with_context(|| format!("reading Sonnet XML: {}", xml_path.display()))?;

    let hints = parse_sonnet_hints(&xml)?;

    let width_m = hints.width_m.unwrap_or(DEFAULT_WIDTH_M);
    let height_m = hints.height_m.unwrap_or(DEFAULT_HEIGHT_M);
    let cells_x_raw = hints.cells_x.unwrap_or(DEFAULT_CELLS_X).max(2);
    let cells_y_raw = hints.cells_y.unwrap_or(DEFAULT_CELLS_Y).max(2);
    let (cells_x, cells_y) = cap_cells(cells_x_raw, cells_y_raw);

    let freq_min = freq_min_override.or(hints.freq_min_hz).unwrap_or(DEFAULT_FREQ_MIN);
    let freq_max = freq_max_override.or(hints.freq_max_hz).unwrap_or(DEFAULT_FREQ_MAX);
    let freq_step = freq_step_override.or(hints.freq_step_hz).unwrap_or(DEFAULT_FREQ_STEP);

    if !(width_m.is_finite() && width_m > 0.0 && height_m.is_finite() && height_m > 0.0) {
        bail!("invalid Sonnet geometry extent; width/height must be positive");
    }
    if !(freq_min > 0.0 && freq_max >= freq_min && freq_step > 0.0) {
        bail!("invalid frequency sweep; require freq_min > 0, freq_max >= freq_min, freq_step > 0");
    }

    if let Some(parent) = out_msh.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("creating output directory: {}", parent.display()))?;
    }
    if let Some(parent) = out_config.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("creating output directory: {}", parent.display()))?;
    }
    if let Some(step_path) = debug_step_out {
        if let Some(parent) = step_path.parent() {
            std::fs::create_dir_all(parent)
                .with_context(|| format!("creating debug STEP directory: {}", parent.display()))?;
        }
    }

    if let Some(step_path) = debug_step_out {
        let geom_mesh = build_debug_geometry_surface_mesh(
            width_m,
            height_m,
            hints.y_direction_negative,
            hints.local_origin_y_m.unwrap_or(height_m),
            &hints.conductor_polygons_m,
        );
        save_step_to_path(step_path, &geom_mesh)
            .with_context(|| format!("writing debug STEP: {}", step_path.display()))?;
    }

    let (mesh, port_tags, pec_tags) = generate_rect_msh_with_port_tags(
        width_m,
        height_m,
        cells_x,
        cells_y,
        &hints.ports,
        hints.y_direction_negative,
        hints.local_origin_y_m.unwrap_or(height_m),
    );
    save_msh_v2_to_path(out_msh, &mesh)
        .with_context(|| format!("writing mesh: {}", out_msh.display()))?;

    let mut pending_ports: Vec<PendingPort> = hints
        .ports
        .iter()
        .enumerate()
        .map(|(i, p)| PendingPort {
            port: p.clone(),
            attr: port_tags.get(i).copied().unwrap_or(BASE_PEC_TAG),
        })
        .collect();
    pending_ports.sort_by_key(|pp| pp.port.number.unwrap_or(i32::MAX));

    let mut ports_json: Vec<Value> = Vec::new();
    for (i, pp) in pending_ports.iter().enumerate() {
        let p = &pp.port;
        let index = p
            .number
            .and_then(|n| if n > 0 { Some(n as u32) } else { None })
            .unwrap_or((i + 1) as u32);
        let attr = pp.attr;
        let direction = p
            .direction_hint
            .clone()
            .or_else(|| infer_direction_from_polygon_vertex(&p.parent_polygon_points_m, p.vertex_1based))
            .unwrap_or_else(|| infer_port_direction(p.center_xy_m, width_m, height_m).to_string());
        let z0 = p.resistance_ohm.unwrap_or(50.0);
        ports_json.push(json!({
            "Index": index,
            "Attributes": [attr],
            "Direction": direction,
            "Impedance": z0
        }));
    }

    let fast_solver = map_fast_solver(hints.matrix_solver.as_deref(), hints.speed_control.as_deref());
    let singular_tol = map_singular_tol(hints.precision_mode.as_deref(), hints.speed_control.as_deref());

    let mut mom = serde_json::Map::new();
    mom.insert("Equation".to_string(), json!("CFIE"));
    mom.insert("Basis".to_string(), json!("RWG"));
    mom.insert("FreqMin".to_string(), json!(freq_min));
    mom.insert("FreqMax".to_string(), json!(freq_max));
    mom.insert("FreqStep".to_string(), json!(freq_step));
    mom.insert("Alpha".to_string(), json!(0.5));
    mom.insert("SingularTol".to_string(), json!(singular_tol));
    mom.insert("FastSolver".to_string(), json!(fast_solver));
    mom.insert("RefImpedance".to_string(), json!(50.0));
    if let Some(deembed) = hints.deembed_on {
        mom.insert("Deembed".to_string(), json!(deembed));
    }
    if !hints.ref_planes.is_empty() {
        let refs: Vec<Value> = hints
            .ref_planes
            .iter()
            .map(|rp| {
                json!({
                    "Side": rp.side,
                    "Type": rp.plane_type,
                    "CalLength": rp.cal_length_m,
                })
            })
            .collect();
        mom.insert("RefPlanes".to_string(), Value::Array(refs));
    }
    if !ports_json.is_empty() {
        mom.insert("Ports".to_string(), Value::Array(ports_json));
    }
    
    // Add substrate layers if present (for future layered-green GreenFunction support)
    if !hints.dielectric_layers.is_empty() {
        let mut substrate = serde_json::Map::new();
        let mut layers_json: Vec<Value> = Vec::new();
        for layer in &hints.dielectric_layers {
            layers_json.push(json!({
                "Permittivity": layer.eps_r,
                "LossTangent": layer.loss_tan,
                "Permeability": layer.mu_r,
                "Thickness": layer.thickness_m,
                "Name": layer.name
            }));
        }
        substrate.insert("Layers".to_string(), Value::Array(layers_json));
        substrate.insert("BottomPec".to_string(), json!(true));
        mom.insert("Substrate".to_string(), Value::Object(substrate));
    }

    let problem_output = hints
        .output_folder
        .as_deref()
        .filter(|s| !s.trim().is_empty())
        .map(|s| if s == "." { "./output/sonnet19".to_string() } else { s.to_string() })
        .unwrap_or_else(|| "./output/sonnet19".to_string());

    let cfg = json!({
        "Problem": {
            "Type": "MoM",
            "Verbose": 1,
            "Output": problem_output
        },
        "Model": {
            "Mesh": out_msh.to_string_lossy(),
            "L0": 1.0
        },
        "Boundaries": {
            "PEC": { "Attributes": pec_tags }
        },
        "Solver": {
            "MoM": mom
        }
    });

    let pretty = serde_json::to_string_pretty(&cfg)?;
    std::fs::write(out_config, pretty)
        .with_context(|| format!("writing config: {}", out_config.display()))?;

    log::info!(
        "Sonnet19 XML converted: {} -> config={}, mesh={} (w={:.6e} m, h={:.6e} m, nx={} (raw {}), ny={} (raw {}), fast_solver={}, singular_tol={:.1e})",
        xml_path.display(),
        out_config.display(),
        out_msh.display(),
        width_m,
        height_m,
        cells_x,
        cells_x_raw,
        cells_y,
        cells_y_raw,
        fast_solver,
        singular_tol,
    );

    Ok(())
}

fn build_debug_geometry_surface_mesh(
    width_m: f64,
    height_m: f64,
    y_direction_negative: bool,
    local_origin_y_m: f64,
    conductor_polygons_m: &[Vec<(f64, f64)>],
) -> Mesh {
    let mut mesh = Mesh::new();

    let mut next_node_id = 1u64;
    let mut next_elem_id = 1u64;

    for poly in conductor_polygons_m {
        let mut nodes: Vec<u64> = Vec::new();
        for &(x, y_raw) in poly {
            let y = if y_direction_negative {
                local_origin_y_m - y_raw
            } else {
                y_raw
            };
            let nid = next_node_id;
            next_node_id += 1;
            mesh.add_node(Node::new(nid, x, y, 0.0));
            nodes.push(nid);
        }

        // Triangulate each polygon with a simple fan around the first vertex.
        for tri_i in 1..nodes.len().saturating_sub(1) {
            let mut tri = Element::new(
                next_elem_id,
                ElementType::Triangle3,
                vec![nodes[0], nodes[tri_i], nodes[tri_i + 1]],
            );
            tri.physical_tag = Some(BASE_PEC_TAG as i32);
            mesh.add_element(tri);
            next_elem_id += 1;
        }
    }

    if mesh.elements.is_empty() {
        let y0 = if y_direction_negative {
            local_origin_y_m
        } else {
            0.0
        };
        let y1 = if y_direction_negative {
            local_origin_y_m - height_m
        } else {
            height_m
        };

        mesh.add_node(Node::new(1, 0.0, y0, 0.0));
        mesh.add_node(Node::new(2, width_m, y0, 0.0));
        mesh.add_node(Node::new(3, width_m, y1, 0.0));
        mesh.add_node(Node::new(4, 0.0, y1, 0.0));

        let mut t1 = Element::new(1, ElementType::Triangle3, vec![1, 2, 3]);
        t1.physical_tag = Some(BASE_PEC_TAG as i32);
        mesh.add_element(t1);

        let mut t2 = Element::new(2, ElementType::Triangle3, vec![1, 3, 4]);
        t2.physical_tag = Some(BASE_PEC_TAG as i32);
        mesh.add_element(t2);
    }

    mesh
}

fn map_fast_solver(matrix_solver: Option<&str>, speed_control: Option<&str>) -> &'static str {
    let ms = matrix_solver.unwrap_or("").to_ascii_uppercase();
    if ms.contains("DIRECT") {
        return "Direct";
    }
    if ms.contains("ITER") || ms.contains("GMRES") {
        return "GMRES";
    }
    if ms.contains("AUTO") {
        return "ACA";
    }

    let sc = speed_control.unwrap_or("").to_ascii_uppercase();
    if sc.contains("MAX_SPEED") {
        "GMRES"
    } else {
        "ACA"
    }
}

fn map_singular_tol(precision_mode: Option<&str>, speed_control: Option<&str>) -> f64 {
    let p = precision_mode.unwrap_or("").to_ascii_uppercase();
    if p.contains("QUAD") {
        return 1.0e-7;
    }
    if p.contains("DOUBLE") {
        return 1.0e-6;
    }
    if p.contains("SINGLE") || p.contains("MIN") {
        return 1.0e-5;
    }
    let sc = speed_control.unwrap_or("").to_ascii_uppercase();
    if sc.contains("MAX_ACCURACY") {
        1.0e-6
    } else if sc.contains("MAX_SPEED") {
        5.0e-5
    } else {
        1.0e-5
    }
}

fn cap_cells(nx: usize, ny: usize) -> (usize, usize) {
    let prod = nx.saturating_mul(ny);
    if prod <= MAX_GENERATED_CELLS_PRODUCT {
        return (nx, ny);
    }

    let scale = (MAX_GENERATED_CELLS_PRODUCT as f64 / prod as f64).sqrt();
    let nx2 = ((nx as f64) * scale).round().max(2.0) as usize;
    let ny2 = ((ny as f64) * scale).round().max(2.0) as usize;
    (nx2, ny2)
}

fn infer_port_direction(center_xy_m: Option<(f64, f64)>, width_m: f64, height_m: f64) -> &'static str {
    if let Some((x, y)) = center_xy_m {
        let dx = x.min((width_m - x).max(0.0));
        let dy = y.min((height_m - y).max(0.0));
        if dx <= dy {
            "x"
        } else {
            "y"
        }
    } else {
        "x"
    }
}

fn generate_rect_msh_with_port_tags(
    w: f64,
    h: f64,
    nx: usize,
    ny: usize,
    ports: &[SonnetPort],
    y_direction_negative: bool,
    local_origin_y_m: f64,
) -> (Mesh, Vec<u32>, Vec<u32>) {
    let mut port_tags: Vec<u32> = Vec::with_capacity(ports.len());
    for i in 0..ports.len() {
        port_tags.push(PORT_TAG_BASE + i as u32);
    }

    let dx = w / nx as f64;
    let band_half = (2.0 * dx).max(w * 0.03);
    let mut col_tags = vec![BASE_PEC_TAG; nx];

    for ix in 0..nx {
        let cx = (ix as f64 + 0.5) * dx;
        let mut best: Option<(usize, f64)> = None;
        for (pi, p) in ports.iter().enumerate() {
            if let Some((px, _)) = p.center_xy_m {
                let d = (cx - px).abs();
                if d <= band_half {
                    match best {
                        Some((_, bd)) if d >= bd => {}
                        _ => best = Some((pi, d)),
                    }
                }
            }
        }
        if let Some((pi, _)) = best {
            col_tags[ix] = port_tags[pi];
        }
    }

    let mut pec_set: BTreeSet<u32> = BTreeSet::new();
    pec_set.insert(BASE_PEC_TAG);
    for t in &col_tags {
        pec_set.insert(*t);
    }
    let pec_tags: Vec<u32> = pec_set.into_iter().collect();

    let node_id = |ix: usize, iy: usize| -> u64 { (iy * (nx + 1) + ix + 1) as u64 };

    let mut mesh = Mesh::new();

    mesh.physical_names.insert((2, BASE_PEC_TAG as i32), "conductor".to_string());
    mesh.physical_names.insert((1, BOTTOM_TAG as i32), "bottom".to_string());
    mesh.physical_names.insert((1, TOP_TAG as i32), "top".to_string());
    mesh.physical_names.insert((1, LEFT_TAG as i32), "left".to_string());
    mesh.physical_names.insert((1, RIGHT_TAG as i32), "right".to_string());
    for (i, tag) in port_tags.iter().enumerate() {
        mesh.physical_names
            .insert((2, *tag as i32), format!("port_{}", i + 1));
    }

    for iy in 0..=ny {
        for ix in 0..=nx {
            let id = node_id(ix, iy);
            let x = ix as f64 * w / nx as f64;
            let y_raw = iy as f64 * h / ny as f64;
            // Apply YDirection flip if needed: y_actual = local_origin_y - y_raw
            let y = if y_direction_negative {
                local_origin_y_m - y_raw
            } else {
                y_raw
            };
            mesh.add_node(Node::new(id, x, y, 0.0));
        }
    }

    let mut eid = 1u64;

    for ix in 0..nx {
        let mut e = Element::new(
            eid,
            ElementType::Line2,
            vec![node_id(ix, 0), node_id(ix + 1, 0)],
        );
        e.physical_tag = Some(BOTTOM_TAG as i32);
        mesh.add_element(e);
        eid += 1;
    }
    for ix in 0..nx {
        let mut e = Element::new(
            eid,
            ElementType::Line2,
            vec![node_id(ix, ny), node_id(ix + 1, ny)],
        );
        e.physical_tag = Some(TOP_TAG as i32);
        mesh.add_element(e);
        eid += 1;
    }
    for iy in 0..ny {
        let mut e = Element::new(
            eid,
            ElementType::Line2,
            vec![node_id(0, iy), node_id(0, iy + 1)],
        );
        e.physical_tag = Some(LEFT_TAG as i32);
        mesh.add_element(e);
        eid += 1;
    }
    for iy in 0..ny {
        let mut e = Element::new(
            eid,
            ElementType::Line2,
            vec![node_id(nx, iy), node_id(nx, iy + 1)],
        );
        e.physical_tag = Some(RIGHT_TAG as i32);
        mesh.add_element(e);
        eid += 1;
    }

    for iy in 0..ny {
        for ix in 0..nx {
            let tag = col_tags[ix];
            let n00 = node_id(ix, iy);
            let n10 = node_id(ix + 1, iy);
            let n01 = node_id(ix, iy + 1);
            let n11 = node_id(ix + 1, iy + 1);

            let mut e1 = Element::new(eid, ElementType::Triangle3, vec![n00, n10, n11]);
            e1.physical_tag = Some(tag as i32);
            mesh.add_element(e1);
            eid += 1;

            let mut e2 = Element::new(eid, ElementType::Triangle3, vec![n00, n11, n01]);
            e2.physical_tag = Some(tag as i32);
            mesh.add_element(e2);
            eid += 1;
        }
    }

    (mesh, port_tags, pec_tags)
}

fn parse_sonnet_hints(xml: &str) -> anyhow::Result<SonnetHints> {
    let mut reader = Reader::from_str(xml);
    reader.config_mut().trim_text(true);

    let mut buf = Vec::new();
    let mut hints = SonnetHints::new();
    let mut state = ParserState::default();
    let mut path: Vec<String> = Vec::new();
    let mut pending_sweep: Option<SweepRaw> = None;
    let mut current_port: Option<SonnetPort> = None;
    let mut current_polygon: Option<PolygonCapture> = None;
    let mut current_level_material: Option<String> = None;

    loop {
        match reader.read_event_into(&mut buf) {
            Ok(Event::Start(e)) => {
                let tag = local_name(e.name().as_ref());
                path.push(tag.clone());
                parse_start_like(
                    &tag,
                    &path,
                    e.attributes(),
                    &mut hints,
                    &mut state,
                    &mut pending_sweep,
                    &mut current_port,
                    &mut current_polygon,
                    &mut current_level_material,
                );
            }
            Ok(Event::Empty(e)) => {
                let tag = local_name(e.name().as_ref());
                path.push(tag.clone());
                parse_start_like(
                    &tag,
                    &path,
                    e.attributes(),
                    &mut hints,
                    &mut state,
                    &mut pending_sweep,
                    &mut current_port,
                    &mut current_polygon,
                    &mut current_level_material,
                );
                path.pop();
            }
            Ok(Event::Text(t)) => {
                let text = String::from_utf8_lossy(t.as_ref()).trim().to_string();
                if !text.is_empty() {
                    parse_text_node(&path, &text, &mut state);
                    if let Some(last) = path.last() {
                        map_kv_hint(&mut hints, last, &text, state);
                        map_text_hint(&mut hints, &path, last, &text);
                    }
                    if path_contains(path.as_slice(), "port") {
                        if let Some(last) = path.last() {
                            if last == "center" {
                                if let Some((x, y)) = parse_center_xy(&text, state.length_unit_scale) {
                                    if let Some(p) = current_port.as_mut() {
                                        p.center_xy_m = Some((x, y));
                                    }
                                }
                            }
                        }
                    }
                    if path_contains(path.as_slice(), "planarpolygon") {
                        if let Some(last) = path.last() {
                            if last == "points" {
                                if let Some(poly) = current_polygon.as_mut() {
                                    poly.points_m.extend(parse_points_list(&text, state.length_unit_scale));
                                }
                            }
                        }
                    }
                }
            }
            Ok(Event::End(e)) => {
                let tag = local_name(e.name().as_ref());
                if tag == "port" {
                    if let Some(p) = current_port.take() {
                        hints.ports.push(p);
                    }
                } else if tag == "planarpolygon" {
                    if let Some(poly) = current_polygon.take() {
                        if poly.points_m.len() >= 3 {
                            hints.conductor_polygons_m.push(poly.points_m);
                        }
                    }
                } else if tag == "level" {
                    current_level_material = None;
                }
                path.pop();
            }
            Ok(Event::Eof) => break,
            Ok(_) => {}
            Err(e) => {
                bail!("failed to parse Sonnet XML: {}", e);
            }
        }
        buf.clear();
    }

    if let Some(sw) = pending_sweep {
        let start_hz = sw.start * state.freq_unit_scale;
        let stop_hz = sw.stop * state.freq_unit_scale;
        let step_hz = if let Some(step) = sw.step {
            step * state.freq_unit_scale
        } else if let Some(target) = sw.target {
            // Sonnet target can be very large; cap generated sweep points for usability.
            let raw_step = ((stop_hz - start_hz).abs() / target.max(1.0)).max(1.0);
            let capped_step = ((stop_hz - start_hz).abs() / 400.0).max(1.0);
            raw_step.max(capped_step)
        } else {
            DEFAULT_FREQ_STEP
        };

        hints.freq_min_hz = Some(start_hz.min(stop_hz));
        hints.freq_max_hz = Some(start_hz.max(stop_hz));
        hints.freq_step_hz = Some(step_hz);
    }

    Ok(hints)
}

fn parse_start_like(
    tag: &str,
    path: &[String],
    attrs: quick_xml::events::attributes::Attributes,
    hints: &mut SonnetHints,
    state: &mut ParserState,
    pending_sweep: &mut Option<SweepRaw>,
    current_port: &mut Option<SonnetPort>,
    current_polygon: &mut Option<PolygonCapture>,
    current_level_material: &mut Option<String>,
) {
    let in_geometry_box = path_contains(path, "geometry") && path_contains(path, "box");

    let mut attrs_vec: Vec<(String, String)> = Vec::new();
    for attr in attrs.flatten() {
        let key = local_name(attr.key.as_ref());
        let value = String::from_utf8_lossy(attr.value.as_ref()).to_string();
        attrs_vec.push((key, value));
    }

    // Parse <Geometry YDirection="Negative"> root attribute
    if tag == "geometry" {
        for (k, v) in &attrs_vec {
            if k == "ydirection" && v.to_ascii_lowercase() == "negative" {
                hints.y_direction_negative = true;
            }
        }
    }

    // Parse <LocalOrigin X="0" Y="1200"/> within Box
    if tag == "localorigin" && in_geometry_box {
        for (k, v) in &attrs_vec {
            if k == "y" {
                hints.local_origin_y_m = parse_f64(v).map(|n| n * state.length_unit_scale);
            }
        }
    }

    if tag == "outputfiles" {
        for (k, v) in &attrs_vec {
            if k == "folder" {
                hints.output_folder = Some(v.clone());
            }
        }
    }

    if tag == "deembed" {
        for (k, v) in &attrs_vec {
            if k == "on" {
                let vv = v.trim().to_ascii_lowercase();
                hints.deembed_on = Some(vv == "true" || vv == "1");
            }
        }
    }

    if tag == "refplane" && in_geometry_box {
        let mut side: Option<String> = None;
        let mut plane_type: Option<String> = None;
        let mut cal_length_m: Option<f64> = None;
        for (k, v) in &attrs_vec {
            if k == "side" {
                side = Some(v.trim().to_string());
            } else if k == "type" {
                plane_type = Some(v.trim().to_string());
            } else if k == "callength" {
                cal_length_m = parse_f64(v).map(|n| n * state.length_unit_scale);
            }
        }
        if side.is_some() || plane_type.is_some() || cal_length_m.is_some() {
            hints.ref_planes.push(SonnetRefPlane {
                side: side.unwrap_or_else(|| "UNKNOWN".to_string()),
                plane_type: plane_type.unwrap_or_else(|| "UNKNOWN".to_string()),
                cal_length_m,
            });
        }
    }

    // Parse <Dielectric MacroID="..." Name="..."> with nested Eps/Tan/MuRel
    if tag == "dielectric" && path_contains(path, "geometry") {
        let mut dielectric = SonnetDielectricLayer {
            name: String::new(),
            eps_r: 1.0,
            loss_tan: 0.0,
            mu_r: 1.0,
            thickness_m: 0.0,
        };
        for (k, v) in &attrs_vec {
            if k == "name" {
                dielectric.name = v.clone();
            }
        }
        hints.dielectric_layers.push(dielectric);
    }

    // Update dielectric layer parameters (Eps, Tan, MuRel have Value attribute)
    if path_contains(path, "geometry") && path_contains(path, "dielectric") {
        if !hints.dielectric_layers.is_empty() {
            let diel = hints.dielectric_layers.last_mut().unwrap();
            if tag == "eps" {
                for (k, v) in &attrs_vec {
                    if k == "value" {
                        if let Some(val) = parse_f64(v) {
                            diel.eps_r = val;
                        }
                    }
                }
            } else if tag == "tan" {
                for (k, v) in &attrs_vec {
                    if k == "value" {
                        if let Some(val) = parse_f64(v) {
                            diel.loss_tan = val;
                        }
                    }
                }
            } else if tag == "murel" {
                for (k, v) in &attrs_vec {
                    if k == "value" {
                        if let Some(val) = parse_f64(v) {
                            diel.mu_r = val;
                        }
                    }
                }
            }
        }
    }

    // Capture current level material so thickness can be mapped to the correct dielectric.
    if tag == "level" && path_contains(path, "geometry") {
        for (k, v) in &attrs_vec {
            if k == "materialname" {
                *current_level_material = Some(v.clone());
            }
        }
    }

    // Parse <DielectricMaterialModel Thickness="10.0"/> for layer thickness
    if tag == "dielectricmaterialmodel" && path_contains(path, "level") && path_contains(path, "geometry") {
        for (k, v) in &attrs_vec {
            if k == "thickness" {
                if let Some(val) = parse_f64(v) {
                    let thickness = val * state.length_unit_scale;
                    if let Some(level_mat) = current_level_material.as_ref() {
                        if let Some(layer) = hints
                            .dielectric_layers
                            .iter_mut()
                            .find(|d| d.name.eq_ignore_ascii_case(level_mat))
                        {
                            layer.thickness_m = thickness;
                            break;
                        }
                    }
                    if let Some(last) = hints.dielectric_layers.last_mut() {
                        if last.thickness_m == 0.0 {
                            last.thickness_m = thickness;
                        }
                    }
                }
            }
        }
    }

    if tag == "size" && in_geometry_box {
        let mut sx = None;
        let mut sy = None;
        for (k, v) in &attrs_vec {
            if k == "x" {
                sx = parse_f64(v).map(|n| n * state.length_unit_scale);
            } else if k == "y" {
                sy = parse_f64(v).map(|n| n * state.length_unit_scale);
            }
        }
        if let Some(v) = sx {
            hints.width_m = Some(v);
        }
        if let Some(v) = sy {
            hints.height_m = Some(v);
        }
    }

    if tag == "numcells" && in_geometry_box {
        for (k, v) in &attrs_vec {
            if k == "x" {
                hints.cells_x = parse_usize(v);
            } else if k == "y" {
                hints.cells_y = parse_usize(v);
            }
        }
    }

    if tag == "port" {
        let mut port = SonnetPort::default();
        if let Some(poly) = current_polygon.as_ref() {
            port.parent_polygon_points_m = poly.points_m.clone();
        }
        for (k, v) in &attrs_vec {
            if k == "number" {
                port.number = v.trim().parse::<i32>().ok();
            } else if k == "vertex" {
                port.vertex_1based = parse_usize(v);
            }
        }
        *current_port = Some(port);
    }

    if tag == "planarpolygon" {
        *current_polygon = Some(PolygonCapture::default());
    }

    if tag == "impedance" {
        if let Some(p) = current_port.as_mut() {
            for (k, v) in &attrs_vec {
                if k == "resistance" {
                    p.resistance_ohm = parse_f64(v);
                }
            }
        }
    }

    if tag == "gndref" {
        if let Some(p) = current_port.as_mut() {
            for (k, v) in &attrs_vec {
                if k == "direction" {
                    let dir = v.trim().to_ascii_lowercase();
                    if dir == "x" || dir == "y" || dir == "z" {
                        p.direction_hint = Some(dir);
                    }
                }
            }
        }
    }

    if tag == "sweep" && path_contains(path, "frequencies") {
        let mut start = None;
        let mut stop = None;
        let mut step = None;
        let mut target = None;
        for (k, v) in &attrs_vec {
            if k == "start" {
                start = parse_f64(v);
            } else if k == "stop" {
                stop = parse_f64(v);
            } else if k == "step" {
                step = parse_f64(v);
            } else if k == "target" {
                target = parse_f64(v);
            }
        }
        if let (Some(s), Some(e)) = (start, stop) {
            *pending_sweep = Some(SweepRaw {
                start: s,
                stop: e,
                step,
                target,
            });
        }
    }

    for (k, v) in attrs_vec {
        map_kv_hint(hints, &k, &v, *state);
    }
}

fn parse_center_xy(raw: &str, scale: f64) -> Option<(f64, f64)> {
    let t = raw.trim();
    let inner = t.strip_prefix('(')?.strip_suffix(')')?;
    let mut parts = inner.split(',').map(|s| s.trim());
    let x = parts.next()?.parse::<f64>().ok()? * scale;
    let y = parts.next()?.parse::<f64>().ok()? * scale;
    Some((x, y))
}

fn parse_points_list(raw: &str, scale: f64) -> Vec<(f64, f64)> {
    let mut pts = Vec::new();
    let mut rest = raw;
    loop {
        let Some(start) = rest.find('(') else { break };
        let rem = &rest[start + 1..];
        let Some(end) = rem.find(')') else { break };
        let pair = &rem[..end];
        let mut it = pair.split(',').map(|s| s.trim());
        if let (Some(xs), Some(ys)) = (it.next(), it.next()) {
            if let (Ok(x), Ok(y)) = (xs.parse::<f64>(), ys.parse::<f64>()) {
                pts.push((x * scale, y * scale));
            }
        }
        rest = &rem[end + 1..];
    }
    pts
}

fn infer_direction_from_polygon_vertex(points: &[(f64, f64)], vertex_1based: Option<usize>) -> Option<String> {
    if points.len() < 2 {
        return None;
    }
    let vi = vertex_1based?.saturating_sub(1);
    if vi >= points.len() {
        return None;
    }
    let n = points.len();
    let prev = points[(vi + n - 1) % n];
    let next = points[(vi + 1) % n];
    let tx = next.0 - prev.0;
    let ty = next.1 - prev.1;
    let dir = if tx.abs() >= ty.abs() { "y" } else { "x" };
    Some(dir.to_string())
}

fn parse_text_node(path: &[String], text: &str, state: &mut ParserState) {
    let Some(last) = path.last() else { return; };
    if path_contains(path, "units") && last == "length" {
        if let Some(scale) = length_unit_scale(text) {
            state.length_unit_scale = scale;
        }
    }
    if path_contains(path, "units") && last == "frequency" {
        if let Some(scale) = freq_unit_scale(text) {
            state.freq_unit_scale = scale;
        }
    }
}

fn map_text_hint(hints: &mut SonnetHints, path: &[String], last: &str, text: &str) {
    let in_control = path_contains(path, "control");
    let in_ufft = path_contains(path, "ufftcontrol");

    if in_control && !in_ufft && last == "sweeptype" {
        hints.sweep_type = Some(text.trim().to_string());
    }
    if in_control && !in_ufft && last == "speedcontrol" {
        hints.speed_control = Some(text.trim().to_string());
    }
    if in_control && !in_ufft && last == "precision" {
        hints.precision_mode = Some(text.trim().to_string());
    }
    if in_ufft && last == "matrixsolver" {
        hints.matrix_solver = Some(text.trim().to_string());
    }
    if in_control && !in_ufft && last == "subsperlambda" {
        hints.subs_per_lambda = parse_f64(text);
    }
}

fn map_kv_hint(hints: &mut SonnetHints, key: &str, raw: &str, state: ParserState) {
    let k = key.to_ascii_lowercase();

    if is_x_extent_key(&k) {
        if let Some(v) = parse_length_m(raw, state.length_unit_scale) {
            hints.width_m = Some(v);
        }
    } else if is_y_extent_key(&k) {
        if let Some(v) = parse_length_m(raw, state.length_unit_scale) {
            hints.height_m = Some(v);
        }
    } else if is_cells_x_key(&k) {
        if let Some(v) = parse_usize(raw) {
            hints.cells_x = Some(v);
        }
    } else if is_cells_y_key(&k) {
        if let Some(v) = parse_usize(raw) {
            hints.cells_y = Some(v);
        }
    } else if is_freq_min_key(&k) {
        if let Some(v) = parse_freq_hz(raw, state.freq_unit_scale) {
            hints.freq_min_hz = Some(v);
        }
    } else if is_freq_max_key(&k) {
        if let Some(v) = parse_freq_hz(raw, state.freq_unit_scale) {
            hints.freq_max_hz = Some(v);
        }
    } else if is_freq_step_key(&k) {
        if let Some(v) = parse_freq_hz(raw, state.freq_unit_scale) {
            hints.freq_step_hz = Some(v);
        }
    }
}

fn parse_usize(raw: &str) -> Option<usize> {
    let s = raw.trim();
    if s.is_empty() {
        return None;
    }
    if let Ok(v) = s.parse::<usize>() {
        return Some(v);
    }
    let int_part: String = s.chars().take_while(|c| c.is_ascii_digit()).collect();
    if int_part.is_empty() {
        None
    } else {
        int_part.parse::<usize>().ok()
    }
}

fn parse_freq_hz(raw: &str, default_scale: f64) -> Option<f64> {
    parse_scaled(raw, UnitKind::Frequency, default_scale)
}

fn parse_length_m(raw: &str, default_scale: f64) -> Option<f64> {
    parse_scaled(raw, UnitKind::Length, default_scale)
}

enum UnitKind {
    Frequency,
    Length,
}

fn parse_scaled(raw: &str, kind: UnitKind, default_scale: f64) -> Option<f64> {
    let s = raw.trim().to_ascii_lowercase();
    if s.is_empty() {
        return None;
    }

    let mut split = 0usize;
    for (i, c) in s.char_indices() {
        if !(c.is_ascii_digit() || c == '.' || c == '-' || c == '+' || c == 'e') {
            split = i;
            break;
        }
    }

    let (num_s, unit_s) = if split == 0 {
        (&s[..], "")
    } else {
        (&s[..split], s[split..].trim())
    };

    let num = num_s.parse::<f64>().ok()?;
    if !num.is_finite() {
        return None;
    }

    let scale = match kind {
        UnitKind::Frequency => match unit_s {
            "" => default_scale,
            "hz" => 1.0,
            "khz" | "k" => 1.0e3,
            "mhz" | "m" => 1.0e6,
            "ghz" | "g" => 1.0e9,
            _ => default_scale,
        },
        UnitKind::Length => match unit_s {
            "" => default_scale,
            "m" => 1.0,
            "mm" => 1.0e-3,
            "um" | "micron" | "microns" => 1.0e-6,
            "nm" => 1.0e-9,
            "cm" => 1.0e-2,
            "mil" => 25.4e-6,
            "in" | "inch" | "inches" => 0.0254,
            _ => default_scale,
        },
    };

    Some(num * scale)
}

fn parse_f64(raw: &str) -> Option<f64> {
    raw.trim().parse::<f64>().ok()
}

fn length_unit_scale(unit: &str) -> Option<f64> {
    match unit.trim().to_ascii_lowercase().as_str() {
        "m" => Some(1.0),
        "cm" => Some(1.0e-2),
        "mm" => Some(1.0e-3),
        "um" => Some(1.0e-6),
        "nm" => Some(1.0e-9),
        "mil" => Some(25.4e-6),
        "in" | "inch" => Some(0.0254),
        _ => None,
    }
}

fn freq_unit_scale(unit: &str) -> Option<f64> {
    match unit.trim().to_ascii_lowercase().as_str() {
        "hz" => Some(1.0),
        "khz" => Some(1.0e3),
        "mhz" => Some(1.0e6),
        "ghz" => Some(1.0e9),
        _ => None,
    }
}

fn local_name(bytes: &[u8]) -> String {
    let raw = String::from_utf8_lossy(bytes).to_string();
    raw.rsplit(':').next().unwrap_or(&raw).to_ascii_lowercase()
}

fn path_contains(path: &[String], seg: &str) -> bool {
    path.iter().any(|p| p == seg)
}

fn is_x_extent_key(k: &str) -> bool {
    matches_any(k, &[
        "xsize", "sizex", "boxx", "boxsizex", "xwidth", "widthx", "xlength", "lengthx", "xdim", "xextent",
    ])
}

fn is_y_extent_key(k: &str) -> bool {
    matches_any(k, &[
        "ysize", "sizey", "boxy", "boxsizey", "ywidth", "widthy", "ylength", "lengthy", "ydim", "yextent",
    ])
}

fn is_cells_x_key(k: &str) -> bool {
    matches_any(k, &["xcells", "cellsx", "xcell", "cellsx", "nx", "nxcells"])
}

fn is_cells_y_key(k: &str) -> bool {
    matches_any(k, &["ycells", "cellsy", "ycell", "ny", "nycells"])
}

fn is_freq_min_key(k: &str) -> bool {
    matches_any(k, &["freqmin", "fmin", "startfreq", "startfrequency", "beginfreq"])
}

fn is_freq_max_key(k: &str) -> bool {
    matches_any(k, &["freqmax", "fmax", "stopfreq", "stopfrequency", "endfreq"])
}

fn is_freq_step_key(k: &str) -> bool {
    matches_any(k, &["freqstep", "df", "stepfreq", "frequencystep"])
}

fn matches_any(hay: &str, keys: &[&str]) -> bool {
    keys.iter().any(|k| hay == *k || hay.contains(k))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_hints_from_typical_sonnet_like_xml() {
        let xml = r#"
<Project>
  <Box XSize="20mm" YSize="10mm" XCells="120" YCells="60" />
  <Sweep FreqMin="1GHz" FreqMax="5GHz" FreqStep="100MHz" />
</Project>
"#;
        let hints = parse_sonnet_hints(xml).unwrap();
        assert_eq!(hints.cells_x, Some(120));
        assert_eq!(hints.cells_y, Some(60));
        assert!((hints.width_m.unwrap() - 0.02).abs() < 1.0e-12);
        assert!((hints.height_m.unwrap() - 0.01).abs() < 1.0e-12);
        assert!((hints.freq_min_hz.unwrap() - 1.0e9).abs() < 1.0);
        assert!((hints.freq_max_hz.unwrap() - 5.0e9).abs() < 1.0);
        assert!((hints.freq_step_hz.unwrap() - 1.0e8).abs() < 1.0);
    }

    #[test]
    fn parse_text_nodes_and_units() {
        let xml = r#"
<Project>
  <XSize>12.5 mm</XSize>
  <YSize>6.25mm</YSize>
  <FreqMin>2.4GHz</FreqMin>
  <FreqMax>2.8GHz</FreqMax>
  <FreqStep>50MHz</FreqStep>
</Project>
"#;
        let hints = parse_sonnet_hints(xml).unwrap();
        assert!((hints.width_m.unwrap() - 0.0125).abs() < 1.0e-12);
        assert!((hints.height_m.unwrap() - 0.00625).abs() < 1.0e-12);
        assert!((hints.freq_min_hz.unwrap() - 2.4e9).abs() < 1.0);
        assert!((hints.freq_max_hz.unwrap() - 2.8e9).abs() < 1.0);
        assert!((hints.freq_step_hz.unwrap() - 50.0e6).abs() < 1.0);
    }

    #[test]
    fn parse_readout_fixture() {
        let xml = include_str!("../../../tests/fixtures/readout.sonx");
        let hints = parse_sonnet_hints(xml).unwrap();
        assert_eq!(hints.cells_x, Some(1600));
        assert_eq!(hints.cells_y, Some(1200));
        assert!((hints.width_m.unwrap() - 1600.0e-6).abs() < 1.0e-12);
        assert!((hints.height_m.unwrap() - 1200.0e-6).abs() < 1.0e-12);
        assert!((hints.freq_min_hz.unwrap() - 0.1e9).abs() < 1.0);
        assert!((hints.freq_max_hz.unwrap() - 10.0e9).abs() < 1.0);
        assert!(hints.freq_step_hz.unwrap() > 0.0);
        assert_eq!(hints.ports.len(), 2);
    }

    #[test]
    fn parse_xy_bit_fixture() {
        let xml = include_str!("../../../tests/fixtures/xy_bit.sonx");
        let hints = parse_sonnet_hints(xml).unwrap();
        assert_eq!(hints.cells_x, Some(800));
        assert_eq!(hints.cells_y, Some(800));
        assert!((hints.width_m.unwrap() - 400.0e-6).abs() < 1.0e-12);
        assert!((hints.height_m.unwrap() - 400.0e-6).abs() < 1.0e-12);
        assert!((hints.freq_min_hz.unwrap() - 4.0e9).abs() < 1.0);
        assert!((hints.freq_max_hz.unwrap() - 8.0e9).abs() < 1.0);
        assert_eq!(hints.ports.len(), 5);
    }

    #[test]
    fn parse_interdigital_fixture() {
        let xml = include_str!("../../../tests/fixtures/interdigital_capacitor.sonx");
        let hints = parse_sonnet_hints(xml).unwrap();
        assert_eq!(hints.cells_x, Some(1536));
        assert_eq!(hints.cells_y, Some(1728));
        assert!((hints.width_m.unwrap() - 1536.0e-6).abs() < 1.0e-12);
        assert!((hints.height_m.unwrap() - 1728.0e-6).abs() < 1.0e-12);
        assert!((hints.freq_min_hz.unwrap() - 0.1e9).abs() < 1.0);
        assert!((hints.freq_max_hz.unwrap() - 1.0e9).abs() < 1.0);
        assert!((hints.freq_step_hz.unwrap() - 0.001e9).abs() < 1.0);
        assert_eq!(hints.ports.len(), 2);
    }

    #[test]
    fn parse_supercond_dielectrics_and_levels() {
        let xml = include_str!("../../../testdata/sonnet/supercond_filter/supercond_filter_actualhousing.sonx");
        let hints = parse_sonnet_hints(xml).unwrap();

        let air = hints
            .dielectric_layers
            .iter()
            .find(|d| d.name == "Air")
            .expect("Air dielectric should be parsed");
        assert!((air.eps_r - 1.0).abs() < 1.0e-12);
        assert!((air.mu_r - 1.0).abs() < 1.0e-12);
        assert!((air.thickness_m - 3.81e-3).abs() < 1.0e-12);

        let mgo = hints
            .dielectric_layers
            .iter()
            .find(|d| d.name == "MgO")
            .expect("MgO dielectric should be parsed");
        assert!((mgo.eps_r - 9.7).abs() < 1.0e-12);
        assert!((mgo.mu_r - 1.0).abs() < 1.0e-12);
        assert!((mgo.thickness_m - 0.507e-3).abs() < 1.0e-12);
    }

    #[test]
    fn parse_supercond_conductor_polygons() {
        let xml = include_str!("../../../testdata/sonnet/supercond_filter/supercond_filter_actualhousing.sonx");
        let hints = parse_sonnet_hints(xml).unwrap();

        assert!(!hints.conductor_polygons_m.is_empty());
        assert!(hints
            .conductor_polygons_m
            .iter()
            .any(|poly| poly.len() >= 4));
    }

        #[test]
        fn parse_supercond_control_hints() {
                let xml = include_str!("../../../testdata/sonnet/supercond_filter/supercond_filter_actualhousing.sonx");
                let hints = parse_sonnet_hints(xml).unwrap();

                assert_eq!(hints.sweep_type.as_deref(), Some("VARSWP"));
                assert_eq!(hints.matrix_solver.as_deref(), Some("AUTO"));
                assert_eq!(hints.speed_control.as_deref(), Some("MAX_ACCURACY"));
                assert_eq!(hints.precision_mode.as_deref(), Some("DOUBLE"));
                assert_eq!(hints.deembed_on, Some(true));
                assert_eq!(hints.output_folder.as_deref(), Some("."));
                assert_eq!(hints.subs_per_lambda, Some(100.0));
                assert!(hints.y_direction_negative);
        }

            #[test]
            fn parse_supercond_refplanes() {
                let xml = include_str!("../../../testdata/sonnet/supercond_resonators/1_QtrWave_SuperCondResonator.sonx");
                let hints = parse_sonnet_hints(xml).unwrap();

                assert_eq!(hints.ref_planes.len(), 2);
                assert_eq!(hints.ref_planes[0].side, "LEFT");
                assert_eq!(hints.ref_planes[0].plane_type, "NONE");
                let left_cal = hints.ref_planes[0].cal_length_m.unwrap();
                assert!((left_cal - 200.0e-6).abs() < 1.0e-12);
            }

        #[test]
        fn singular_tol_prefers_top_level_control_precision() {
                let xml = r#"
<SonnetProject>
    <Control>
        <Precision>DOUBLE</Precision>
        <UFFTControl>
            <Precision>MIN</Precision>
            <MatrixSolver>AUTO</MatrixSolver>
        </UFFTControl>
    </Control>
    <Geometry YDirection="Negative">
        <Box>
            <Size X="10" Y="5"/>
            <NumCells X="10" Y="5"/>
            <LocalOrigin Y="5"/>
        </Box>
    </Geometry>
    <Sweeps>
        <Set>
            <Frequencies>
                <Sweep Start="1" Stop="2" Step="0.1"/>
            </Frequencies>
        </Set>
    </Sweeps>
</SonnetProject>
"#;
                let hints = parse_sonnet_hints(xml).unwrap();
                let tol = map_singular_tol(hints.precision_mode.as_deref(), hints.speed_control.as_deref());
                assert!((tol - 1.0e-6).abs() < f64::EPSILON);
                assert!(hints.y_direction_negative);
        }

    #[test]
    fn cap_cells_preserves_ratio_and_limit() {
        let (nx, ny) = cap_cells(1600, 1200);
        assert!(nx * ny <= MAX_GENERATED_CELLS_PRODUCT);
        assert!(nx >= 2 && ny >= 2);
    }

    #[test]
    fn center_parser_works() {
        let c = parse_center_xy("(196.5,202)", 1.0e-6).unwrap();
        assert!((c.0 - 196.5e-6).abs() < 1.0e-12);
        assert!((c.1 - 202.0e-6).abs() < 1.0e-12);
    }

    #[test]
    fn points_parser_works() {
        let pts = parse_points_list("(0,1) (2,3) (4,5)", 1.0);
        assert_eq!(pts.len(), 3);
        assert_eq!(pts[1], (2.0, 3.0));
    }

    #[test]
    fn infer_direction_from_polygon_vertex_works() {
        let poly = vec![(0.0, 0.0), (10.0, 0.0), (10.0, 5.0), (0.0, 5.0)];
        let d = infer_direction_from_polygon_vertex(&poly, Some(2)).unwrap();
        assert!(d == "x" || d == "y");
    }
}
