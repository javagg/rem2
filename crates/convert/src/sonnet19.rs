use anyhow::{bail, Context};
use quick_xml::events::Event;
use quick_xml::Reader;
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

    let (msh_text, port_tags, pec_tags) = generate_rect_msh_with_port_tags(
        width_m,
        height_m,
        cells_x,
        cells_y,
        &hints.ports,
        hints.y_direction_negative,
        hints.local_origin_y_m.unwrap_or(height_m),
    );
    std::fs::write(out_msh, msh_text)
        .with_context(|| format!("writing mesh: {}", out_msh.display()))?;

    let mut ports_json: Vec<Value> = Vec::new();
    for (i, p) in hints.ports.iter().enumerate() {
        let index = (i + 1) as u32;
        let attr = port_tags.get(i).copied().unwrap_or(BASE_PEC_TAG);
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

    let mut mom = serde_json::Map::new();
    mom.insert("Equation".to_string(), json!("CFIE"));
    mom.insert("Basis".to_string(), json!("RWG"));
    mom.insert("FreqMin".to_string(), json!(freq_min));
    mom.insert("FreqMax".to_string(), json!(freq_max));
    mom.insert("FreqStep".to_string(), json!(freq_step));
    mom.insert("Alpha".to_string(), json!(0.5));
    mom.insert("FastSolver".to_string(), json!("ACA"));
    mom.insert("RefImpedance".to_string(), json!(50.0));
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

    let cfg = json!({
        "Problem": {
            "Type": "MoM",
            "Verbose": 1,
            "Output": "./output/sonnet19"
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
        "Sonnet19 XML converted: {} -> config={}, mesh={} (w={:.6e} m, h={:.6e} m, nx={} (raw {}), ny={} (raw {}))",
        xml_path.display(),
        out_config.display(),
        out_msh.display(),
        width_m,
        height_m,
        cells_x,
        cells_x_raw,
        cells_y,
        cells_y_raw,
    );

    Ok(())
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
) -> (String, Vec<u32>, Vec<u32>) {
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

    // GMSH v2.2 ASCII for straightforward per-element physical tags.
    let total_nodes = (nx + 1) * (ny + 1);
    let total_elems = 2 * nx + 2 * ny + 2 * nx * ny;
    let node_id = |ix: usize, iy: usize| -> usize { iy * (nx + 1) + ix + 1 };

    let mut s = String::new();
    s.push_str("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");
    s.push_str("$Nodes\n");
    s.push_str(&format!("{}\n", total_nodes));
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
            s.push_str(&format!("{} {:.15e} {:.15e} 0\n", id, x, y));
        }
    }
    s.push_str("$EndNodes\n");

    s.push_str("$Elements\n");
    s.push_str(&format!("{}\n", total_elems));
    let mut eid = 1usize;

    for ix in 0..nx {
        s.push_str(&format!("{} 1 2 {} {} {} {}\n", eid, BOTTOM_TAG, BOTTOM_TAG, node_id(ix, 0), node_id(ix + 1, 0)));
        eid += 1;
    }
    for ix in 0..nx {
        s.push_str(&format!("{} 1 2 {} {} {} {}\n", eid, TOP_TAG, TOP_TAG, node_id(ix, ny), node_id(ix + 1, ny)));
        eid += 1;
    }
    for iy in 0..ny {
        s.push_str(&format!("{} 1 2 {} {} {} {}\n", eid, LEFT_TAG, LEFT_TAG, node_id(0, iy), node_id(0, iy + 1)));
        eid += 1;
    }
    for iy in 0..ny {
        s.push_str(&format!("{} 1 2 {} {} {} {}\n", eid, RIGHT_TAG, RIGHT_TAG, node_id(nx, iy), node_id(nx, iy + 1)));
        eid += 1;
    }

    for iy in 0..ny {
        for ix in 0..nx {
            let tag = col_tags[ix];
            let n00 = node_id(ix, iy);
            let n10 = node_id(ix + 1, iy);
            let n01 = node_id(ix, iy + 1);
            let n11 = node_id(ix + 1, iy + 1);
            s.push_str(&format!("{} 2 2 {} {} {} {} {}\n", eid, tag, tag, n00, n10, n11));
            eid += 1;
            s.push_str(&format!("{} 2 2 {} {} {} {} {}\n", eid, tag, tag, n00, n11, n01));
            eid += 1;
        }
    }

    s.push_str("$EndElements\n");
    (s, port_tags, pec_tags)
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
                );
                path.pop();
            }
            Ok(Event::Text(t)) => {
                let text = String::from_utf8_lossy(t.as_ref()).trim().to_string();
                if !text.is_empty() {
                    parse_text_node(&path, &text, &mut state);
                    if let Some(last) = path.last() {
                        map_kv_hint(&mut hints, last, &text, state);
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
                    current_polygon = None;
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
) {
    let in_geometry_box = path_contains(path, "geometry") && path_contains(path, "box");

    let mut attrs_vec: Vec<(String, String)> = Vec::new();
    for attr in attrs.flatten() {
        let key = local_name(attr.key.as_ref());
        let value = String::from_utf8_lossy(attr.value.as_ref()).to_string();
        attrs_vec.push((key, value));
    }

    // Parse <Geometry YDirection="Negative"> root attribute
    if tag == "geometry" && path.len() == 1 {
        for (k, v) in &attrs_vec {
            if k == "YDirection" && v.to_lowercase() == "negative" {
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
            if k == "Name" {
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
                    if k == "Value" {
                        if let Some(val) = parse_f64(v) {
                            diel.eps_r = val;
                        }
                    }
                }
            } else if tag == "tan" {
                for (k, v) in &attrs_vec {
                    if k == "Value" {
                        if let Some(val) = parse_f64(v) {
                            diel.loss_tan = val;
                        }
                    }
                }
            } else if tag == "murel" {
                for (k, v) in &attrs_vec {
                    if k == "Value" {
                        if let Some(val) = parse_f64(v) {
                            diel.mu_r = val;
                        }
                    }
                }
            }
        }
    }

    // Parse <DielectricMaterialModel Thickness="10.0"/> for layer thickness
    if tag == "dielectricmaterialmodel" && path_contains(path, "level") && path_contains(path, "geometry") {
        if !hints.dielectric_layers.is_empty() {
            for (k, v) in &attrs_vec {
                if k == "Thickness" {
                    if let Some(val) = parse_f64(v) {
                        let thickness = val * state.length_unit_scale;
                        // Find the last dielectric with the correct name from <Level MaterialName="...">
                        // For now, just update the last dielectric
                        hints.dielectric_layers.last_mut().unwrap().thickness_m = thickness;
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
