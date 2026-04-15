use crate::examples;
use wasm_bindgen::JsCast;
use web_sys::{CanvasRenderingContext2d, HtmlCanvasElement, MouseEvent, WheelEvent};
use yew::prelude::*;

use rmsh_geo::{classify::classify, extract::{extract_surface, extract_surface_colored, extract_wireframe}};
use rmsh_renderer::{background_rgb, build_preview_frame, CameraExt, OrbitCamera, RenderConfig};

#[derive(Clone, PartialEq, Properties)]
pub struct MeshSceneProps {
    pub example_key: String,
}

#[derive(Clone)]
struct PreviewGeometry {
    surface_positions: Vec<[f32; 3]>,
    surface_colors: Vec<[f32; 3]>,
    surface_indices: Vec<u32>,
    wire_positions: Vec<[f32; 3]>,
    wire_indices: Vec<u32>,
    bounds_center: [f32; 3],
    bounds_diag: f32,
    node_count: usize,
    elem_count: usize,
    vol_count: usize,
    surf_count: usize,
    edge_count: usize,
}

#[derive(Clone, Copy)]
enum DragMode {
    Rotate,
    Pan,
}

fn compute_bounds(points: &[[f32; 3]]) -> ([f32; 3], f32) {
    if points.is_empty() {
        return ([0.0, 0.0, 0.0], 1.0);
    }

    let mut min_v = points[0];
    let mut max_v = points[0];
    for p in points.iter().skip(1) {
        min_v[0] = min_v[0].min(p[0]);
        min_v[1] = min_v[1].min(p[1]);
        min_v[2] = min_v[2].min(p[2]);
        max_v[0] = max_v[0].max(p[0]);
        max_v[1] = max_v[1].max(p[1]);
        max_v[2] = max_v[2].max(p[2]);
    }

    let center = [
        (min_v[0] + max_v[0]) * 0.5,
        (min_v[1] + max_v[1]) * 0.5,
        (min_v[2] + max_v[2]) * 0.5,
    ];
    let dx = max_v[0] - min_v[0];
    let dy = max_v[1] - min_v[1];
    let dz = max_v[2] - min_v[2];
    let diag = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-4);

    (center, diag)
}

fn draw_mesh(
    canvas: &HtmlCanvasElement,
    ctx: &CanvasRenderingContext2d,
    geom: Option<&PreviewGeometry>,
    camera: &OrbitCamera,
    cfg: &RenderConfig,
    status: &str,
) {
    let width = canvas.client_width().max(1) as u32;
    let height = canvas.client_height().max(1) as u32;
    if canvas.width() != width {
        canvas.set_width(width);
    }
    if canvas.height() != height {
        canvas.set_height(height);
    }

    let w = width as f64;
    let h = height as f64;

    let grad = ctx.create_linear_gradient(0.0, 0.0, 0.0, h);
    let (top, bot) = background_rgb(cfg);
    let top_color = format!("rgb({},{},{})", top[0], top[1], top[2]);
    let bot_color = format!("rgb({},{},{})", bot[0], bot[1], bot[2]);
    let _ = grad.add_color_stop(0.0, &top_color);
    let _ = grad.add_color_stop(1.0, &bot_color);
    ctx.set_fill_style_canvas_gradient(&grad);
    ctx.fill_rect(0.0, 0.0, w, h);

    if let Some(g) = geom {
        let frame = build_preview_frame(
            camera,
            cfg,
            &g.surface_positions,
            &g.surface_colors,
            &g.surface_indices,
            &g.wire_positions,
            &g.wire_indices,
            g.node_count,
            g.elem_count,
            g.vol_count,
            g.surf_count,
            g.edge_count,
            Some(status),
            width,
            height,
        );

        for overlay in &frame.overlay_lines {
            let style = format!(
                "rgba({},{},{},{})",
                overlay.rgba[0],
                overlay.rgba[1],
                overlay.rgba[2],
                (overlay.rgba[3] as f32 / 255.0).clamp(0.0, 1.0)
            );
            ctx.set_stroke_style_str(&style);
            ctx.set_line_width(overlay.width as f64);
            ctx.begin_path();
            ctx.move_to(overlay.start[0] as f64, overlay.start[1] as f64);
            ctx.line_to(overlay.end[0] as f64, overlay.end[1] as f64);
            ctx.stroke();
        }

        if cfg.show_faces {
            for tri in &frame.triangles {
                let face_style = format!(
                    "rgba({},{},{},{})",
                    tri.color_rgb[0],
                    tri.color_rgb[1],
                    tri.color_rgb[2],
                    cfg.surface_opacity.clamp(0.05, 1.0)
                );
                ctx.set_fill_style_str(&face_style);
                ctx.begin_path();
                ctx.move_to(tri.points[0][0] as f64, tri.points[0][1] as f64);
                ctx.line_to(tri.points[1][0] as f64, tri.points[1][1] as f64);
                ctx.line_to(tri.points[2][0] as f64, tri.points[2][1] as f64);
                ctx.close_path();
                ctx.fill();
            }
        }

        if cfg.show_edges {
            let edge_style = format!(
                "rgb({},{},{})",
                frame.edge_rgb[0], frame.edge_rgb[1], frame.edge_rgb[2]
            );
            ctx.set_stroke_style_str(&edge_style);
            ctx.set_line_width(1.0);

            for line in &frame.lines {
                ctx.begin_path();
                ctx.move_to(line.start[0] as f64, line.start[1] as f64);
                ctx.line_to(line.end[0] as f64, line.end[1] as f64);
                ctx.stroke();
            }
        }

        for text in &frame.overlay_texts {
            let style = format!(
                "rgba({},{},{},{})",
                text.rgba[0],
                text.rgba[1],
                text.rgba[2],
                (text.rgba[3] as f32 / 255.0).clamp(0.0, 1.0)
            );
            ctx.set_fill_style_str(&style);
            ctx.set_font(&format!("{}px ui-monospace, SFMono-Regular, Menlo, Consolas, monospace", text.font_px));
            let _ = ctx.fill_text(&text.text, text.position[0] as f64, text.position[1] as f64);
        }
    } else {
        let frame = build_preview_frame(
            camera,
            cfg,
            &[],
            &[],
            &[],
            &[],
            &[],
            0,
            0,
            0,
            0,
            0,
            Some(status),
            width,
            height,
        );

        for text in &frame.overlay_texts {
            let style = format!(
                "rgba({},{},{},{})",
                text.rgba[0],
                text.rgba[1],
                text.rgba[2],
                (text.rgba[3] as f32 / 255.0).clamp(0.0, 1.0)
            );
            ctx.set_fill_style_str(&style);
            ctx.set_font(&format!("{}px ui-monospace, SFMono-Regular, Menlo, Consolas, monospace", text.font_px));
            let _ = ctx.fill_text(&text.text, text.position[0] as f64, text.position[1] as f64);
        }
    }
}

#[function_component(MeshScene)]
pub fn mesh_scene(props: &MeshSceneProps) -> Html {
    let canvas_ref = use_node_ref();
    let geom = use_state(|| None::<PreviewGeometry>);
    let status = use_state(|| "Preparing mesh preview...".to_string());
    let redraw_tick = use_state(|| 0u64);
    let color_by_topology = use_state(|| true);
    let show_faces = use_state(|| true);
    let show_edges = use_state(|| true);
    let show_axes = use_state(|| true);

    let camera_ref = use_mut_ref(|| {
        let mut cam = OrbitCamera::new();
        cam.set_isometric();
        cam
    });
    let drag_state = use_mut_ref(|| false);
    let last_pos = use_mut_ref(|| (0f32, 0f32));
    let drag_mode = use_mut_ref(|| DragMode::Rotate);

    {
        let key = props.example_key.clone();
        let geom = geom.clone();
        let status = status.clone();
        let camera_ref = camera_ref.clone();
        let redraw_tick = redraw_tick.clone();
        let color_by_topology = color_by_topology.clone();
        use_effect_with((key, *color_by_topology), move |(k, topo_colored)| {
            status.set(format!("Loading mesh for {}...", k));
            let bytes = examples::get_mesh_bytes(k);
            match rmsh_io::load_msh_from_bytes(&bytes) {
                Ok(mesh) => {
                    let surface = if *topo_colored {
                        let topo = classify(&mesh, 40.0);
                        extract_surface_colored(&mesh, &topo)
                    } else {
                        extract_surface(&mesh)
                    };
                    let wire = extract_wireframe(&mesh, &[1, 2, 3]);
                    let positions = if !surface.positions.is_empty() {
                        surface.positions.clone()
                    } else {
                        wire.positions.clone()
                    };
                    let (center, diag) = compute_bounds(&positions);

                    {
                        let mut cam = camera_ref.borrow_mut();
                        cam.fit_to_bbox(center, diag);
                        cam.set_isometric();
                    }

                    let built = PreviewGeometry {
                        surface_positions: surface.positions,
                        surface_colors: surface.colors,
                        surface_indices: surface.indices,
                        wire_positions: wire.positions,
                        wire_indices: wire.indices,
                        bounds_center: center,
                        bounds_diag: diag,
                        node_count: mesh.node_count(),
                        elem_count: mesh.element_count(),
                        vol_count: mesh.elements_by_dimension(3).len(),
                        surf_count: mesh.elements_by_dimension(2).len(),
                        edge_count: mesh.elements_by_dimension(1).len(),
                    };
                    geom.set(Some(built));
                    status.set("Mesh preview ready".to_string());
                    redraw_tick.set(js_sys::Date::now() as u64);
                }
                Err(err) => {
                    geom.set(None);
                    status.set(format!("Mesh preview failed: {}", err));
                }
            }
            || ()
        });
    }

    {
        let canvas_ref = canvas_ref.clone();
        let geom = geom.clone();
        let status = status.clone();
        let camera_ref = camera_ref.clone();
        let show_faces = show_faces.clone();
        let show_edges = show_edges.clone();
        let show_axes = show_axes.clone();
        let _redraw_tick = *redraw_tick;
        use_effect(move || {
            if let Some(canvas) = canvas_ref.cast::<HtmlCanvasElement>() {
                if let Ok(Some(raw)) = canvas.get_context("2d") {
                    if let Ok(ctx) = raw.dyn_into::<CanvasRenderingContext2d>() {
                        let mut cfg = RenderConfig::default();
                        cfg.show_faces = *show_faces;
                        cfg.show_edges = *show_edges;
                        cfg.show_axes = *show_axes;
                        let cam = camera_ref.borrow();
                        draw_mesh(&canvas, &ctx, geom.as_ref(), &cam, &cfg, &status);
                    }
                }
            }
            || ()
        });
    }

    let on_mouse_down = {
        let drag_state = drag_state.clone();
        let last_pos = last_pos.clone();
        let drag_mode = drag_mode.clone();
        Callback::from(move |e: MouseEvent| {
            *drag_state.borrow_mut() = true;
            *last_pos.borrow_mut() = (e.offset_x() as f32, e.offset_y() as f32);
            *drag_mode.borrow_mut() = if e.shift_key() || e.button() == 1 {
                DragMode::Pan
            } else {
                DragMode::Rotate
            };
        })
    };

    let on_mouse_up = {
        let drag_state = drag_state.clone();
        Callback::from(move |_e: MouseEvent| {
            *drag_state.borrow_mut() = false;
        })
    };

    let on_mouse_leave = {
        let drag_state = drag_state.clone();
        Callback::from(move |_e: MouseEvent| {
            *drag_state.borrow_mut() = false;
        })
    };

    let on_mouse_move = {
        let drag_state = drag_state.clone();
        let last_pos = last_pos.clone();
        let camera_ref = camera_ref.clone();
        let drag_mode = drag_mode.clone();
        let redraw_tick = redraw_tick.clone();
        Callback::from(move |e: MouseEvent| {
            if !*drag_state.borrow() {
                return;
            }
            let x = e.offset_x() as f32;
            let y = e.offset_y() as f32;
            let (lx, ly) = *last_pos.borrow();
            let dx = x - lx;
            let dy = y - ly;
            *last_pos.borrow_mut() = (x, y);

            {
                let mut cam = camera_ref.borrow_mut();
                match *drag_mode.borrow() {
                    DragMode::Rotate => cam.rotate(dx * 0.01, -dy * 0.01),
                    DragMode::Pan => cam.pan(dx, dy),
                }
            }
            redraw_tick.set(js_sys::Date::now() as u64);
        })
    };

    let on_wheel = {
        let camera_ref = camera_ref.clone();
        let redraw_tick = redraw_tick.clone();
        Callback::from(move |e: WheelEvent| {
            e.prevent_default();
            let delta = (e.delta_y() as f32 / 1000.0).clamp(-0.2, 0.2);
            {
                let mut cam = camera_ref.borrow_mut();
                cam.zoom(delta);
            }
            redraw_tick.set(js_sys::Date::now() as u64);
        })
    };

    let on_reset_view = {
        let camera_ref = camera_ref.clone();
        let geom = geom.clone();
        let redraw_tick = redraw_tick.clone();
        Callback::from(move |_e: MouseEvent| {
            let mut cam = camera_ref.borrow_mut();
            if let Some(g) = geom.as_ref() {
                cam.fit_to_bbox(g.bounds_center, g.bounds_diag);
            }
            cam.set_isometric();
            redraw_tick.set(js_sys::Date::now() as u64);
        })
    };

    let on_toggle_topology_color = {
        let color_by_topology = color_by_topology.clone();
        Callback::from(move |e: Event| {
            let input: web_sys::HtmlInputElement = e.target_unchecked_into();
            color_by_topology.set(input.checked());
        })
    };

    let on_toggle_faces = {
        let show_faces = show_faces.clone();
        let redraw_tick = redraw_tick.clone();
        Callback::from(move |e: Event| {
            let input: web_sys::HtmlInputElement = e.target_unchecked_into();
            show_faces.set(input.checked());
            redraw_tick.set(js_sys::Date::now() as u64);
        })
    };

    let on_toggle_edges = {
        let show_edges = show_edges.clone();
        let redraw_tick = redraw_tick.clone();
        Callback::from(move |e: Event| {
            let input: web_sys::HtmlInputElement = e.target_unchecked_into();
            show_edges.set(input.checked());
            redraw_tick.set(js_sys::Date::now() as u64);
        })
    };

    let on_toggle_axes = {
        let show_axes = show_axes.clone();
        let redraw_tick = redraw_tick.clone();
        Callback::from(move |e: Event| {
            let input: web_sys::HtmlInputElement = e.target_unchecked_into();
            show_axes.set(input.checked());
            redraw_tick.set(js_sys::Date::now() as u64);
        })
    };

    html! {
        <div class="mesh-scene-wrap">
            <div class="mesh-scene-toolbar">
                <label class="mesh-scene-toggle">
                    <input type="checkbox" checked={*color_by_topology} onchange={on_toggle_topology_color} />
                    {"Topo Color"}
                </label>
                <label class="mesh-scene-toggle">
                    <input type="checkbox" checked={*show_faces} onchange={on_toggle_faces} />
                    {"Faces"}
                </label>
                <label class="mesh-scene-toggle">
                    <input type="checkbox" checked={*show_edges} onchange={on_toggle_edges} />
                    {"Edges"}
                </label>
                <label class="mesh-scene-toggle">
                    <input type="checkbox" checked={*show_axes} onchange={on_toggle_axes} />
                    {"Axes"}
                </label>
                <button type="button" class="collapse-btn" onclick={on_reset_view}>{"Reset View"}</button>
            </div>
            <canvas
                class="mesh-canvas"
                ref={canvas_ref}
                onmousedown={on_mouse_down}
                onmouseup={on_mouse_up}
                onmouseleave={on_mouse_leave}
                onmousemove={on_mouse_move}
                onwheel={on_wheel}
            />
        </div>
    }
}
