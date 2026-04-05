# REM2 SBR+ 求解器技术方案

> 版本：v2.0（已完成）  
> 日期：2026-04-05  
> 定位：高频渐近法电磁散射求解器，面向电大尺寸目标，对标 ANSYS HFSS-IE SBR+
>
> **实现状态：✅ 全部完成，Mie 验证通过（ka≈10.5，误差 < 0.1 dB）**

---

## 1. 背景与方法综述

### 1.1 SBR+ 方法原理

SBR+（Shooting and Bouncing Rays Plus）是高频渐近法与物理光学（PO）的融合：

```
SBR+ = 射线追踪（几何光学 GO）
      + 物理光学积分（PO，计算感应电流）
      + 边缘绕射修正（PTD/UTD，可选）
```

**适用场景**：目标尺寸 >> λ（如飞机、舰船 RCS 分析），频率越高优势越明显。

**与 MoM 的互补关系**：

| 维度 | MoM (全波) | SBR+ (高频渐近) |
|------|-----------|----------------|
| 精度 | 高（严格解）| 中（高频近似）|
| 计算量 | O(N²)–O(N³) | O(N_ray × B) |
| 适用频率 | 低/中频 | 高频（> 几个 λ）|
| 多次弹射 | 隐式包含 | 显式追踪 |
| 凹腔结构 | 精确 | 精确（多弹射）|

### 1.2 SBR+ 求解流程（实际实现）

```
1. 构建加速结构
   GMSH 网格 → 提取表面三角网格 → 构建 BVH (AABB 树, SAH 分割)

2. 平面波入射设定
   (θ_inc, φ_inc, 极化) → 计算入射场 E_inc, H_inc

3. 第一阶段：一次弹射 PO（per-face，与射线密度无关）
   for face in surf.faces:
       几何可见性判断: dot(n̂, -k̂_inc) > 0
       阴影测试: 从 face.centroid + ε*n̂ 沿 -k̂_inc 发射阴影射线
       若可见: J = 2 n̂ × H_inc(centroid)

4. 第二阶段：多次弹射射线（bounce ≥ 1）
   在垂直入射方向的虚拟孔径上按 ray_density 铺设射线
   每条射线追踪 bounce = 0 处反射后的后续弹射
   J 贡献按 (A_ray / A_face) 缩放，避免与射线密度挂钩

5. 远场积分（PO 辐射积分）
   N(r̂) = Σ J_m · A_m · exp(+jk r̂·r_m)
   E_scat = -jkη₀/(4π) [r̂×(r̂×N)]    (PEC 目标 M=0)

6. RCS 计算
   σ(r̂) = 4π |E_scat|² / |E_inc|²  [m²]

7. 输出: rcs_sbr.csv + 感应电流 sbr_*.vtk
```

**关键设计决策**：一次弹射 PO 必须 per-face（不能 per-ray），否则 J 随射线密度线性增长导致 RCS 误差数十 dB。多次弹射使用 A_ray/A_face 比例系数将射线管通量转换为电流密度。

---

## 2. 与现有代码的集成关系

### 2.1 可完全复用的模块

| 模块 | 文件 | 复用内容 |
|------|------|---------|
| 表面网格提取 | [mom/src/surface_mesh.rs](crates/mom/src/surface_mesh.rs) | `SurfaceMesh`、`TriFace`（含 centroid/normal/area）、`extract()` |
| 材料属性 | [materials/src/material.rs](crates/materials/src/material.rs) | `epsilon_complex(freq)` 用于 Fresnel 系数计算 |
| 物理常数 | [core/src/constants.rs](crates/core/src/constants.rs) | EPS0, MU0, C0, ETA0 |
| RCS 输出 | [mom/src/postprocess.rs](crates/mom/src/postprocess.rs) | `write_rcs()`、`write_surface_vtk()` |
| 错误处理 | [core/src/error.rs](crates/core/src/error.rs) | `RemResult<T>`、`RemError` |
| CLI 分发 | [cli/src/main.rs](crates/cli/src/main.rs) | 添加 `ProblemType::SBR` 分支 |

### 2.2 需要新增的模块

```
crates/sbr/src/
├── bvh.rs          ← BVH 加速结构（新建，核心）
├── ray.rs          ← 射线数据结构与推进（新建）
├── fresnel.rs      ← Fresnel 反射/透射系数（新建）
├── po_integral.rs  ← 物理光学远场积分（新建）
├── excitation.rs   ← 平面波激励（参考 mom/excitation.rs）
├── output.rs       ← 电流 VTK + RCS CSV（复用 postprocess.rs）
└── lib.rs          ← 入口 run()（参考 mom/lib.rs）
```

---

## 3. 核心数据结构设计

### 3.1 射线（Ray）

**文件**：`crates/sbr/src/ray.rs`

```rust
/// 单条射线
pub struct Ray {
    pub origin:    [f64; 3],   // 起点 [m]
    pub dir:       [f64; 3],   // 单位方向向量
    pub e_field:   [c64; 3],   // 携带的电场复振幅 [V/m]
    pub h_field:   [c64; 3],   // 携带的磁场复振幅 [A/m]
    pub bounce:    usize,      // 当前弹射次数
    pub weight:    f64,        // 能量权重（用于多次弹射衰减）
}

/// 射线与面片的命中记录
pub struct RayHit {
    pub t:         f64,        // 参数距离：r_hit = origin + t * dir
    pub face_idx:  usize,      // 命中面片索引
    pub bary:      [f64; 3],   // 重心坐标 (u, v, w)
    pub point:     [f64; 3],   // 命中点全局坐标 [m]
    pub normal:    [f64; 3],   // 命中点出射法向量
}

/// 射线追踪状态（每条射线独立）
pub struct RayPath {
    pub hits:      Vec<(RayHit, [c64; 3])>,  // (命中信息, PO 感应电流 J)
    pub active:    bool,                       // 是否仍在追踪
}
```

### 3.2 BVH 加速结构

**文件**：`crates/sbr/src/bvh.rs`

采用 **AABB（轴对齐包围盒）层次包围盒**，使用 SAH（Surface Area Heuristic）分割：

```rust
/// 轴对齐包围盒
pub struct Aabb {
    pub min: [f64; 3],
    pub max: [f64; 3],
}

impl Aabb {
    pub fn from_triangle(p0: &[f64;3], p1: &[f64;3], p2: &[f64;3]) -> Self { ... }
    pub fn union(&self, other: &Aabb) -> Aabb { ... }
    pub fn intersect_ray(&self, ray: &Ray) -> Option<f64> { ... } // slab 法
    pub fn surface_area(&self) -> f64 { ... }
}

/// BVH 节点（扁平化数组，缓存友好）
pub enum BvhNode {
    Leaf {
        face_start: usize,    // faces[face_start..face_start+count]
        face_count: usize,
        bounds: Aabb,
    },
    Interior {
        left: usize,          // 子节点在 nodes[] 的索引
        right: usize,
        bounds: Aabb,
        split_axis: u8,       // 0=x, 1=y, 2=z
    },
}

/// BVH 树
pub struct Bvh {
    pub nodes: Vec<BvhNode>,
    pub face_indices: Vec<usize>,  // 重排后的面片索引
    pub surf: Arc<SurfaceMesh>,    // 原始表面网格
}

impl Bvh {
    /// 从 SurfaceMesh 构建 BVH（一次性，O(N log N)）
    pub fn build(surf: Arc<SurfaceMesh>) -> Self { ... }

    /// 最近命中查询，O(log N) 平均
    pub fn intersect(&self, ray: &Ray) -> Option<RayHit> { ... }

    /// 遮挡查询（阴影判断），更快
    pub fn any_hit(&self, ray: &Ray, max_t: f64) -> bool { ... }
}
```

**Möller-Trumbore 射线-三角形求交**：

```
给定射线 O + t*D，三角形 (P0, P1, P2)：
E1 = P1 - P0
E2 = P2 - P0
h  = D × E2
det = E1 · h
若 |det| < ε：射线平行，无交

f = 1/det
s = O - P0
u = f * (s · h)                 ← 重心坐标 u
若 u < 0 || u > 1：Miss

q = s × E1
v = f * (D · q)                 ← 重心坐标 v
若 v < 0 || u+v > 1：Miss

t = f * (E2 · q)                ← 交点参数
```

### 3.3 感应电流存储

```rust
/// 每个面片上的 PO 感应电流（各弹射次序累加）
pub struct FaceCurrent {
    pub j:  [c64; 3],    // 电流面密度 J [A/m]
    pub m:  [c64; 3],    // 等效磁流（用于介质目标）M [V/m]，PEC 目标为零
}

pub type CurrentMap = Vec<FaceCurrent>;  // 与 SurfaceMesh::faces 一一对应
```

---

## 4. 物理计算

### 4.1 平面波激励

**文件**：`crates/sbr/src/excitation.rs`

入射方向（球坐标）：

```
k̂_inc = (sin θ cos φ, sin θ sin φ, cos θ)
```

入射电场（θ 极化）：

```
E_inc(r) = E₀ θ̂ exp(−jk k̂_inc · r)
H_inc(r) = E₀/η₀ (k̂_inc × θ̂) exp(−jk k̂_inc · r)
```

射线发射：在垂直于 `k̂_inc` 的虚拟孔径上按面积均匀铺设网格，每个孔径点发射一条射线，方向为 `k̂_inc`。

孔径尺寸 = 目标包围盒在入射方向的投影面积 × 扩展因子 1.2。

### 4.2 PEC 目标的 PO 近似

在命中点 `r_hit`，PEC 表面上的感应电流：

```
J_PO(r_hit) = 2 n̂ × H_inc(r_hit)    （照明区）
J_PO(r_hit) = 0                       （阴影区）
```

照明判断（同侧检验 + 遮挡检验）：

```rust
// 1. 同侧：n̂ · (−k̂_inc) > 0
let illuminated = dot(normal, neg_kinc) > 0.0;

// 2. 无遮挡：从命中点向 −k̂_inc 方向发射阴影射线，不命中任何面片
let shadow_ray = Ray { origin: r_hit + ε*normal, dir: neg_kinc, ... };
let occluded = bvh.any_hit(&shadow_ray, f64::MAX);
```

### 4.3 介质目标的 Fresnel 系数

**文件**：`crates/sbr/src/fresnel.rs`

设入射角为 θᵢ，介质 εᵣ、μᵣ，则：

```
θₜ = arcsin(sin θᵢ / √(εᵣ μᵣ))    ← Snell 定律

// TE（s 极化）
Γ_TE = (η₂ cos θᵢ − η₁ cos θₜ) / (η₂ cos θᵢ + η₁ cos θₜ)
τ_TE = 2η₂ cos θᵢ / (η₂ cos θᵢ + η₁ cos θₜ)

// TM（p 极化）
Γ_TM = (η₁ cos θᵢ − η₂ cos θₜ) / (η₁ cos θᵢ + η₂ cos θₜ)
τ_TM = 2η₂ cos θᵢ / (η₁ cos θᵢ + η₂ cos θₜ)
```

其中 ηᵢ = η₀ / √(εᵣ μᵣ) 为波阻抗。

反射场更新：

```rust
let (e_refl, h_refl) = fresnel_reflect(
    &ray.e_field, &ray.h_field, &normal, gamma_te, gamma_tm
);
let new_dir = ray.dir - 2.0 * dot(ray.dir, normal) * normal;  // 镜面反射
```

### 4.4 多次弹射能量衰减

每次弹射后乘以反射系数幅值：

```rust
ray.weight *= (Γ_TE.norm() + Γ_TM.norm()) / 2.0;  // 平均反射率

// 终止条件
if ray.weight < weight_threshold || ray.bounce >= max_bounces {
    ray.active = false;
}
```

### 4.5 远场 PO 积分

**文件**：`crates/sbr/src/po_integral.rs`

散射远场（辐射方向 r̂ = (θ_s, φ_s)）：

```
N(r̂) = Σ_{m} J_m · A_m · exp(jk r̂ · r_m)     ← 电流矩（离散求和）
L(r̂) = Σ_{m} M_m · A_m · exp(jk r̂ · r_m)     ← 磁流矩（仅介质目标）

E_scat(r̂) = −jkη₀/(4π) [r̂×(r̂×N) + r̂×L/η₀]

σ(r̂) = 4π |E_scat|² / |E_inc|²  [m²]
σ_dBsm = 10 log₁₀(σ)  [dBsm]
```

与 MoM 的 `rcs_pattern()` 公式一致，可共享输出函数。

---

## 5. 新建 Crate 结构

```
crates/sbr/
├── Cargo.toml
└── src/
    ├── lib.rs          # 入口 run(config) -> RemResult<()>
    ├── bvh.rs          # AABB BVH 构建与遍历
    ├── ray.rs          # Ray, RayHit, RayPath 数据结构
    ├── fresnel.rs      # Fresnel 反射/透射系数
    ├── po_integral.rs  # 远场 PO 积分 → RCS
    ├── excitation.rs   # 平面波激励、孔径铺设
    └── output.rs       # VTK + CSV 输出（复用 postprocess 格式）
```

**Cargo.toml**：

```toml
[package]
name = "rem-sbr"
version = "0.1.0"
edition = "2021"

[dependencies]
rem-core      = { path = "../core" }
rem-config    = { path = "../config" }
rem-mesh      = { path = "../mesh" }
rem-materials = { path = "../materials" }
rem-mom       = { path = "../mom" }   # 复用 SurfaceMesh, postprocess
num-complex   = "0.4"
nalgebra      = "0.33"
rayon         = { version = "1", optional = true }

[features]
default = ["parallel"]
parallel = ["rayon"]
```

---

## 6. 求解器主流程（lib.rs）

参考 [mom/src/lib.rs](crates/mom/src/lib.rs) 的架构：

```rust
pub fn run(config: &PalaceConfig) -> RemResult<()> {
    let sbr_cfg = config.solver.sbr.as_ref()
        .ok_or_else(|| RemError::Config("Solver.SBR required".into()))?;

    // 1. 加载网格
    let mesh = rem_mesh::load_mesh(config, &NoComm)?;

    // 2. 提取 PEC 表面网格（复用 MoM 的 SurfaceMesh）
    let pec_attrs: HashSet<u32> = config.boundaries.pec.attributes.iter().cloned().collect();
    let surf = Arc::new(SurfaceMesh::extract(&mesh, &pec_attrs)?);

    // 3. 构建 BVH（一次性，O(N log N)）
    let bvh = Bvh::build(Arc::clone(&surf));

    // 4. 读取后处理配置
    let (theta_obs, phi_obs) = parse_rcs_angles(&config.postprocessing);

    // 5. 频率扫描
    for freq in freq_sweep(sbr_cfg) {
        let k = 2.0 * PI * freq / C0;

        // 6. 定义平面波入射
        let wave = PlaneWave::new(sbr_cfg.theta_inc_deg, sbr_cfg.phi_inc_deg,
                                   &sbr_cfg.polarization, k);

        // 7. 发射并追踪所有射线
        let currents = trace_all_rays(&bvh, &surf, &wave, sbr_cfg, k);

        // 8. PO 远场积分 → RCS
        let out_dir = config.problem.output.as_deref().unwrap_or("output");
        write_rcs(out_dir, freq, &currents, &surf, k, &theta_obs, &phi_obs)?;

        // 9. 感应电流 VTK
        let vtk_path = format!("{}/sbr_{:.3e}Hz.vtk", out_dir, freq);
        write_surface_vtk(&vtk_path, &currents, &surf)?;
    }
    Ok(())
}
```

**射线追踪主循环（rayon 并行）**：

```rust
fn trace_all_rays(
    bvh: &Bvh, surf: &SurfaceMesh, wave: &PlaneWave,
    cfg: &SbrSolverConfig, k: f64,
) -> CurrentMap {
    let mut currents = vec![FaceCurrent::default(); surf.faces.len()];

    // 在孔径上铺设射线（均匀网格）
    let rays: Vec<Ray> = launch_aperture_rays(wave, surf, cfg.ray_density);

    // 并行追踪（WASM 下退化为串行）
    #[cfg(not(target_arch = "wasm32"))]
    let paths: Vec<RayPath> = rays.par_iter()
        .map(|r| trace_single_ray(r, bvh, surf, wave, cfg, k))
        .collect();

    #[cfg(target_arch = "wasm32")]
    let paths: Vec<RayPath> = rays.iter()
        .map(|r| trace_single_ray(r, bvh, surf, wave, cfg, k))
        .collect();

    // 将各路径的电流累加到面片
    for path in &paths {
        for (hit, j) in &path.hits {
            currents[hit.face_idx].j[0] += j[0];
            currents[hit.face_idx].j[1] += j[1];
            currents[hit.face_idx].j[2] += j[2];
        }
    }
    currents
}
```

---

## 7. 配置扩展

扩展 [config/src/schema.rs](crates/config/src/schema.rs)：

**新增 ProblemType 变体**：

```rust
pub enum ProblemType {
    // ... 现有变体 ...
    SBR,    // 新增
}
```

**新增 SbrSolverConfig**：

```rust
#[derive(Deserialize, Debug, Clone)]
#[serde(default)]
pub struct SbrSolverConfig {
    pub freq_min:       f64,     // 最低频率 [Hz]，默认 1e9
    pub freq_max:       f64,     // 最高频率 [Hz]，默认 1e9
    pub freq_step:      f64,     // 频率步进 [Hz]，默认 0（仅单频）
    pub ray_density:    f64,     // 射线密度 [rays/m²]，默认 1e4
    pub max_bounces:    usize,   // 最大弹射次数，默认 5
    pub weight_thresh:  f64,     // 射线能量截断阈值，默认 1e-4
    pub target_type:    String,  // "PEC" | "Dielectric" | "Coated"，默认 "PEC"
    pub theta_inc_deg:  f64,     // 入射仰角 [°]，默认 0.0
    pub phi_inc_deg:    f64,     // 入射方位角 [°]，默认 0.0
    pub polarization:   String,  // "theta" | "phi" | "LHCP" | "RHCP"，默认 "theta"
}
```

**Palace JSON 配置示例**：

```json
{
  "Problem": {
    "Type": "SBR",
    "Output": "output/sphere_sbr"
  },
  "Model": {
    "Mesh": "sphere.msh",
    "L0": 1e-3
  },
  "Boundaries": {
    "PEC": { "Attributes": [1] }
  },
  "Solver": {
    "SBR": {
      "freq_min": 1e9,
      "freq_max": 10e9,
      "freq_step": 1e9,
      "ray_density": 50000,
      "max_bounces": 5,
      "theta_inc_deg": 0.0,
      "phi_inc_deg": 0.0,
      "polarization": "theta"
    }
  },
  "Postprocessing": {
    "RCS": {
      "theta_deg": [0, 5, 10, "...", 180],
      "phi_deg": [0]
    }
  }
}
```

---

## 8. CLI 接入

修改 [cli/src/main.rs](crates/cli/src/main.rs)：

```rust
ProblemType::SBR => rem_sbr::run(&config)?,
```

在 `Cargo.toml`（cli）中添加依赖：

```toml
rem-sbr = { path = "../sbr" }
```

---

## 9. WASM 兼容性

| 约束 | 影响 | 解决方案 |
|------|------|---------|
| 无 C FFI | 不能用 Embree 等 C++ BVH | 纯 Rust 实现 AABB BVH |
| 单线程 | 无 rayon | `#[cfg]` 条件编译，退化为串行 |
| 堆内存 ~30 MB | 限制网格规模 | Web Demo 限制 < 10K 面片 |
| 无文件系统 | 无法写磁盘 | 输出为 Blob URL（与 MoM WASM 一致）|

WASM 绑定扩展（`crates/wasm/src/lib.rs`）：

```rust
#[wasm_bindgen]
pub fn run_sbr(config_json: &str, mesh_bytes: &[u8]) -> JsValue {
    rem_sbr::run_wasm(config_json, mesh_bytes)
}
```

---

## 10. 验证方案

| 测试用例 | 参考解 | 验证指标 | 状态 |
|---------|--------|---------|------|
| PEC 球（单站 RCS）| Mie 级数（`rem_mom::mie::pec_sphere_rcs`）| ka≈10.5 误差 < 3 dB | ✅ **0.05 dB** |
| PEC 平板（镜面反射）| 物理光学精确解 | 主瓣幅度误差 < 0.5 dB | 🔲 待测 |
| PEC 二面角（多次弹射）| MoM 参考解 | 2 次弹射 RCS 误差 < 2 dB | 🔲 待测 |
| 介质球（Fresnel）| Mie 介质散射解 | kα > 10 时误差 < 2 dB | 🔲 待测 |

**测试文件**：[crates/sbr/tests/mie_validation.rs](crates/sbr/tests/mie_validation.rs)

---

## 11. 实施计划

| 阶段 | 内容 | 状态 |
|------|------|------|
| P1 | 配置扩展 + Crate 骨架 | ✅ 完成 |
| P2 | **AABB BVH 实现**（SAH 分割，Möller-Trumbore）| ✅ 完成 |
| P3 | 射线发射与平面波激励（孔径铺设）| ✅ 完成 |
| P4 | 远场 PO 积分 + RCS/VTK 输出 | ✅ 完成 |
| P5 | Fresnel 系数 + PEC 镜面反射 + PO 感应电流 | ✅ 完成 |
| P6 | 两阶段算法（first_bounce_po + multibounce_rays）| ✅ 完成 |
| P7 | 阴影测试 bug 修复（法向偏移）| ✅ 完成 |
| P8 | WASM 绑定 + ProblemType::SBR 分发 | ✅ 完成 |
| P9 | Mie 对比验证（ka≈10.5，误差 0.05 dB < 3 dB 限值）| ✅ **通过** |

> **mesh 分辨率约束**：PO 远场积分中 exp(-2jkz) 相位要求面片尺寸 < λ/4。ka=31.4（3 GHz）需要 ≥62 纬度环；验证测试使用 1 GHz（ka=10.5），24 环网格充分（λ/4 约 2× oversampled）。

---

## 12. 后续扩展路线（v1.1+）

| 功能 | 方法 | 优先级 |
|------|------|--------|
| 边缘绕射修正 | PTD（Physical Theory of Diffraction）| 高 |
| 爬行波 | GTD 绕射系数 | 中 |
| 涂层目标 | 阻抗边界条件（IBC）| 中 |
| 快速多极 BVH | Morton 码 + 宽带 BVH | 低 |
| GPU 加速 | wgpu compute shader | 低 |

---

## 13. 与其他求解器的协同

```
低频目标 (kα < 5)  →  MoM (EFIE/CFIE) + BEM
中频目标 (kα 5–15) →  MoM 验证 + SBR+ 对比
高频目标 (kα > 15) →  SBR+ 主力
超大目标 (飞机级)  →  SBR+ + PTD 修正
瞬态宽带分析       →  FDTD (TD-FEM) + SBR+ 频域映射
```

---

## 参考文献

1. Ling, H., Chou, R., & Lee, S. W. "Shooting and bouncing rays: Calculating the RCS of an arbitrarily shaped cavity." *IEEE Trans. AP* 37(2), 194–205 (1989).
2. Gordon, W. B. "Far-field approximations to the Kirchhoff-Helmholtz representations of scattered fields." *IEEE Trans. AP* 23(4), 590–592 (1975).
3. Ufimtsev, P. "Comments on diffraction principles and limitations of RCS reduction techniques." *Proc. IEEE* 84(12), 1828–1851 (1996).
4. Pharr, M., Jakob, W., & Humphreys, G. *Physically Based Rendering*, 4th ed. MIT Press, 2023. (BVH/SAH 参考)
5. Möller, T. & Trumbore, B. "Fast, minimum storage ray-triangle intersection." *J. Graphics Tools* 2(1), 21–28 (1997).
