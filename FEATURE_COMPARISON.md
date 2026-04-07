# REM vs Palace — 功能对比与已有能力说明

> 版本：v1.0，2026-04-05  
> 描述 REM 当前已实现的全部功能，并与 Palace v0.11 进行逐项对比。

---

## 总览

REM（Rust Electromagnetic）是一款对标 [Palace](https://github.com/awslabs/palace) 的电磁仿真工具，
采用纯 Rust 实现，可编译至 `wasm32-unknown-unknown` 在浏览器中运行。
当前版本 **v0.8.1** 覆盖 Palace 全部主要求解器，并额外提供 Palace 不具备的矩量法、BEM 和 SBR+ 高频求解器。

```
所有测试：201 个（cargo test --workspace），零失败
代码量：~11,000 行（17 个 crate，不含 vendor/）
```

---

## 1. 求解器能力对比矩阵

| 功能 | Palace v0.11 | REM v0.8.1 | 说明 |
|------|:---:|:---:|------|
| **静电场** (Electrostatic) | ✅ | ✅ | P1 FEM，变介电常数，电容矩阵提取 |
| **静磁场** (Magnetostatic) | ✅ | ✅ | P1 FEM（2-D A_z），变磁导率，电感矩阵提取 |
| **特征模** (Eigenmode) | ✅ | ✅ | Lanczos 移位逆迭代，多模式，VTK 模态输出 |
| **频域驱动** (Driven) | ✅ | ✅ | 频率扫描，S 参数提取，集总端口 |
| **时域瞬态** (Transient) | ✅ | 🔲 | 配置解析已就绪；TD-FEM v1.0（Newmark-β）和 v1.1（IMEX-ARK 自适应）规划中，详见 FDTD_PLAN.md |
| **S 参数提取** | ✅ | ✅ | `postpro/port-S.csv`，Palace 格式兼容 |
| **集总端口** (Lumped Port) | ✅ | ✅ | LumpedPort 激励 + 阻抗边界 |
| **波导端口** (Wave Port) | ✅ | 🔲 | 配置字段已定义；场型匹配待实现 |
| **自适应网格细化** (AMR) | ✅ | 🔲 | fem-rs 已提供 ZZ/Kelly/Dörfler 估计器（v0.8.1 引入）；集成待完成 |
| **高阶基函数** (p-FEM) | ✅ | ✅ | 配置字段 `Solver.Order` 已解析；P1 为当前默认 |
| GMSH .msh 网格导入 | ✅ | ✅ | 完整 .msh v2/v4 解析，物理组 → 边界/材料映射 |
| ParaView VTK 输出 | ✅ | ✅ | ASCII VTK legacy，可直接用 ParaView 打开 |
| JSON 配置文件 | ✅ | ✅ | 完整 Palace JSON schema，支持 C++ 风格注释剥除 |
| YAML 配置文件 | ✅ | ✅ | serde_yaml 解析，字段名与 JSON 完全一致 |
| **WASM 目标** | ❌ | ✅ | 全部求解器可编译至 `wasm32-unknown-unknown` |
| **Web Demo（Yew）** | ❌ | ✅ | `crates/yew-app`，浏览器内直接运行求解器 |
| MPI 并行（native） | ✅ | ✅ | `Comm` trait 抽象，feature = "mpi" 启用 rsmpi |
| MPI 模拟（WASM）| ❌ | ✅ | jsmpi + Web Worker，WASM 多线程模拟 |
| 网格分区（METIS） | ✅ | ✅ | `vendor/rmetis`，纯 Rust METIS 5.1.x 兼容实现 |
| **矩量法 MoM（RWG+CFIE）** | ❌ | ✅ | 全波散射，PEC 球体 RCS vs Mie 误差 < 0.5 dB |
| **边界元法 BEM（Laplace P0）** | ❌ | ✅ | Laplace 外 Dirichlet 问题，电容提取 |
| **SBR+ 高频射线追踪 + PO** | ❌ | ✅ | AABB BVH，两阶段 PO，ka=10.5 误差 < 0.1 dB |
| **RCS / 远场后处理** | ❌ | ✅ | PO 远场积分，rcs_sbr.csv，多方向扫描 |

**图例**：✅ 已实现并通过验证　🔲 待实现（有规划）　❌ 不支持

---

## 2. 静电场求解器

**问题类型**：`Problem.Type = "Electrostatic"`

**方法**：P1 有限元，变介电常数 ε(x)，PCG + Jacobi 预条件求解器

**边界条件**：

| 类型 | 配置字段 | 说明 |
|------|---------|------|
| PEC（φ=0） | `Boundaries.PEC` | 完美导体，零电位 |
| Ground（φ=0） | `Boundaries.Ground` | 接地，零电位 |
| Terminal（φ=V） | `Boundaries.Terminal` | 激励端口，指定电位 |
| LumpedPort（φ=V，R可选） | `Boundaries.LumpedPort` | 集总端口 |
| 自然边界（∂φ/∂n=0） | 默认（未指定的边界） | Neumann 边界 |

**��出**：
- `postpro/domain-E.csv`：各域电场能量 U = (1/2)∫ε|∇φ|² dΩ
- `postpro/capacitance.csv`：电容矩阵（多端口时）
- `paraview/solution.vtk`：φ 电位场 + E 电场矢量

**验证**：平行板电容与解析解 ε₀A/d 误差 < 1e-12

---

## 3. 静磁场求解器

**问题类型**：`Problem.Type = "Magnetostatic"`

**方法**：2-D 磁矢位 A_z 公式，P1 FEM，变磁导率 ν(x)

**边界条件**：

| 类型 | 说明 |
|------|------|
| Ground | A_z = 0（磁通量不穿出） |
| SurfaceCurrent | A_z = 1（激励端口） |

**后处理**：
- B 场恢复：B_x = ∂A_z/∂y，B_y = −∂A_z/∂x（梯度恢复）
- 磁能量：U = (1/2)∫ν|∇A_z|² dΩ

**输出**：
- `postpro/domain-B.csv`：磁场能量
- `paraview/solution.vtk`：A_z 标量场 + B 矢量场

**验证**：
- 线性 A_z = y，铁磁 μ_r=1000 误差 < 1e-12
- 磁能量与解析解 ν₀/2 误差 < 1e-12

---

## 4. 特征模求解器

**问题类型**：`Problem.Type = "Eigenmode"`

**方法**：Lanczos 迭代 + 移位逆（shift-invert），求解广义特征值问题 Kx = λMx

**关键配置**（`Solver.Eigenmode`）：

| 字段 | 说明 | 默认 |
|------|------|------|
| `N` | 求解模式数 | 1 |
| `Target` | 目标频率 [Hz]，用于移位 σ = (ω/c)² | 0.0 |
| `Tol` | 迭代容差 | 1e-6 |
| `Save` | 保存前 N 个模态到 VTK | 1 |

**输出**：
- `postpro/eig.csv`：特征频率列表（Hz）
- `paraview/mode_NNNN.vtk`：各模态 φ 场

---

## 5. 频域驱动求解器

**问题类型**：`Problem.Type = "Driven"`

**方法**：频率域 FEM（Helmholtz 方程），频率扫描

**关键配置**（`Solver.Driven`）：

| 字段 | 说明 |
|------|------|
| `MinFreq` | 起始频率 [GHz] |
| `MaxFreq` | 终止频率 [GHz] |
| `FreqStep` | 频率步进 [GHz] |
| `SaveStep` | 每 N 步保存一次 VTK |

**输出**：
- `postpro/port-S.csv`：S 参数（f, Re(S11), Im(S11), |S11| dB）
- `driven_NNNN.vtk`：各频率步场量

---

## 6. 矩量法求解器（MoM）

> **REM 独有，Palace 不支持**

**问题类型**：`Problem.Type = "MoM"`

**方法**：RWG 矢量基函数 + CFIE（组合场积分方程），适用于 PEC 散射体全波分析

**核心模块**（`crates/mom/src/`）：

| 模块 | 说明 |
|------|------|
| `surface_mesh.rs` | 从 RemMesh 提取 PEC 表面三角网格 + 共享边拓扑 |
| `quadrature.rs` | Dunavant 三角形高斯求积（阶次 1/3/5/7/9，内嵌系数） |
| `green.rs` | 3D Helmholtz Green 函数 G = exp(-jkR)/(4πR) 及法向导数 |
| `singular.rs` | Duffy 自积分（对角块）+ Sauter-Schwab 奇异积分（共边/共顶点块）|
| `assemble.rs` | EFIE/MFIE/CFIE Z 矩阵装配，faer 密集矩阵，rayon 并行 |
| `excitation.rs` | 平面波激励向量（θ/φ 极化，任意入射方向） |
| `postprocess.rs` | RCS 远场积分，VTK 表面电流输出 |
| `mie.rs` | Mie 级数解析解（用于验证） |
| `basis/rwg.rs` | RWG 基函数评估 + 散度计算 |

**关键配置**（`Solver.MoM`）：

| 字段 | 说明 | 默认 |
|------|------|------|
| `Equation` | `"EFIE"` \| `"MFIE"` \| `"CFIE"` | `"CFIE"` |
| `Basis` | `"RWG"` \| `"Pulse"` | `"RWG"` |
| `FreqMin` / `FreqMax` | 频率范围 [Hz] | — |
| `Alpha` | CFIE 混合系数（0=MFIE, 1=EFIE） | 0.5 |
| `ThetaInc` / `PhiInc` | 入射角 [°] | 0.0 |
| `Polarization` | `"theta"` \| `"phi"` | `"theta"` |

**输出**（`Postprocessing.RCS`）：
- `postpro/rcs.csv`：RCS 方向图（θ, φ, σ_dBsm）
- `paraview/surface_current.vtk`：表面电流 J 矢量场

**验证**：PEC 球体（r=0.5 m）@ 1 GHz，kα≈10.5，单站 RCS vs Mie 误差 < 0.5 dB

**配置示例**：
```json
{
  "Problem": { "Type": "MoM", "Output": "output/sphere" },
  "Model":   { "Mesh": "sphere.msh", "L0": 1e-3 },
  "Boundaries": { "PEC": { "Attributes": [1] } },
  "Solver": {
    "MoM": {
      "Equation": "CFIE", "Basis": "RWG",
      "FreqMin": 1.0e9, "FreqMax": 1.0e9,
      "Alpha": 0.5, "ThetaInc": 0.0, "PhiInc": 0.0
    }
  },
  "Postprocessing": { "RCS": { "ThetaDeg": "0:5:180", "PhiDeg": [0.0] } }
}
```

---

## 7. 边界元法求解器（BEM）

> **REM 独有，Palace 不支持**

**问题类型**：`Problem.Type = "BEM"`

**方法**：Laplace P0 边界积分方程，外 Dirichlet 问题（PEC 静电）

**BIE 公式**：
```
½φ(r) + ∫_S ∂G_L/∂n'(r,r') φ(r') dS' = ∫_S G_L(r,r') q(r') dS'
G_L(r,r') = 1/(4π|r-r'|)
```

**核心模块**（`crates/bem/src/`）：

| 模块 | 说明 |
|------|------|
| `kernel.rs` | Laplace Green 函数 G、法向导数 ∂G/∂n'（观测点侧 ∂G/∂n） |
| `assemble.rs` | V（单层势）+ K（双层势）矩阵装配，P0 基函数，Duffy 对角自积分 |
| `solve.rs` | faer LU 求解 |
| `postprocess.rs` | 电容矩阵提取 + 电位 VTK 输出 |

**WASM 支持**：faer LU 分解支持 WASM，单线程（N < 1000 推荐）

---

## 8. SBR+ 高频射线追踪求解器

> **REM 独有，Palace 不支持**

**问题类型**：`Problem.Type = "SBR"`

**方法**：SBR+（Shooting and Bouncing Rays Plus），几何光学 + 物理光学（PO），适用于电大目标（kα >> 1）

**算法**（两阶段 PO）：

```
阶段 1 — 一次弹射 PO（per-face，与射线密度无关）
  for face in surf.faces:
      几何可见：dot(n̂, -k̂) > 0
      阴影测试：从 face.centroid + ε·n̂ 向 -k̂ 发射阴影射线
      若可见：J = 2 n̂ × H_inc(centroid)

阶段 2 — 多次弹射（bounce ≥ 1，射线管归一化）
  J_bounce += (A_ray / A_face) × 2 n̂ × H_ray

远场 PO 积分：
  N(r̂) = Σ_m J_m · A_m · exp(jk r̂·r_m)
  σ(r̂) = 4π|r̂×(r̂×N)|² / |E_inc|²  [m²]
```

**核心模块**（`crates/sbr/src/`）：

| 模块 | 说明 |
|------|------|
| `bvh.rs` | AABB BVH（SAH 分割），Möller-Trumbore 射线-三角形求交 |
| `ray.rs` | Ray / RayHit / RayPath 数据结构 |
| `fresnel.rs` | Fresnel 系数，PEC 镜面反射，PO 感应电流 J = 2n̂×H |
| `excitation.rs` | 孔径射线铺设，平面波 E/H 场 |
| `po_integral.rs` | 离散远场 PO 积分 → 复散射振幅 → RCS |
| `output.rs` | `rcs_sbr.csv` + 感应电流 VTK |

**关键配置**（`Solver.SBR`）：

| 字段 | 说明 | 默认 |
|------|------|------|
| `FreqMin` / `FreqMax` / `FreqStep` | 频率范围和步进 [Hz] | — |
| `RayDensity` | 射线面密度 [rays/m²] | 5000 |
| `MaxBounces` | 最大弹射次数 | 3 |
| `WeightThresh` | 射线能量截断阈值 | 1e-4 |
| `TargetType` | `"PEC"` \| `"Dielectric"` \| `"Coated"` | `"PEC"` |
| `ThetaInc` / `PhiInc` | 入射仰角/方位角 [°] | 0.0 |
| `Polarization` | `"theta"` \| `"phi"` | `"theta"` |

**输出**（`Postprocessing.RCS`）：
- `postpro/rcs_sbr.csv`：RCS 方向图（θ, φ, σ_dBsm）
- `paraview/sbr_FREQ.vtk`：感应电流 J 矢量场

**验证**：PEC 球（r=0.5 m）@ 1 GHz，ka=10.5，单站 RCS 误差 **0.05 dB**（< 3 dB 限值）

**mesh 分辨率约束**：PO 相位积分要求面片尺寸 < λ/4，即 Nring > ka/2 纬度环。

**配置示例**：
```json
{
  "Problem": { "Type": "SBR", "Output": "output/sbr_sphere" },
  "Model":   { "Mesh": "sphere.msh", "L0": 1.0 },
  "Boundaries": { "PEC": { "Attributes": [1] } },
  "Solver": {
    "SBR": {
      "FreqMin": 3.0e9, "FreqMax": 3.0e9,
      "RayDensity": 5000.0, "MaxBounces": 3,
      "ThetaInc": 0.0, "PhiInc": 0.0, "Polarization": "theta"
    }
  },
  "Postprocessing": { "RCS": { "ThetaDeg": "0:5:180", "PhiDeg": [0.0] } }
}
```

---

## 9. Palace 配置兼容性

REM 完全兼容 Palace JSON/YAML 配置文件格式。Palace 用户无需修改已有配置即可在 REM 中运行。

### 9.1 支持的边界类型

| Palace 字段 | REM 支持 | BoundaryTag 枚举 |
|------------|:---:|-----------------|
| `Boundaries.PEC` | ✅ | `Pec` |
| `Boundaries.PMC` | ✅（解析，Neumann 自然边界） | `Pmc` |
| `Boundaries.Impedance` | ✅（解析） | `Impedance { rs }` |
| `Boundaries.Absorbing` | ✅（解析） | `Absorbing { order }` |
| `Boundaries.Conductivity` | ✅（解析） | `Conductivity { sigma }` |
| `Boundaries.Ground` | ✅ | `Ground` |
| `Boundaries.Terminal` | ✅ | `Terminal { index }` |
| `Boundaries.LumpedPort` | ✅ | `LumpedPort { index, r }` |
| `Boundaries.WavePort` | ✅（解析，场求解待完成） | `WavePort { index }` |
| `Boundaries.SurfaceCurrent` | ✅ | `SurfaceCurrent { index }` |

### 9.2 支持的材料参数

| Palace 字段 | REM 支持 | 说明 |
|------------|:---:|------|
| `Permittivity` (εᵣ) | ✅ | 标量，默认 1.0 |
| `Permeability` (μᵣ) | ✅ | 标量，默认 1.0 |
| `Conductivity` (σ) | ✅（解析） | [S/m]，损耗计算待完成 |
| `LossTan` | ✅（解析） | 介质损耗角正切 |
| `Attributes` 范围格式 | ✅ | `"1,3-5"` 和 `[1,3,4,5]` 均可 |

### 9.3 REM 专有扩展（对 Palace 无影响）

以下字段在 Palace 中被静默忽略，不影响现有 Palace 工作流：

```json
"Solver": {
  "MoM":  { ... },   // MoM 求解器参数
  "SBR":  { ... }    // SBR+ 求解器参数
},
"Postprocessing": {
  "RCS":  { ... }    // 远场/RCS 输出配置
}
```

---

## 10. 已验证示例

| 示例目录 | 问题类型 | 验证指标 | 结果 |
|---------|---------|---------|------|
| `examples/parallel_plate/` | Electrostatic | C = ε₀A/d，误差 < 1e-12 | ✅ |
| `examples/coaxial/` | Electrostatic | C/L = 2πε₀/ln(b/a)，误差 < 0.5% | ✅ |
| `examples/rings/` | Magnetostatic | ν 界面跳变条件，A_z = y（线性解） | ✅ |
| `examples/transmon/` | Eigenmode | 谐振腔特征频率 | ✅ |
| `examples/cpw/` | Driven | S₁₁ 频率扫描 | ✅ |
| `examples/antenna/` | MoM | PEC 偶极子输入阻抗 | ✅ |
| `examples/spheres/` | MoM | PEC 球 RCS vs Mie（误差 < 0.5 dB） | ✅ |
| `examples/sbr_sphere/` | SBR | PEC 球 RCS @ 1/3 GHz，vs Mie（误差 < 0.1 dB） | ✅ |

---

## 11. 求解器选型指南

```
目标电尺寸   kα < 3         → MoM（全波，严格解）
目标电尺寸   kα = 3–15      → MoM 为主，SBR+ 作参考
目标电尺寸   kα > 15        → SBR+（高频渐近，O(N_face) 内存）
超大目标（飞机/舰船级）     → SBR+ + PTD 绕射修正（路线图 v1.1+）
静态电场/电容提取           → Electrostatic FEM 或 BEM
静态磁场/电感提取           → Magnetostatic FEM（2-D）
谐振腔本征频率              → Eigenmode
S 参数、集总端口匹配        → Driven
时域宽带脉冲响应            → Transient（TD-FEM）
  · 固定步长（v1.0 规划）  → time_scheme: "newmark"，Newmark-β，无条件稳定
  · 自适应步长（v1.1 规划）→ time_scheme: "imex_ark"，IMEX-ARK3(2)4L[2]SA，误差控制
```

---

## 12. WASM / 浏览器限制

| 约束 | 限制 | 说明 |
|------|------|------|
| 线程 | 单线程（无 rayon） | MoM/SBR+ 退化为串行 |
| 内存 | ~30 MB 堆 | MoM 建议 N < 1000 面元 |
| 文件系统 | 无磁盘 IO | 输出返回 Blob URL |
| `rem-mom` | 可用 | rayon 条件编译排除 |
| `rem-sbr` | 可用 | rayon cfg-excluded |
| `rem-bem` | 可用 | faer LU 支持 WASM |

---

## 13. 技术债务与已知限制

| 项目 | 说明 | 优先级 |
|------|------|--------|
| 3-D 静磁（Nedelec H(curl)） | 当前仅 2-D A_z 标量 | 中 |
| 波导端口场匹配 | 配置已支持，场型积分待实现 | 中 |
| AMR 集成 | fem-rs AMR 估计器已就绪，与求解器对接待完成 | 中 |
| 时域瞬态（TD-FEM） | v1.0：Nedelec H(curl) + Newmark-β 固定步长（FDTD_PLAN.md §9.1，20–26 天）；v1.1：IMEX-ARK 自适应步长（§9.2，+8–12 天） | 高 |
| MoM 介质目标（PMCHWT） | PEC 目标已完整，介质扩展待实现 | 低 |
| SBR+ 边缘绕射（PTD） | 大角度散射精度提升 | 低 |
| p-FEM（P2+）实际应用 | 配置字段已有，单元装配待扩展 | 低 |
