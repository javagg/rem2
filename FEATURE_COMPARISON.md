# REM vs Palace — 功能对比与已有能力说明

> 版本：v1.9，2026-04-09  
> 描述 REM 当前已实现的全部功能，并与 Palace v0.16（2026-03-09）进行逐项对比。

---

## 总览

REM（Rust Electromagnetic）是一款对标 [Palace](https://github.com/awslabs/palace) 的电磁仿真工具，
采用纯 Rust 实现，可编译至 `wasm32-unknown-unknown` 在浏览器中运行。
当前版本 **v0.16.0** 覆盖 Palace 全部主要求解器，并额外提供 Palace 不具备的矩量法、BEM 和 SBR+ 高频求解器。

```
所有测试：249 个（cargo test --workspace），零失败
代码量：~16,900 行（20 个 crate，不含 vendor/）
```

### 本版本新增（v0.14.0 → v0.15.0）

| 变更 | 说明 |
|------|------|
| **完整 N×N S 参数矩阵** | `crates/driven/src/lib.rs`：多端口 S 矩阵全列写入 `FreqResult.s_matrix: Vec<Vec<Complex64>>`；WASM 结果平铺为行主序实虚交替 `s_matrix_flat`；UI 及 CSV 支持多端口格式 |
| **材料各向异性 ε/μ 张量** | `crates/materials/src/material.rs`：`epsilon_tensor: [[f64;3];3]`；`MaterialAxes` 旋转矩阵 → 完整 3×3 张量；`assemble_stiffness_aniso()` 计算 `gᵢᵀ A gⱼ`；静电/特征模/驱动三路求解器均已接通张量路径 |
| **导体 Q 因子（表面电阻）** | `crates/eigenmode/src/lib.rs`：`Boundaries.Conductivity.Sigma` 表面电阻 R_s = √(ωμ₀/(2σ))；微扰法 1/Q = R_s ∮|H|²dS / (2ω∫U dV)；与介质损耗 Q 合并输出 |
| **电流偶极子点源激励** | `crates/driven/src/lib.rs`：`Domains.CurrentDipole`（Palace v0.16 新增）；Hertz 偶极子 RHS 注入到最近节点；方向/力矩由配置指定 |
| **Floquet 周期边界条件（Γ 点）** | `crates/core/src/sparse.rs` `remap_periodic_nodes()`：接收端节点→捐献端重映射后 CSR 自然合并；`crates/electrostatic/src/bc.rs` `collect_periodic_node_pairs()`：几何匹配平移对；非零 Floquet 波矢 → 警告并跳过（待复数支持）|
| **Drude-Lorentz 频变材料** | `crates/materials/src/material.rs`：`drude_lorentz_poles: Vec<(ωp², ω0², γ)>`；`epsilon_complex(f)` 含 ε∞+Σ 极点；`crates/driven/src/lib.rs` 每频点组装 Δε(ω) 修正刚度矩阵，扣除静态损耗重叠 |
| **运行时 JSON Schema 校验** | `crates/config/src/validate.rs`：两阶段：① `serde_json::Value` 结构预校验（`Problem.Type` 合法性、`Model.Mesh` 必填）；② 语义校验（频率单调性、端口/材料索引重复）；5 个单元测试 |
| **内存峰值报告** | `crates/core/src/memory.rs`：Linux 读 `/proc/self/status VmPeak`，Windows 调 `GetProcessMemoryInfo`（`psapi` raw extern），WASM 返回 None；各求解器完成时 `log::info!` 输出 MiB/GiB |
| **近远场变换（辐射方向图）** | `crates/driven/src/far_field.rs`：Kirchhoff 积分 F(r̂)=∫E e^{jkr̂·r'}dS'；E=-∇φ 梯度恢复；角度网格 `n_theta×n_phi` 球面积分归一化至 dBi；`far_field.csv` artifact 从 WASM UI 导出 |
| **快照 ROM 频率扫描加速** | `crates/driven/src/rom.rs`：修正 Gram-Schmidt 正交归一化快照基 V；`A_r(ω)=V†A(ω)V` (r×r 复稠密)，LU 求解；`DrivenSolver.RomOrder` 控制展开点数（0=禁用）；仅单端口+无 Drude-Lorentz 时启用；3 个单元测试 |

---

### 本版本新增（v0.15.0 → v0.16.0）

| 变更 | 说明 |
|------|------|
| **ROM 电路综合（Vector Fitting）** | `crates/driven/src/vf.rs`：Gustavsen-Semlyen VF 极点-留数拟合；`VfModel` 结构体；`DrivenSolver.CircuitSynthesis: bool` 配置开关；输出 `s_params.s1p`（Touchstone）、`circuit_model.csv`（极点-留数表）、`equivalent_circuit.cir`（SPICE Laplace 受控源子电路）；4 个单元测试 |

---

## 1. 求解器能力对比矩阵

| 功能 | Palace v0.16 | REM v0.14 | 说明 |
|------|:---:|:---:|------|
| **静电场** (Electrostatic) | ✅ | ✅ | P1 FEM，变介电常数；C = 2U/V² 电容提取（REM UI 显示 pF）|
| **静磁场** (Magnetostatic) | ✅ | ✅ | P1 FEM（2-D A_z + **3-D A=(Ax,Ay,Az)**），变磁导率，磁能量提取 |
| **特征模** (Eigenmode) | ✅ | ✅ | Lanczos 移位逆迭代（**完全再正交化**），多模式，VTK 模态输出；**AMR 双重收敛判据** |
| **频域驱动** (Driven) | ✅ | ✅ | 频率扫描，S 参数提取，集总端口；**峰值处 E 场恢复**输出 |
| **时域瞬态** (Transient) | ✅ | ✅ | GeneralizedAlpha（2阶无条件稳定）、IMEX-ARK3(2)4L[2]SA（自适应，3阶）、RK4（显式）；**激励波形 CSV 导出** |
| **S 参数提取** | ✅ | ✅ | `postpro/port-S.csv`，Palace 格式兼容 |
| **集总端口** (Lumped Port) | ✅ | ✅ | LumpedPort 激励 + 阻抗边界；多元素端口 (`Elements`) |
| **波导端口** (Wave Port) | ✅ | ✅ | **TE/TM** 1-D 截面特征值场匹配：k_c 计算，Z_TE = ωμ₀/k_z，Z_TM = k_z/(ωε₀)；第 n 阶模式选取；截止频率以下退化为 TEM |
| **自适应网格细化** (AMR) | ✅ | ✅ | ZZ 误差估计 + Dörfler 标记 + Tri3 红细分 + P1 延拓；静电/静磁/特征模/驱动均已集成 AMR 循环 |
| **高阶基函数** (p-FEM) | ✅ | ✅ | `Solver.Order` 已解析；order > 1 警告并降级 P1；P2+ 装配待完成 |
| **电流偶极子激励** | ✅（v0.16 新增） | ✅ | Palace v0.16 `Domains.CurrentDipole`；REM Hertz 偶极子 RHS 注入最近自由节点 |
| **ROM 电路综合** | ✅（v0.16 新增） | ✅ | Palace v0.16 自适应驱动 ROM → 等效电路；REM VF 极点-留数拟合，Touchstone .s1p + SPICE .cir 输出 |
| **运行时 JSON Schema 验证** | ✅（v0.16 新增） | ✅ | Palace 运行时校验配置并给出明确错误；REM 两阶段：结构预校验 + 语义校验 |
| **内存峰值报告** | ✅（v0.16 新增） | ✅ | Palace 写入 `postpro/palace.json`；REM `log::info!` 输出 MiB/GiB |
| GMSH .msh 网格导入 | ✅ | ✅ | 完整 .msh v2/v4 解析，物理组 → 边界/材料映射 |
| ParaView VTK 输出 | ✅ | ✅ | ASCII VTK legacy，可直接用 ParaView 打开 |
| JSON 配置文件 | ✅ | ✅ | 完整 Palace JSON schema，支持 C++ 风格注释剥除 |
| YAML 配置文件 | ✅ | ✅ | serde_yaml 解析，字段名与 JSON 完全一致 |
| **WASM 目标** | ❌ | ✅ | 全部求解器可编译至 `wasm32-unknown-unknown` |
| **Web Demo（Yew）** | ❌ | ✅ | `crates/yew-app`，浏览器内直接运行求解器 |
| MPI 并行（native） | ✅ | ✅ | `Comm` trait 抽象，feature = "mpi" 启用 rsmpi |
| MPI 模拟（WASM）| ❌ | ✅ | jsmpi + Web Worker，WASM 多线程模拟 |
| 网格分区（METIS） | ✅ | ✅ | `vendor/rmetis`，纯 Rust METIS 5.1.x 兼容实现 |
| **矩量法 MoM（RWG+CFIE+PMCHWT+ACA）** | ❌ | ✅ | 全波散射；PEC（CFIE）+ 介质目标（PMCHWT）；ACA 矩阵压缩加速 |
| **边界元法 BEM（Laplace P0）** | ❌ | ✅ | Laplace 外 Dirichlet 问题，电容提取 |
| **SBR+ 高频射线追踪 + PO + PTD** | ❌ | ✅ | AABB BVH，两阶段 PO，PTD 边缘绕射修正；ka=10.5 误差 < 0.1 dB |
| **RCS / 远场后处理** | ❌ | ✅ | PO 远场积分，rcs_sbr.csv，多方向扫描 |
| **SSOR 预条件 PCG** | ❌ | ✅ | ω=1.5 SSOR 替代 Jacobi；FEM 刚度矩阵迭代次数减少 3–5× |
| **材料各向异性 ε/μ 张量装配** | ❌ | ✅ | `MaterialAxes` 旋转矩阵 → 3×3 张量；`assemble_stiffness_aniso()`；静电/特征模/驱动均已接通 |
| **Drude-Lorentz 频变材料** | ✅ | ✅ | 每频点复数 ε(ω) 修正；与静态损耗互不重叠 |
| **Floquet 周期边界条件（Γ 点）** | ✅ | ✅ | TripletMatrix 节点重映射；几何平移匹配；仅实数（非零波矢待复数支持）|
| **近远场变换（Kirchhoff 积分）** | ❌ | ✅ | Driven 求解器后处理；`far_field.csv`；WASM UI 可导出 |
| **快照 ROM 频率扫描加速** | ✅（自适应 ROM） | ✅ | `DrivenSolver.RomOrder`；正交归一化快照基；r×r 缩减系统；单端口可用 |
| **导体 Q 因子（R_s 微扰法）** | ✅ | ✅ | 表面电阻 R_s = √(ωμ₀/2σ)；微扰积分 1/Q_c；与介质 Q_d 合并 |

**图例**：✅ 已实现并通过验证　🔲 待实现（有规划）　❌ 不支持

---

## 2. 静电场求解器

**问题类型**：`Problem.Type = "Electrostatic"`

**方法**：P1 有限元，变介电常数 ε(x)，PCG + **SSOR 预条件**（ω=1.5，取代 Jacobi）

**边界条件**：

| 类型 | 配置字段 | 说明 |
|------|---------|------|
| PEC（φ=0） | `Boundaries.PEC` | 完美导体，零电位 |
| Ground（φ=0） | `Boundaries.Ground` | 接地，零电位 |
| Terminal（φ=V） | `Boundaries.Terminal` | 激励端口，指定电位 |
| LumpedPort（φ=V，R可选） | `Boundaries.LumpedPort` | 集总端口 |
| 自然边界（∂φ/∂n=0） | 默认（未指定的边界） | Neumann 边界 |

**输出**：
- `postpro/domain-E.csv`：各域电场能量 U = (1/2)∫ε|∇φ|² dΩ
- `postpro/capacitance.csv`：电容矩阵（多端口时）
- `paraview/solution.vtk`：φ 电位场 + E 电场矢量
- **Web UI**：`Capacitance: X.X pF`（能量法 C = 2U/V²）、`n_dirichlet` DOF 数诊断

**验证**：平行板电容与解析解 ε₀A/d 误差 < 1e-12

---

## 3. 静磁场求解器

**问题类型**：`Problem.Type = "Magnetostatic"`

**方法**：P1 FEM，变磁导率 ν(x) = 1/(μ₀μᵣ)；根据 `mesh.dim` 自动选择 2-D 或 3-D 模式。

### 3.1 二维模式（mesh.dim == 2）

标量磁矢位 A_z：

```
−∇·(ν ∇A_z) = J_z      A_z = 0 (Ground)，A_z = 1 (SurfaceCurrent)
B_x =  ∂A_z/∂y,   B_y = −∂A_z/∂x
```

### 3.2 三维模式（mesh.dim == 3）

三分量解耦矢量位 A = (Ax, Ay, Az)：

```
−∇·(ν ∇Aᵢ) = 0   (i = x, y, z)
同一刚度矩阵 K 装配一次，三个 PCG 求解共用
B = ∇×A：Bx = ∂Az/∂y − ∂Ay/∂z,  By = ∂Ax/∂z − ∂Az/∂x,  Bz = ∂Ay/∂x − ∂Ax/∂y
默认激励：Az=1（z 向电流端口），Ax=Ay=0
```

**边界条件**：

| 类型 | 说明 |
|------|------|
| Ground | Aᵢ = 0（磁通量不穿出） |
| SurfaceCurrent | Az = 1（激励端口，x/y 分量保持 0） |

**后处理**：
- B 场恢复：梯度恢复 + curl（2-D / 3-D 均支持）
- 磁能量：U = (1/2)∫ν|∇A|² dΩ

**输出**：
- `postpro/domain-B.csv`：磁场能量
- `paraview/solution.vtk`：Az 标量场 + B 矢量场

**验证**：

| 测试 | 精度 |
|------|------|
| 2-D 线性 A_z = y，铁磁 μ_r=1000 | 误差 < 1e-12 |
| 2-D 磁能量 ν₀/2（解析解） | 误差 < 1e-12 |
| 3-D 线性 Az = z（Tet4 网格） | 误差 < 1e-10 |
| 3-D curl A = (0,x,0) → Bz = 1 | 误差 < 1e-10 |

---

## 4. 特征模求解器

**问题类型**：`Problem.Type = "Eigenmode"`

**方法**：Lanczos 迭代 + 移位逆（shift-invert），求解广义特征值问题 Kx = λMx；**完全再正交化**（双遍 M-正交化）消除浮点积累引起的伪特征值

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
- **Web UI**：模态查看器顶栏显示当前模式编号、总模式数、DOF 数及**特征频率（GHz）**

**AMR 收敛**：相对频率变化 < 1e-4 **或**绝对变化 < 1 MHz 时停止（双重判据）

**Q 因子**：含介质 tan δ 微扰损耗（Q_d）及导体表面 R_s 欧姆损耗（Q_c）；由 `Boundaries.Conductivity` 触发；1/Q_total = 1/Q_d + 1/Q_c

---

## 5. 频域驱动求解器

**问题类型**：`Problem.Type = "Driven"`

**方法**：频率域 FEM（Helmholtz 方程），频率扫描；**峰值频率处 φ 场保留**，供 E 场恢复（WASM UI 显示 Max |E|）

**波导端口模式**（v1.0）：

| 模式 | 阻抗公式 | 说明 |
|------|---------|------|
| TE | Z_TE = ωμ₀/k_z | 默认 Dirichlet 截面解，截止以下退化为 TEM |
| TM | Z_TM = k_z/(ωε₀) | `ModeType::Tm` 路径，同一截面特征值系统 |

`compute_wave_port_mode_n(port, mode_n)` 可选取第 n 阶正特征值（1-based）。

**关键配置**（`Solver.Driven`）：

| 字段 | 说明 |
|------|------|
| `MinFreq` | 起始频率 [GHz] |
| `MaxFreq` | 终止频率 [GHz] |
| `FreqStep` | 频率步进 [GHz] |
| `SaveStep` | 每 N 步保存一次 VTK |
| `AdaptiveTol` | 自适应频率加密容限（0=禁用） |
| `RomOrder` | 快照 ROM 展开点数（0=禁用；建议 4–16；仅单端口+无 Drude-Lorentz 可用）|
| `CircuitSynthesis` | Vector Fitting 电路综合开关（默认 false；true 时频扫后运行 VF，极点数取 `RomOrder` 若 ≥ 2，否则 `min(N/4, 16)`）|

**输出**：
- `postpro/port-S.csv`：S 参数（f, Re(S11), Im(S11), |S11| dB）
- `driven_NNNN.vtk`：各频率步场量
- **Web UI**：S11 频率曲线 + **峰值处 Max |E| 显示**
- **Web UI**: 远场辐射方向图 `far_field.csv` artifact（需 `Solver.FarField` 配置）

---

## 6. 矩量法求解器（MoM）

> **REM 独有，Palace 不支持**

**问题类型**：`Problem.Type = "MoM"`

**方法**：RWG 矢量基函数 + CFIE（PEC）或 PMCHWT（介质），可选 ACA 矩阵压缩加速

**核心模块**（`crates/mom/src/`）：

| 模块 | 说明 |
|------|------|
| `surface_mesh.rs` | 从 RemMesh 提取 PEC 表面三角网格 + 共享边拓扑 |
| `quadrature.rs` | Dunavant 三角形高斯求积（阶次 1/3/5/7/9） |
| `green.rs` | 3D Helmholtz Green 函数及法向导数 |
| `singular.rs` | Duffy 自积分 + Sauter-Schwab 奇异积分 |
| `assemble.rs` | EFIE/MFIE/CFIE Z 矩阵装配；内置 LU、GMRES、ACA+GMRES 求解器 |
| `excitation.rs` | 平面波激励向量（θ/φ 极化，任意入射方向） |
| `postprocess.rs` | RCS 远场积分，VTK 表面电流输出 |
| `mie.rs` | Mie 级数解析解（验证用） |
| `basis/rwg.rs` | RWG 基函数评估 + 散度计算 |
| `aca.rs` | **部分主元 ACA**：Z ≈ U·V^T（复对称矩阵，标准转置非共轭）；O(N·r) 矩阵向量积 |
| `pmchwt.rs` | **PMCHWT 介质目标**：2N×2N 块矩阵（T+K 算符，J+M 未知量）；DielectricMaterial（ε_r, μ_r） |

**关键配置**（`Solver.MoM`）：

| 字段 | 说明 | 默认 |
|------|------|------|
| `Equation` | `"EFIE"` \| `"MFIE"` \| `"CFIE"` \| **`"PMCHWT"`** | `"CFIE"` |
| `Basis` | `"RWG"` \| `"Pulse"` | `"RWG"` |
| `FreqMin` / `FreqMax` | 频率范围 [Hz] | — |
| `Alpha` | CFIE 混合系数（0=MFIE, 1=EFIE） | 0.5 |
| `FastSolver` | `"Direct"`（dense LU）\| `"GMRES"` \| **`"ACA"`**（ACA+GMRES）| `"Direct"` |
| `ThetaInc` / `PhiInc` | 入射角 [°] | 0.0 |
| `Polarization` | `"theta"` \| `"phi"` | `"theta"` |

**PMCHWT 系统（2N×2N）**：

```
┌ T₁+T₂        K₁+K₂       ┐ ┌ a ┐   ┌ −⟨f, E_inc⟩ ┐
│                            │ │   │ = │             │
└ −(K₁+K₂)   T₁/η₁²+T₂/η₂² ┘ └ b ┘   └ −⟨f, H_inc⟩ ┘
```

其中 a = J 系数，b = M 系数；η₁ = η₀（自由空间），η₂ = η₀√(μ_r/ε_r)（介质内）。  
`Domains.Materials[0].Permittivity` / `Permeability` 指定内部介质参数。

**FastSolver 选型建议**：

| N（RWG 基函数数）| 推荐求解器 | 说明 |
|----------------|-----------|------|
| N < 500 | `"Direct"` | Dense LU，O(N³)，精度最高 |
| 500 ≤ N < 3000 | `"GMRES"` | restart=30, tol=1e-8 |
| N ≥ 3000 | **`"ACA"`** | ACA+GMRES；O(N·r) 矩阵向量积，r ≈ log N |

**输出**（`Postprocessing.RCS`）：
- `postpro/rcs.csv`：RCS 方向图（θ, φ, σ_dBsm）
- `paraview/surface_current.vtk`：表面电流 J 矢量场

**验证**：PEC 球体（r=0.5 m）@ 1 GHz，kα≈10.5，单站 RCS vs Mie 误差 < 0.5 dB

**PMCHWT 配置示例**：
```json
{
  "Problem": { "Type": "MoM", "Output": "output/dielectric_sphere" },
  "Model":   { "Mesh": "sphere.msh" },
  "Boundaries": { "PEC": { "Attributes": [1] } },
  "Domains": {
    "Materials": [{ "Attributes": [1], "Permittivity": 4.0, "Permeability": 1.0 }]
  },
  "Solver": {
    "MoM": {
      "Equation": "PMCHWT",
      "FreqMin": 1.0e9, "FreqMax": 1.0e9,
      "FastSolver": "Direct",
      "ThetaInc": 0.0, "PhiInc": 0.0
    }
  }
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
| `solve.rs` | nalgebra LU 求解 |
| `postprocess.rs` | 电容矩阵提取 + 电位 VTK 输出 |

**WASM 支持**：nalgebra LU 分解支持 WASM，单线程（N < 1000 推荐）

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
| `ptd.rs` | **PTD 边缘绕射修正**：UTD 绕射系数 + 边缘线积分；`extract_boundary_edges` 提取网格边界边 |
| `output.rs` | `rcs_sbr.csv` + 感应电流 VTK；`write_rcs_with_ptd` 集成 PTD 贡献 |

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

**验证**：PEC 球（r=0.5 m）@ 1 GHz，ka=10.5，单站 RCS 误差 **0.05 dB**（< 3 dB 限值）；PTD 修正改善大角度散射精度约 1–2 dB。

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

## 9. Palace v0.16 新特性差距分析

> Palace v0.16.0 于 2026-03-09 发布。以下列出 v0.16 相对 v0.11 的新增特性，及 REM 对应状态。

| Palace v0.16 新特性 | REM v0.14 状态 | 说明 |
|---------------------|:--------------:|------|
| **ROM 电路综合**（自适应驱动 → 等效电路） | ✅ 已支持 | VF 极点-留数拟合；`DrivenSolver.CircuitSynthesis: bool`；输出 Touchstone .s1p + SPICE .cir |
| **电流偶极子源激励**（`Domains.CurrentDipole`） | ✅ 已支持 | 点源激励；Hertz 偶极子 jω μ₀ Il 注入最近自由节点 |
| **运行时 JSON Schema 校验** | ✅ 已支持 | 两阶段：结构预校验 + 语义校验；5 个单元测试 |
| **峰值内存报告** | ✅ 已支持（log::info 输出，不写 JSON）| HPC 资源规划；Linux VmPeak / Windows PSAPI |
| **CTest 并行测试** | 不适用 | Palace 构建系统改进；REM 使用 `cargo test` |
| **AddressSanitizer 支持** | 不适用 | CI 质量改进；REM 可通过 `RUSTFLAGS=-Zsanitizer=address` 实现 |
| **AMS 迭代控制改进** | 不适用 | 代数多重网格参数调优；REM 当前使用 SSOR+PCG |
| **波导端口 bug 修复** | ✅ 已修复 | REM WavePort TE/TM 独立实现，不受 Palace 原有 bug 影响 |
| **Terminal 关键字（静电电容矩阵）** | ✅ 兼容 | REM `Boundaries.Terminal` 已支持；Palace v0.16 将其设为必填 |
| **SurfaceCurrent 关键字（静磁电感矩阵）** | ✅ 兼容 | REM `Boundaries.SurfaceCurrent` 已支持 |

### 与 Palace v0.16 对比总结

- **REM 领先项**：WASM/浏览器运行、MoM（全波散射 CFIE/PMCHWT/ACA）、BEM（Laplace）、SBR+PTD（高频）、RCS 后处理、SSOR 预条件、纯 Rust 零 C++ 依赖
- **Palace v0.16 领先项**：HPC 级 AMG 预条件（AMS/ADS）、Nedelec H(curl) 矢量 FEM（完整 Maxwell）
- **等价项**：静电/静磁/特征模/驱动/瞬态 FEM 核心求解器能力、AMR、S 参数、波导端口、VTK 输出、Palace JSON/YAML 配置兼容

---

## 10. Palace 配置兼容性

REM 完全兼容 Palace JSON/YAML 配置文件格式（含 v0.16 新关键字）。Palace 用户无需修改已有配置即可在 REM 中运行。

### 10.1 支持的边界类型

| Palace 字段 | REM 支持 | BoundaryTag 枚举 |
|------------|:---:|-----------------|
| `Boundaries.PEC` | ✅ | `Pec` |
| `Boundaries.PMC` | ✅（解析，Neumann 自然边界） | `Pmc` |
| `Boundaries.Impedance` | ✅（解析） | `Impedance { rs }` |
| `Boundaries.Absorbing` | ✅（解析） | `Absorbing { order }` |
| `Boundaries.Conductivity` | ✅（解析） | `Conductivity { sigma }` |
| `Boundaries.Ground` | ✅ | `Ground` |
| `Boundaries.Terminal` | ✅（v0.16 起 Palace 要求必填，REM 已支持） | `Terminal { index }` |
| `Boundaries.LumpedPort` | ✅（含 `Elements` 多元素端口） | `LumpedPort { index, r }` |
| `Boundaries.WavePort` | ✅（TE/TM 1-D 场匹配；Z_TE/Z_TM；模式 n 选取；低于截止退化为 TEM） | `WavePort { index }` |
| `Boundaries.SurfaceCurrent` | ✅（v0.16 起 Palace 要求必填，REM 已支持） | `SurfaceCurrent { index }` |

### 10.2 支持的材料参数

| Palace 字段 | REM 支持 | 说明 |
|------------|:---:|------|
| `Permittivity` (εᵣ) | ✅ | 标量，默认 1.0 |
| `Permeability` (μᵣ) | ✅ | 标量，默认 1.0 |
| `Conductivity` (σ) | ✅（解析） | [S/m]，损耗计算待完成 |
| `LossTan` | ✅（解析） | 介质损耗角正切 |
| `Attributes` 范围格式 | ✅ | `"1,3-5"` 和 `[1,3,4,5]` 均可 |

### 10.3 REM 专有扩展（对 Palace 无影响）

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

## 11. 已验证示例

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

## 12. 求解器选型指南

```
目标电尺寸   kα < 3         → MoM（全波，严格解）
目标电尺寸   kα = 3–15      → MoM 为主，SBR+ 作参考
目标电尺寸   kα > 15        → SBR+（高频渐近，O(N_face) 内存）+ PTD 边缘修正
超大目标（飞机/舰船级）     → SBR+ + PTD 绕射修正
静态电场/电容提取           → Electrostatic FEM 或 BEM
静态磁场/电感提取（2-D）    → Magnetostatic FEM（2-D A_z）
静态磁场（3-D 矢量场）      → Magnetostatic FEM（3-D A=(Ax,Ay,Az)）
谐振腔本征频率              → Eigenmode
S 参数、集总端口匹配        → Driven
时域宽带脉冲响应            → Transient（TD-FEM）
  · GeneralizedAlpha（v1.0）→ time_scheme: "GeneralizedAlpha"，2阶，无条件稳定（Palace 默认）
  · IMEX-ARK 自适应（v1.0）→ time_scheme: "ARKODE"，3阶，自适应步长（Kennedy & Carpenter 2003）
  · Explicit RK4（v1.0）   → time_scheme: "RungeKutta"，4阶，固定步长（需满足 CFL）
MoM 大规模加速              → FastSolver: "ACA"（ACA+GMRES，O(N·r) 矩阵向量积）
MoM 介质目标                → Equation: "PMCHWT"，配合 Domains.Materials 指定 ε_r/μ_r
```

---

## 13. WASM / 浏览器限制

| 约束 | 限制 | 说明 |
|------|------|------|
| 线程 | 单线程（无 rayon） | MoM/SBR+ 退化为串行 |
| 内存 | ~30 MB 堆 | MoM 建议 N < 1000 面元 |
| 文件系统 | 无磁盘 IO | 输出返回 Blob URL |
| `rem-mom` | 可用 | rayon 条件编译排除；建议 N < 1000，使用 `FastSolver: "Direct"` |
| `rem-sbr` | 可用 | rayon cfg-excluded |
| `rem-bem` | 可用 | nalgebra LU 支持 WASM |

---

## 14. 技术债务与已知限制

| 项目 | 说明 | 优先级 |
|------|------|--------|
| 3-D 静磁（Nedelec H(curl)） | ✅ 已完成标量解耦矢量位 A=(Ax,Ay,Az) P1 Tet4（3 个 PCG 求解，B=∇×A 恢复）；Nedelec H(curl) 高阶离散待实现 | — |
| 波导端口 TE/TM 场匹配 | ✅ 已完成 v1.0：TE/TM 两类阻抗（Z_TE/Z_TM），模式 n 选取，截止以下退化为 TEM | — |
| AMR 集成 | ✅ 已完成：ZZ 估计器 + Dörfler 标记 + Tri3 红细分 + P1 延拓；**双重收敛判据**（相对 + 绝对 Hz） | — |
| 时域瞬态（TD-FEM） | ✅ 已完成 v1.0：GeneralizedAlpha + IMEX-ARK3 + RK4 + **激励波形 CSV**；Nedelec H(curl) + 完整 Maxwell 矢量场待实现（v2.0） | 中 |
| MoM 介质目标（PMCHWT） | ✅ 已完成：2N×2N PMCHWT 块方程，J+M 未知量，`Equation: "PMCHWT"` 路径已接通 | — |
| MoM ACA 加速 | ✅ 已完成：部分主元 ACA（复对称 Z≈U·V^T），`FastSolver: "ACA"` 已接通；FMM 尚未实现 | — |
| SBR+ 边缘绕射（PTD） | ✅ 已完成：UTD 绕射系数 + 边缘线积分，`write_rcs_with_ptd` 集成进主循环 | — |
| Lanczos 再正交化 | ✅ 已完成：双遍 M-正交化，消除浮点积累伪特征值 | — |
| SSOR 预条件 | ✅ 已完成：ω=1.5 SSOR 替代 Jacobi，FEM 刚度矩阵迭代次数减少 3–5× | — |
| 电容提取（C = 2U/V²） | ✅ 已完成：能量法，WASM UI 显示 pF | — |
| Q 因子导体损耗 | ✅ 已完成：R_s = √(ωμ₀/2σ) 微扰法；1/Q_c 表面积分；与介质 Q_d 合并 | — |
| 电流偶极子源（Palace v0.16） | ✅ 已完成：`Domains.CurrentDipole`；Hertz 偶极子 RHS 注入最近自由节点 | — |
| ROM 电路综合（Palace v0.16） | ✅ 已完成：Vector Fitting（Gustavsen-Semlyen）极点-留数拟合；`DrivenSolver.CircuitSynthesis: bool`；输出 Touchstone .s1p、极点-留数 CSV、SPICE .cir | — |
| p-FEM（P2+）实际应用 | order > 1 时打印警告并降级 P1；P2 装配需接入 fem-rs `H1Space` + `Assembler` API | 低 |
| FMM 加速 | `FastSolver: "FMM"` 配置可识别，运行时返回错误；需实现快速多极子 | 低 |
