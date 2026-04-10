# REM 电磁仿真能力文档

> 版本：v1.0，2026-04-10  
> REM 版本：**v0.17.0**（1496 个测试全部通过，~23,800 行代码，19 个 crate）  
> 对标基准：Palace v0.16（2026-03-09）、Sonnet Suite 19

---

## 总览

REM（Rust Electromagnetic）是一款对标 [Palace](https://github.com/awslabs/palace) 的全波电磁仿真工具，
采用纯 Rust 实现，可编译至 `wasm32-unknown-unknown` 在浏览器中运行。
v0.17.0 覆盖 Palace 全部主要求解器，并额外提供 Palace 不具备的矩量法（MoM）、边界元法（BEM）和 SBR+ 高频求解器，同时提供与商业工具 Sonnet Suite 19 的平面电路分析对齐路线图。

```
测试数：1496（cargo test --workspace），零失败
代码量：~23,800 行（19 个 crate，不含 vendor/）
```

---

## 一、与 Palace v0.16 功能对比

### 1.1 求解器能力总矩阵

| 功能 | Palace v0.16 | REM v0.17 | 说明 |
|------|:---:|:---:|------|
| **静电场** (Electrostatic) | ✅ | ✅ | P1 FEM，变介电常数；C = 2U/V² 电容提取（REM UI 显示 pF）|
| **静磁场** (Magnetostatic) | ✅ | ✅ | P1 FEM（2-D A_z + **3-D A=(Ax,Ay,Az)**），变磁导率，磁能量提取 |
| **特征模** (Eigenmode) | ✅ | ✅ | Lanczos 移位逆迭代（**完全再正交化**），多模式，VTK 模态输出；**AMR 双重收敛判据** |
| **频域驱动** (Driven) | ✅ | ✅ | 频率扫描，S 参数提取，集总端口；**峰值处 E 场恢复**输出 |
| **时域瞬态** (Transient) | ✅ | ✅ | GeneralizedAlpha（2阶无条件稳定）、IMEX-ARK3(2)4L[2]SA（自适应，3阶）、RK4（显式）；**激励波形 CSV 导出** |
| **S 参数提取** | ✅ | ✅ | `postpro/port-S.csv`，Palace 格式兼容 |
| **集总端口** (Lumped Port) | ✅ | ✅ | LumpedPort 激励 + 阻抗边界；多元素端口 (`Elements`) |
| **波导端口** (Wave Port) | ✅ | ✅ | **TE/TM** 1-D 截面特征值场匹配：k_c 计算，Z_TE = ωμ₀/k_z，Z_TM = k_z/(ωε₀)；第 n 阶模式选取 |
| **自适应网格细化** (AMR) | ✅ | ✅ | ZZ 误差估计 + Dörfler 标记 + Tri3 红细分 + P1 延拓；静电/静磁/特征模/驱动均已集成 |
| **高阶基函数** (p-FEM) | ✅ | ✅ | `Solver.Order` 已解析；order > 1 警告并降级 P1；P2+ 装配待完成 |
| **电流偶极子激励** | ✅（v0.16 新增） | ✅ | `Domains.CurrentDipole`；REM Hertz 偶极子 RHS 注入最近自由节点 |
| **ROM 电路综合** | ✅（v0.16 新增） | ✅ | Vector Fitting 极点-留数拟合；`DrivenSolver.CircuitSynthesis: bool`；输出 Touchstone .s1p + SPICE .cir |
| **运行时 JSON Schema 验证** | ✅（v0.16 新增） | ✅ | 两阶段：结构预校验 + 语义校验；5 个单元测试 |
| **内存峰值报告** | ✅（v0.16 新增） | ✅ | Linux VmPeak / Windows PSAPI；`log::info!` 输出 MiB/GiB |
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
| **MoM 集总端口 + S 参数** | ❌ | ✅ | `MomLumpedPort`；face_attrs 端口标签；N×N S 矩阵扫频；Touchstone `.sNp` + Palace CSV 输出 |
| **边界元法 BEM（Laplace P0）** | ❌ | ✅ | Laplace 外 Dirichlet 问题，电容提取 |
| **SBR+ 高频射线追踪 + PO + PTD** | ❌ | ✅ | AABB BVH，两阶段 PO，PTD 边缘绕射修正；ka=10.5 误差 < 0.1 dB |
| **RCS / 远场后处理** | ❌ | ✅ | PO 远场积分，rcs_sbr.csv，多方向扫描 |
| **SSOR 预条件 PCG** | ❌ | ✅ | ω=1.5 SSOR 替代 Jacobi；FEM 刚度矩阵迭代次数减少 3–5× |
| **材料各向异性 ε/μ 张量装配** | ❌ | ✅ | `MaterialAxes` 旋转矩阵 → 3×3 张量；静电/特征模/驱动均已接通 |
| **Drude-Lorentz 频变材料** | ✅ | ✅ | 每频点复数 ε(ω) 修正；与静态损耗互不重叠 |
| **Floquet 周期边界条件（Γ 点）** | ✅ | ✅ | TripletMatrix 节点重映射；几何平移匹配；仅实数（非零波矢待复数支持）|
| **近远场变换（Kirchhoff 积分）** | ❌ | ✅ | Driven 求解器后处理；`far_field.csv`；WASM UI 可导出 |
| **快照 ROM 频率扫描加速** | ✅（自适应 ROM） | ✅ | `DrivenSolver.RomOrder`；正交归一化快照基；r×r 缩减系统；单端口可用 |
| **导体 Q 因子（R_s 微扰法）** | ✅ | ✅ | 表面电阻 R_s = √(ωμ₀/2σ)；微扰积分 1/Q_c；与介质 Q_d 合并 |
| **Touchstone I/O（独立 crate）** | ❌ | ✅ | `rem-touchstone`：N 端口 RI/MA/DB 格式，`write_snp()`，v0.17.0 新增 |

**图例**：✅ 已实现并通过验证　❌ 不支持

### 1.2 与 Palace v0.16 对比总结

- **REM 领先项**：WASM/浏览器运行、MoM（全波散射 CFIE/PMCHWT/ACA）、BEM（Laplace）、SBR+PTD（高频）、RCS 后处理、SSOR 预条件、纯 Rust 零 C++ 依赖、独立 Touchstone crate
- **Palace v0.16 领先项**：HPC 级 AMG 预条件（AMS/ADS）、Nedelec H(curl) 矢量 FEM（完整 Maxwell）
- **等价项**：静电/静磁/特征模/驱动/瞬态 FEM 核心求解器、AMR、S 参数、波导端口、VTK 输出、Palace JSON/YAML 配置兼容

---

## 二、各求解器详细说明

### 2.1 静电场求解器

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
- **Web UI**：`Capacitance: X.X pF`（能量法 C = 2U/V²）

**验证**：平行板电容与解析解 ε₀A/d 误差 < 1e-12

---

### 2.2 静磁场求解器

**问题类型**：`Problem.Type = "Magnetostatic"`

**方法**：P1 FEM，变磁导率 ν(x) = 1/(μ₀μᵣ)；根据 `mesh.dim` 自动选择 2-D 或 3-D 模式。

**2-D 模式（mesh.dim == 2）**：标量磁矢位 A_z

```
−∇·(ν ∇A_z) = J_z      A_z = 0 (Ground)，A_z = 1 (SurfaceCurrent)
B_x =  ∂A_z/∂y,   B_y = −∂A_z/∂x
```

**3-D 模式（mesh.dim == 3）**：三分量解耦矢量位 A = (Ax, Ay, Az)

```
−∇·(ν ∇Aᵢ) = 0   (i = x, y, z)
同一刚度矩阵 K 装配一次，三个 PCG 求解共用
B = ∇×A：Bx = ∂Az/∂y − ∂Ay/∂z,  By = ∂Ax/∂z − ∂Az/∂x,  Bz = ∂Ay/∂x − ∂Ax/∂y
```

**验证**：

| 测试 | 精度 |
|------|------|
| 2-D 线性 A_z = y，铁磁 μ_r=1000 | 误差 < 1e-12 |
| 3-D 线性 Az = z（Tet4 网格） | 误差 < 1e-10 |
| 3-D curl A = (0,x,0) → Bz = 1 | 误差 < 1e-10 |

---

### 2.3 特征模求解器

**问题类型**：`Problem.Type = "Eigenmode"`

**方法**：Lanczos 迭代 + 移位逆（shift-invert），求解广义特征值问题 Kx = λMx；**完全再正交化**（双遍 M-正交化）消除浮点积累引起的伪特征值

**关键配置**（`Solver.Eigenmode`）：

| 字段 | 说明 | 默认 |
|------|------|------|
| `N` | 求解模式数 | 1 |
| `Target` | 目标频率 [Hz] | 0.0 |
| `Tol` | 迭代容差 | 1e-6 |
| `Save` | 保存前 N 个模态 | 1 |

**Q 因子**：含介质 tan δ 微扰损耗（Q_d）及导体表面 R_s 欧姆损耗（Q_c）；1/Q_total = 1/Q_d + 1/Q_c

**AMR 收敛**：相对频率变化 < 1e-4 **或**绝对变化 < 1 MHz 时停止（双重判据）

---

### 2.4 频域驱动求解器

**问题类型**：`Problem.Type = "Driven"`

**方法**：频率域 FEM（Helmholtz 方程），频率扫描；**峰值频率处 φ 场保留**，供 E 场恢复

**波导端口模式**：

| 模式 | 阻抗公式 | 说明 |
|------|---------|------|
| TE | Z_TE = ωμ₀/k_z | 默认，截止以下退化为 TEM |
| TM | Z_TM = k_z/(ωε₀) | `ModeType::Tm` 路径 |

**关键配置**（`Solver.Driven`）：

| 字段 | 说明 |
|------|------|
| `MinFreq` / `MaxFreq` / `FreqStep` | 频率范围和步进 [GHz] |
| `RomOrder` | 快照 ROM 展开点数（0=禁用；建议 4–16；仅单端口+无 Drude-Lorentz 可用）|
| `CircuitSynthesis` | VF 电路综合开关（true 时输出 Touchstone .s1p + SPICE .cir）|

**输出**：
- `postpro/port-S.csv`：S 参数（f, Re(S11), Im(S11), |S11| dB）
- `driven_NNNN.vtk`：各频率步场量
- `far_field.csv`：近远场变换辐射方向图（需 `Solver.FarField` 配置）

---

### 2.5 时域瞬态求解器

**问题类型**：`Problem.Type = "Transient"`

**三种时间积分方案**：

| 方案 | 阶数 | 稳定性 | 适用场景 |
|------|------|--------|---------|
| GeneralizedAlpha | 2阶 | 无条件稳定 | 一般时域仿真（Palace 默认方案） |
| IMEX-ARK3(2)4L[2]SA | 3阶 | 自适应步长 | Kennedy & Carpenter 2003 |
| Explicit RK4 | 4阶 | 显式（需满足 CFL）| 高精度固定步长 |

---

### 2.6 矩量法求解器（MoM）

> **REM 独有，Palace 不支持**

**问题类型**：`Problem.Type = "MoM"`

**方法**：RWG 矢量基函数 + CFIE（PEC）或 PMCHWT（介质），可选 ACA 矩阵压缩加速

**核心模块**（`crates/mom/src/`）：

| 模块 | 说明 |
|------|------|
| `surface_mesh.rs` | 从 RemMesh 提取 PEC 表面三角网格 + 共享边拓扑；`face_attrs` 物理组标签 |
| `quadrature.rs` | Dunavant 三角形高斯求积（阶次 1/3/5/7/9） |
| `green.rs` | 3D Helmholtz Green 函数及法向导数 |
| `singular.rs` | Duffy 自积分 + Sauter-Schwab 奇异积分 |
| `assemble.rs` | EFIE/MFIE/CFIE Z 矩阵装配；内置 LU、GMRES、ACA+GMRES 求解器 |
| `excitation.rs` | 平面波激励向量（θ/φ 极化，任意入射方向） |
| `postprocess.rs` | RCS 远场积分，VTK 表面电流输出 |
| `mie.rs` | Mie 级数解析解（验证用） |
| `basis/rwg.rs` | RWG 基函数评估 + 散度计算 |
| `aca.rs` | **部分主元 ACA**：Z ≈ U·V^T（复对称矩阵）；O(N·r) 矩阵向量积 |
| `pmchwt.rs` | **PMCHWT 介质目标**：2N×2N 块矩阵（T+K 算符，J+M 未知量） |
| `port.rs` | **集总端口（v0.17）**：`MomLumpedPort`，按 face_attrs 标签过滤 RWG；激励 RHS；电流提取 |
| `sparams.rs` | **S 参数扫频（v0.17）**：N×N S 矩阵，LU 扫频；Touchstone `.sNp` + `port-S.csv` 输出 |

**PMCHWT 系统（2N×2N）**：

```
┌ T₁+T₂        K₁+K₂       ┐ ┌ a ┐   ┌ −⟨f, E_inc⟩ ┐
│                            │ │   │ = │             │
└ −(K₁+K₂)   T₁/η₁²+T₂/η₂² ┘ └ b ┘   └ −⟨f, H_inc⟩ ┘
```

**FastSolver 选型建议**：

| N（RWG 基函数数）| 推荐求解器 | 说明 |
|----------------|-----------|------|
| N < 500 | `"Direct"` | Dense LU，O(N³)，精度最高 |
| 500 ≤ N < 3000 | `"GMRES"` | restart=30, tol=1e-8 |
| N ≥ 3000 | **`"ACA"`** | ACA+GMRES；O(N·r) 矩阵向量积，r ≈ log N |

**验证**：PEC 球体（r=0.5 m）@ 1 GHz，kα≈10.5，单站 RCS vs Mie 误差 < 0.5 dB

---

### 2.7 边界元法求解器（BEM）

> **REM 独有，Palace 不支持**

**问题类型**：`Problem.Type = "BEM"`

**方法**：Laplace P0 边界积分方程，外 Dirichlet 问题（PEC 静电）

**BIE 公式**：
```
½φ(r) + ∫_S ∂G_L/∂n'(r,r') φ(r') dS' = ∫_S G_L(r,r') q(r') dS'
G_L(r,r') = 1/(4π|r-r'|)
```

---

### 2.8 SBR+ 高频射线追踪求解器

> **REM 独有，Palace 不支持**

**问题类型**：`Problem.Type = "SBR"`

**方法**：SBR+（Shooting and Bouncing Rays Plus），几何光学 + 物理光学（PO），适用于电大目标（kα >> 1）

**算法（两阶段 PO）**：

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
| `fresnel.rs` | Fresnel 系数，PEC 镜面反射，PO 感应电流 J = 2n̂×H |
| `ptd.rs` | **PTD 边缘绕射修正**：UTD 绕射系数 + 边缘线积分 |
| `output.rs` | `rcs_sbr.csv` + 感应电流 VTK；`write_rcs_with_ptd` 集成 PTD 贡献 |

**验证**：PEC 球（r=0.5 m）@ 1 GHz，ka=10.5，单站 RCS 误差 **0.05 dB**；PTD 修正改善大角度散射精度约 1–2 dB

---

## 三、REM MoM vs Sonnet Suite 19 对比

### 3.1 技术定位

| 维度 | REM MoM | Sonnet Suite 19 |
|------|---------|-----------------|
| **定位** | 通用全波三维表面积分 MoM，适用于任意封闭/开放 3D 散射体 | 平面化三维（2.5D）MoM，专用于多层基板上的平面导体结构 |
| **几何适用** | 任意三维封闭曲面（球体、飞机、天线等任意曲面网格） | 平面导体层叠（MMIC、PCB、微带、槽线、天线阵列平面结构） |
| **目标用户** | 雷达散射截面（RCS）计算、3D 目标散射与辐射 | 微波/毫米波电路 S 参数提取、无源器件建模、RFIC/MMIC 版图验证 |
| **参考工具** | 对标 FEKO（表面积分方程）、WIPL-D | 商业"标准"，常作为 ADS、Cadence AWR 的 EM 联合仿真后端 |

### 3.2 核心算法对比

| 特性 | REM MoM | Sonnet Suite 19 |
|------|---------|-----------------|
| **EFIE** | ✅ RWG 矢量基，3D Helmholtz 核 | ✅ RWG 平面截断版 |
| **MFIE** | ✅ 完整实现 | 内部使用，不直接暴露 |
| **CFIE（EFIE+MFIE 混合）** | ✅ α 参数可配，消除内谐振 | Sonnet 对封闭结构默认采用无内谐振公式 |
| **PMCHWT 介质目标** | ✅ 完整 2N×2N 块矩阵，ε_r/μ_r 可配 | ❌ 不支持任意三维介质目标 |
| **Green 函数** | 均匀自由空间 3D Helmholtz Green 函数 | **分层介质 Green 函数**（Sommerfeld 积分，精确建模多层基板） |
| **ACA 矩阵压缩** | ✅ 部分主元 ACA，O(N·r) | ❌ 不支持 ACA |
| **FFT 加速 MoM** | ❌ 未实现（路线图中） | ✅ **核心优势**：FFT 加速矩阵填充，O(N log N) |
| **平面波激励** | ✅ θ/φ 极化，任意入射角 | ❌ 不支持平面波（仅端口激励） |
| **集总端口 + S 参数** | ✅（v0.17 新增） | ✅ 集总端口，去嵌入，完整 S/Y/Z 矩阵 |
| **RCS 输出** | ✅ θ/φ 扫描，全球面 dBsm | ❌ 不支持 RCS |

### 3.3 材料与物理建模对比

| 特性 | REM MoM | Sonnet Suite 19 |
|------|---------|-----------------|
| **PEC 导体** | ✅ CFIE 完整实现 | ✅ 理想 PEC |
| **介质目标（均匀）** | ✅ PMCHWT（2N×2N） | ❌ 不支持任意三维介质散射 |
| **分层介质基板** | ❌ 仅自由空间（路线图 v0.18） | ✅ **核心优势**：Sommerfeld 积分精确建模 |
| **有损导体（表面阻抗）** | ❌（路线图 v0.19） | ✅ 有限电导率 σ，R_s 建模 |
| **各向异性介质** | ❌ MoM 无；FEM 求解器支持 3×3 张量 | ❌ 各向同性基板 |

### 3.4 平台与工程集成对比

| 特性 | REM MoM | Sonnet Suite 19 |
|------|---------|-----------------|
| **编程语言** | 纯 Rust（零 C++ 依赖，内存安全） | C++（专有），闭源 |
| **授权** | 开源（MIT/Apache 协议） | 商业授权（年费/节点锁定） |
| **WASM 支持** | ✅ `wasm32-unknown-unknown`，浏览器运行 | ❌ 仅本地 Windows/Linux |
| **EDA 工具集成** | ❌ 独立工具 | ✅ ADS、AWR、Cadence 无缝集成 |
| **参数化扫描 / 优化** | ❌ | ✅ 内置参数扫描 + 梯度优化 |
| **版本控制友好** | ✅ 纯文本 JSON/YAML 配置 | 部分版本可比较的 `.son` 格式 |
| **工程成熟度** | v0.17.0 早期版本 | 40 年商业产品历史 |

### 3.5 场景选型建议

```
┌─────────────────────────────────────────────────────────────────┐
│  仿真场景                          推荐工具                     │
├─────────────────────────────────────────────────────────────────┤
│  任意三维目标 RCS（飞机、导弹、球）  → REM MoM（CFIE/PMCHWT）  │
│  三维均匀介质散射体                 → REM MoM（PMCHWT）         │
│  平面波入射、雷达截面方向图         → REM MoM                   │
│  微带/CPW 滤波器/耦合器 S 参数      → Sonnet 19                 │
│  MMIC/RFIC 版图 EM 验证             → Sonnet 19                 │
│  多层 PCB 串扰/匹配分析             → Sonnet 19                 │
│  贴片天线 S11 + 方向图              → Sonnet 19（S11）+ REM FEM │
│  超大规模平面电路（N > 50,000）      → Sonnet 19（FFT MoM）     │
│  无授权/开源/嵌入浏览器仿真         → REM MoM                   │
│  三维介质与 RCS 联合分析            → REM MoM + SBR+ PTD       │
└─────────────────────────────────────────────────────────────────┘
```

---

## 四、MoM → Sonnet 对齐路线图

### 4.1 执行摘要

Sonnet 19 的核心竞争力集中在四个维度，REM 当前差距及优先级：

| 差距维度 | Sonnet 核心能力 | REM 现状 | 优先级 |
|---------|---------------|---------|--------|
| **Green 函数** | Sommerfeld 分层介质积分 | 仅自由空间 | P0（物理根基） |
| **端口 + S 参数** | 集总/波导端口，Touchstone 输出 | **v0.17 已完成集总端口** | P0（已交付） |
| **快速算法** | FFT O(N log N) 加速 | ACA O(N·r)（三维散射） | P1（规模扩展） |
| **导体损耗** | 有限电导率 SIBC | 仅 PEC | P1（材料完整性） |

### 4.2 版本路线图

```
v0.17.0（当前）── MoM 集总端口激励 + N×N S 参数 + Touchstone .sNp  ✅ 已完成
    │
    ▼
v0.18.0 ── 分层介质 Green 函数（Sommerfeld DCIM）
    │       新 crate：crates/layered_green
    │       GreenFunction trait 抽象（FreeSpaceGreen / LayeredGreen）
    │       FR4 单层基板 + 贴片天线谐振验证（误差 < 1%）
    ▼
v0.19.0 ── 有损导体 SIBC + FFT 加速平面 MoM
    │       Z_s = (1+j)/(σδ_s)，SIBC 修正 CFIE 对角块
    │       FFT 加速：O(N log N)（仅平面结构，运行时自动检测）
    ▼
v0.20.0 ── MoM AMR + 快照 ROM + Touchstone 完整兼容
            MoM 内 Dörfler 标记 + Tri3 中线分割
            MoM ROM（复用 Driven FEM ROM 模式）
            Touchstone 2.0 完整兼容（RI/MA/DB，MA/DB 注释行）
```

### 4.3 阶段 17 详情：分层介质 Green 函数（v0.18.0）

**物理模型**：Sommerfeld 积分 Green 函数适用于平面分层介质（基板叠层）：

```
G_A(r, r') = (1/4π) ∫₀^∞ g_A(k_ρ, z, z') J₀(k_ρ ρ) k_ρ dk_ρ
```

其中 g_A 通过传输矩阵法（TMM）计算。

**选用离散复像法（DCIM）**：一次 GPOF 拟合后 O(M_pole) per pair，比直接数值 Sommerfeld 积分快 100–1000×。

**新 Crate**：`crates/layered_green/`（`layer.rs`, `transfer_matrix.rs`, `sommerfeld.rs`, `discrete_image.rs`）

**与 MoM 集成**：抽象 `GreenFunction` trait，`assemble_efie_rwg` 接受 `&dyn GreenFunction`，无 Substrate 配置时回退至 `FreeSpaceGreen`。

**配置扩展示例**：
```json
"Solver": {
  "MoM": {
    "Substrate": {
      "BottomPec": true,
      "Layers": [
        { "Permittivity": 4.4, "LossTan": 0.02, "Thickness": 1.6e-3 },
        { "Permittivity": 1.0, "LossTan": 0.0,  "Thickness": 1.0    }
      ]
    }
  }
}
```

### 4.4 阶段 18 详情：SIBC + FFT 加速（v0.19.0）

**SIBC 物理模型**：

```
E_tan = Z_s · (n̂ × H_tan)
Z_s = (1+j) / (σ · δ_s)，δ_s = √(2/(ωμσ)) 为趋肤深度
```

**FFT 加速适用条件**：所有 RWG 基函数源/观测点均在同一水平层（z≈常数），Green 函数退化为卷积形式；运行时检测平面性，非平面时自动降级 ACA。

### 4.5 阶段 19 详情：AMR + ROM + Touchstone（v0.20.0）

**MoM AMR**：表面电流密度梯度误差指示器（η_m = ||∇J_s||_{T_m} × h_m），Dörfler 标记 + Tri3 中线分割 + RWG 重映射。

**MoM ROM**：锚点频率快照基 V（N×r），在其余频率用低维 r×r 近似（复用 Driven FEM ROM 模式）。

**完整 Touchstone 2.0**：`.s1p`/`.s2p`/`.sNp` 自动选后缀；`GHz`/`MHz`/`Hz` 频率单位；`MA`/`RI`/`DB` 数据格式；标准 `# GHz S MA R 50` 选项行；N×N 矩阵按 Touchstone 2.0 规范排列。

---

## 五、Palace 配置兼容性

REM 完全兼容 Palace JSON/YAML 配置文件格式（含 v0.16 新关键字）。Palace 用户无需修改已有配置即可在 REM 中运行。

### 5.1 支持的边界类型

| Palace 字段 | REM 支持 | 说明 |
|------------|:---:|------|
| `Boundaries.PEC` | ✅ | `Pec` |
| `Boundaries.PMC` | ✅ | Neumann 自然边界 |
| `Boundaries.Impedance` | ✅ | 解析 |
| `Boundaries.Absorbing` | ✅ | 解析 |
| `Boundaries.Conductivity` | ✅ | 解析，Q 因子 R_s 微扰法 |
| `Boundaries.Ground` | ✅ | `Ground` |
| `Boundaries.Terminal` | ✅ | v0.16 要求必填，REM 已支持 |
| `Boundaries.LumpedPort` | ✅ | 含 `Elements` 多元素端口 |
| `Boundaries.WavePort` | ✅ | TE/TM，Z_TE/Z_TM，模式 n 选取 |
| `Boundaries.SurfaceCurrent` | ✅ | v0.16 要求必填，REM 已支持 |

### 5.2 REM 专有扩展（Palace 静默忽略）

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

## 六、求解器选型指南

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
  · GeneralizedAlpha（v1.0）→ time_scheme: "GeneralizedAlpha"，2阶，无条件稳定
  · IMEX-ARK 自适应（v1.0）→ time_scheme: "ARKODE"，3阶，自适应步长
  · Explicit RK4（v1.0）   → time_scheme: "RungeKutta"，4阶，固定步长（需满足 CFL）
平面电路 MoM（嵌入基板）   → 等待 v0.18 分层 Green 函数；当前建议 Sonnet 19
MoM 大规模三维散射加速      → FastSolver: "ACA"（ACA+GMRES，O(N·r)）
MoM 介质目标                → Equation: "PMCHWT"，配合 Domains.Materials 指定 ε_r/μ_r
MoM S 参数提取              → Problem.Type = "MoM"，配 Ports 列表（v0.17+）
```

---

## 七、已验证示例

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

## 八、WASM / 浏览器限制

| 约束 | 限制 | 说明 |
|------|------|------|
| 线程 | 单线程（无 rayon） | MoM/SBR+ 退化为串行 |
| 内存 | ~30 MB 堆 | MoM 建议 N < 1000 面元 |
| 文件系统 | 无磁盘 IO | 输出返回 Blob URL |
| `rem-mom` | 可用 | rayon 条件编译排除；建议 N < 1000，使用 `FastSolver: "Direct"` |
| `rem-sbr` | 可用 | rayon cfg-excluded |
| `rem-bem` | 可用 | nalgebra LU 支持 WASM |

---

## 九、已知限制与技术债务

| 项目 | 状态 | 优先级 |
|------|------|--------|
| Nedelec H(curl) 矢量 FEM | 尚未实现；当前静磁用标量解耦矢量位 A=(Ax,Ay,Az) P1 | 中 |
| 时域完整 Maxwell 矢量场（TD-FEM v2.0） | 待实现；当前为标量 P1 FEM 三方案 | 中 |
| p-FEM（P2+ 实际应用） | order > 1 打印警告并降级 P1；P2 装配需接入 fem-rs H1Space | 低 |
| FMM 加速 | `FastSolver: "FMM"` 配置项已预留，运行时返回错误 | 低 |
| MoM 分层介质 Green 函数 | 路线图 v0.18 | P0 |
| MoM 有损导体 SIBC | 路线图 v0.19 | P1 |
| MoM FFT 加速 | 路线图 v0.19 | P1 |
| MoM AMR | 路线图 v0.20 | 低 |
| Driven solver 复数 PCG | 当前使用实数 PCG；高频复数问题收敛需 VF ROM 路径 | 中 |
| Floquet 非零波矢 | 非零 Floquet 波矢警告并跳过；待复数矩阵支持 | 低 |

---

## 十、版本历史

| 版本 | 亮点 |
|------|------|
| v0.17.0 | `rem-touchstone` 独立 crate；MoM 集总端口 + N×N S 参数 + Touchstone；子模块 fem-rs/rmsh 更新 |
| v0.16.0 | ROM Vector Fitting 电路综合（SPICE .cir）；完整 N×N S 矩阵；近远场变换；快照 ROM |
| v0.15.0 | 各向异性 ε/μ 张量；导体 Q 因子；电流偶极子激励；Floquet 周期边界；Drude-Lorentz 频变材料；JSON Schema 校验 |
| v0.14.0 | 时域瞬态（GeneralizedAlpha + IMEX-ARK3 + RK4）；激励波形 CSV 导出 |
| v0.13.0 | 3-D 静磁矢量位 A=(Ax,Ay,Az)；MoM PMCHWT 介质目标；ACA 矩阵压缩；SBR+ PTD 边缘绕射 |
| v0.12.0 | WavePort TE/TM 场匹配；AMR（ZZ+Dörfler+红细分） |

---

## 参考资料

- [Palace（AWS）](https://github.com/awslabs/palace) — REM 对标的开源 EM 仿真工具
- [Sonnet Software](https://www.sonnetsoftware.com) — 平面电路 2.5D MoM 商业工具
- Rao, Wilton, Glisson, "Electromagnetic Scattering by Surfaces of Arbitrary Shape," IEEE TAP, 1982
- Harrington, *Field Computation by Moment Methods*, IEEE Press, 1993
- Fang, *Analytical and Numerical Methods in Electromagnetic Wave Theory*, 1993（DCIM 基础）
- Chow, Mosig et al., "Discrete complex image representation," IEEE TAP, 1994
- Gustavsen & Semlyen, "Rational approximation of frequency domain responses," IEEE TPWRD, 1999
- Touchstone File Format Specification, Version 2.0, IPC-2141, 2009
