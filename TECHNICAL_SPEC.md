# REM — Rust Electromagnetic Solver
## Technical Specification v0.4

> **目标**: 用纯 Rust（可编译至 `wasm32-unknown-unknown`）实现对标 Palace 的全波电磁仿真工具，
> 基于 [fem-rs](https://github.com/javagg/fem-rs) 通用有限元库，兼容 Palace JSON/YAML 配置格式。

---

## 1. 项目范围与对标目标

### 1.1 Palace 功能覆盖矩阵

| 功能领域 | Palace | REM v0.1 | REM v0.2 | REM v1.0 (当前) |
|----------|--------|----------|----------|----------|
| 静电场 (Electrostatic) | ✅ | ✅ | ✅ | ✅ |
| 静磁场 (Magnetostatic) | ✅ | ✅ | ✅ | ✅ |
| 特征模 (Eigenmode) | ✅ | 🔲 | ✅ | ✅ |
| 频域驱动 (Driven) | ✅ | 🔲 | ✅ | ✅ |
| 时域瞬态 (Transient) | ✅ | 🔲 | 🔲 | ✅ |
| S 参数提取 | ✅ | 🔲 | ✅ | ✅ |
| 集总端口 (Lumped Port) | ✅ | 🔲 | ✅ | ✅ |
| 波导端口 (Wave Port) | ✅ | 🔲 | 🔲 | ✅ |
| 自适应网格细化 (AMR) | ✅ | 🔲 | 🔲 | ✅ |
| 高阶基函数 (p-FEM) | ✅ | ✅ | ✅ | ✅ |
| GMSH 网格导入 | ✅ | ✅ | ✅ | ✅ |
| ParaView/VTK 输出 | ✅ | ✅ | ✅ | ✅ |
| JSON 配置文件 | ✅ | ✅ | ✅ | ✅ |
| YAML 配置文件 | ✅ | ✅ | ✅ | ✅ |
| WASM 目标 | ❌ | ✅ | ✅ | ✅ |
| MPI 并行（native rsmpi） | ✅ | 🔲 | ✅ | ✅ |
| MPI 模拟（jsmpi + Web Worker） | ❌ | ✅ | ✅ | ✅ |
| **矩量法 MoM (EFIE/CFIE/PMCHWT+ACA)** | ❌ | ❌ | ❌ | ✅ **已实现** |
| **边界元法 BEM (Laplace P0)** | ❌ | ❌ | ❌ | ✅ **已实现** |
| **SBR+ 高频射线追踪 + PO + PTD** | ❌ | ❌ | ❌ | ✅ **已实现** |
| RCS / 远场后处理 | ❌ | ❌ | ❌ | ✅ **已实现** |

### 1.2 核心差异化特性

- **纯 Rust + WASM**: 无 C/C++ 依赖，可在浏览器运行
- **Palace 配置兼容**: 直接读取 Palace JSON/YAML 配置文件（详见第 7 节）
- **fem-rs 驱动**: 复用已验证的 FEM 基础设施
- **MoM/BEM/SBR+ 扩展**: 矩量法（`crates/mom`，支持 CFIE/PMCHWT/ACA）、Laplace P0 BEM（`crates/bem`）以及 SBR+PTD 高频射线追踪求解器（`crates/sbr`）均已实现，与 FEM 共享网格/配置层
- **超集定位**: Palace 完全兼容 + MoM/BEM/SBR+/PTD/RCS 额外能力，不破坏 Palace 用户现有工作流

---

## 2. 配置文件格式规范（Palace 兼容）

### 2.1 顶层结构

```json
{
  "Problem": { ... },
  "Model": { ... },
  "Domains": { ... },
  "Boundaries": { ... },
  "Solver": { ... }
}
```

YAML 等价：
```yaml
Problem:
  Type: Electrostatic
  Verbose: 1
  Output: ./output
```

### 2.2 Problem 节

```json
{
  "Problem": {
    "Type": "Electrostatic",
    "Verbose": 1,
    "Output": "./output"
  }
}
```

| 字段 | 类型 | 默认 | 说明 |
|------|------|------|------|
| `Type` | string | 必填 | `"Electrostatic"` \| `"Magnetostatic"` \| `"Eigenmode"` \| `"Driven"` \| `"Transient"` |
| `Verbose` | int | 1 | 日志详细度 0-3 |
| `Output` | string | `"."` | 结果输出目录 |

### 2.3 Model 节

```json
{
  "Model": {
    "Mesh": "path/to/mesh.msh",
    "L0": 1.0e-3,
    "Refinement": {
      "MaxIter": 0,
      "Tol": 1.0e-2,
      "Nonconformal": false
    }
  }
}
```

| 字段 | 类型 | 默认 | 说明 |
|------|------|------|------|
| `Mesh` | string | 必填 | 网格文件路径（GMSH .msh） |
| `L0` | float | 1.0 | 全局长度单位（米） |
| `Refinement.MaxIter` | int | 0 | AMR 最大迭代次数，0 禁用 |
| `Refinement.Tol` | float | 1e-2 | AMR 误差容限 |
| `Refinement.Nonconformal` | bool | false | 非协调细化 |

### 2.4 Domains 节（材料定义）

```json
{
  "Domains": {
    "Materials": [
      {
        "Attributes": [1, 2],
        "Permeability": 1.0,
        "Permittivity": 4.5,
        "LossTan": 0.02,
        "Conductivity": 0.0
      },
      {
        "Attributes": [3],
        "Permeability": 1000.0,
        "Permittivity": 1.0
      }
    ]
  }
}
```

| 字段 | 类型 | 默认 | 说明 |
|------|------|------|------|
| `Attributes` | int[] | 必填 | GMSH 物理组 ID 列表 |
| `Permittivity` | float | 1.0 | 相对介电常数 εᵣ |
| `Permeability` | float | 1.0 | 相对磁导率 μᵣ |
| `LossTan` | float | 0.0 | 介质损耗角正切 tan δ |
| `Conductivity` | float | 0.0 | 电导率 σ [S/m] |

### 2.5 Boundaries 节

```json
{
  "Boundaries": {
    "PEC": {
      "Attributes": [1, 2, 3]
    },
    "PMC": {
      "Attributes": [4]
    },
    "Impedance": [
      {
        "Attributes": [5],
        "Rs": 377.0,
        "Ls": 0.0,
        "Cs": 0.0
      }
    ],
    "LumpedPort": [
      {
        "Index": 1,
        "Attributes": [10],
        "Direction": "+X",
        "R": 50.0,
        "Excitation": true
      }
    ],
    "Ground": {
      "Attributes": [20]
    },
    "ZeroCharge": {
      "Attributes": [21]
    },
    "Absorbing": {
      "Attributes": [30],
      "Order": 1
    }
  }
}
```

**边界条件类型映射**：

| Palace 边界 | 数学表达 | 适用问题 |
|-------------|----------|---------|
| `PEC` | **n** × **E** = 0 | 频域、时域、静磁 |
| `PMC` | **n** × **H** = 0 | 频域、时域 |
| `Impedance` | **n** × **H** = Y(**n** × **E** × **n**) | 频域、时域 |
| `LumpedPort` | 集总激励端口 | 频域驱动 |
| `WavePort` | 波导模式激励 | 频域驱动 |
| `Ground` | φ = 0 | 静电、静磁 |
| `ZeroCharge` | ∂φ/∂n = 0 | 静电 |
| `Absorbing` | 一阶/二阶吸收 BC | 频域、时域 |

### 2.6 Solver 节

```json
{
  "Solver": {
    "Order": 1,
    "Eigenmode": {
      "N": 10,
      "Tol": 1.0e-6,
      "MaxIter": 200,
      "Target": 5.0e9,
      "Save": 2
    },
    "Driven": {
      "MinFreq": 1.0e9,
      "MaxFreq": 10.0e9,
      "FreqStep": 0.1e9,
      "SaveStep": 1,
      "AdaptiveTol": 1.0e-2
    },
    "Transient": {
      "Type": "GeneralizedAlpha",
      "MaxTime": 10.0e-9,
      "TimeStep": 1.0e-11,
      "SaveStep": 100
    },
    "Linear": {
      "Type": "GMRES",
      "Tol": 1.0e-6,
      "MaxIter": 200,
      "KSPMGCycleIter": 1,
      "MGLevels": 10,
      "MGCoarsenType": "Logarithmic",
      "PCType": "JACOBI"
    }
  }
}
```

---

## 3. 物理方程与数值方法

### 3.1 静电场

**强形式**:
```
−∇·(ε₀εᵣ ∇φ) = ρ    in Ω
φ = φ_D               on Γ_D  (PEC/Ground)
ε₀εᵣ ∂φ/∂n = 0       on Γ_N  (ZeroCharge/PMC)
```

**弱形式（Galerkin）**:
```
∫_Ω ε₀εᵣ ∇φ·∇v dΩ = ∫_Ω ρv dΩ + ∫_Γ_N g_N v ds
```

**后处理**:
- 电场: **E** = −∇φ
- 表面电荷密度: σ = ε₀εᵣ ∂φ/∂n
- 电容矩阵: Cᵢⱼ = −∂Qᵢ/∂Vⱼ

### 3.2 静磁场（A-φ 公式）

**2D 强形式**（A_z 分量）:
```
−∇·(ν ∇A_z) = J_z    in Ω
A_z = 0               on Γ_D  (PEC)
ν ∂A_z/∂n = 0        on Γ_N  (PMC)
```

**3D 强形式**（矢量位 **A**）:
```
∇×(ν ∇×**A**) = **J**    in Ω
**A** = 0                  on Γ_D
ν(∇×**A**) × **n** = 0   on Γ_N
```

**后处理**:
- 磁通密度: **B** = ∇×**A**
- 磁场强度: **H** = ν**B** = (1/μ)**B**
- 电感矩阵: Lᵢⱼ = Φᵢ/Iⱼ

### 3.3 频域全波（Driven/Eigenmode）

**矢量波动方程**（时谐 e^{jωt}）:
```
∇×(μᵣ⁻¹ ∇×**E**) − k₀²(εᵣ − jσ/(ωε₀))**E** = −jωμ₀**J**_s
```

其中 k₀ = ω√(μ₀ε₀)

**离散化**: Nedelec (H(curl)) 棱元，需 fem-rs Phase 5（待实现）

### 3.4 时域（Transient）

**FDTD 等价连续形式**:
```
ε ∂**E**/∂t = ∇×**H** − σ**E**
μ ∂**H**/∂t = −∇×**E**
```

**时间积分**: Generalized-α 法（隐式，无条件稳定）

---

## 4. 模块架构

```
rem2/
├── Cargo.toml                   # workspace
├── crates/
│   ├── core/                    # 公共类型、常量、错误
│   │   └── src/
│   │       ├── lib.rs
│   │       ├── constants.rs     # ε₀, μ₀, η₀
│   │       └── error.rs
│   │
│   ├── config/                  # Palace 配置解析
│   │   └── src/
│   │       ├── lib.rs
│   │       ├── schema.rs        # 完整 Palace JSON schema 映射
│   │       ├── json_parser.rs   # serde_json 解析（含注释剥除）
│   │       └── yaml_parser.rs   # serde_yaml 解析
│   │
│   ├── mesh/                    # 网格适配层（包装 fem-io）
│   │   └── src/
│   │       ├── lib.rs
│   │       ├── gmsh.rs          # GMSH .msh 读取（复用 fem-io）
│   │       └── mesh_data.rs     # RemMesh 结构（含物理组映射）
│   │
│   ├── materials/               # 材料模型
│   │   └── src/
│   │       ├── lib.rs
│   │       ├── material.rs      # 材料属性结构
│   │       └── domain_map.rs    # 物理组 → 材料 映射
│   │
│   ├── bc/                      # 边界条件处理
│   │   └── src/
│   │       ├── lib.rs
│   │       ├── types.rs         # PEC/PMC/Impedance/Port 枚举
│   │       └── applicator.rs    # 将 BC 应用到 fem-space DOF
│   │
│   ├── electrostatic/           # 静电场求解器
│   │   └── src/
│   │       ├── lib.rs
│   │       ├── solver.rs        # 主求解流程
│   │       └── postproc.rs      # E 场、电容矩阵
│   │
│   ├── magnetostatic/           # 静磁场求解器
│   │   └── src/
│   │       ├── lib.rs
│   │       ├── solver.rs
│   │       └── postproc.rs      # B/H 场、电感矩阵
│   │
│   ├── eigenmode/               # 特征模求解器（v0.2）
│   │   └── src/
│   │       ├── lib.rs
│   │       └── arpack_rs.rs     # 固有值求解器接口
│   │
│   ├── driven/                  # 频域驱动求解器（v0.2）
│   │   └── src/
│   │       ├── lib.rs
│   │       ├── solver.rs
│   │       └── sparams.rs       # S 参数计算
│   │
│   ├── result/                  # 结果输出
│   │   └── src/
│   │       ├── lib.rs
│   │       ├── vtk.rs           # VTK/VTU 输出（复用 fem-io）
│   │       └── csv.rs           # Palace 兼容 CSV 输出
│   │
│   ├── cli/                     # 命令行入口
│   │   └── src/
│   │       └── main.rs
│   │
│   ├── parallel/                # 并行通信抽象层
│   │   └── src/
│   │       ├── lib.rs           # Comm trait + build_comm()
│   │       ├── serial.rs        # SerialComm（单进程 fallback）
│   │       ├── mpi_comm.rs      # MpiComm（feature = "mpi"，rsmpi）
│   │       ├── worker_comm.rs   # WorkerComm（target = wasm32，Web Worker 模拟）
│   │       └── partition.rs     # 网格分区（行分区 / METIS 接口）
│   │
│   ├── mom/                     # 矩量法求解器 ✅ 已实现
│   │   └── src/
│   │       ├── lib.rs           # run() / run_with_mesh() 入口
│   │       ├── surface_mesh.rs  # RWG 面网格提取 + 边拓扑
│   │       ├── quadrature.rs    # Dunavant 三角形高斯求积规则（阶次 1/3/5/7/9）
│   │       ├── green.rs         # 3D Helmholtz Green 函数
│   │       ├── singular.rs      # Duffy 自积分 + Sauter-Schwab 奇异积分
│   │       ├── assemble.rs      # EFIE/CFIE Z 矩阵装配（nalgebra dense + rayon）
│   │       ├── excitation.rs    # 平面波激励向量
│   │       ├── postprocess.rs   # RCS CSV + VTK 输出
│   │       ├── mie.rs           # Mie 级数解析解（验证用）
│   │       └── basis/
│   │           ├── mod.rs
│   │           └── rwg.rs       # RWG 矢量基函数（内部边拓扑）
│   │
│   ├── bem/                     # Laplace P0 BEM ✅ 已实现
│   │   └── src/
│   │       ├── lib.rs           # run() 入口
│   │       ├── kernel.rs        # Laplace Green 函数及法向导数
│   │       ├── assemble.rs      # V/K 矩阵装配（P0 基函数 + Duffy 对角）
│   │       ├── solve.rs         # nalgebra LU 求解
│   │       └── postprocess.rs   # 电容矩阵 + 电位 VTK 输出
│   │
│   ├── sbr/                     # SBR+ 高频射线追踪 + PO ✅ 已实现
│   │   ├── src/
│   │   │   ├── lib.rs           # run() / run_with_mesh() 入口，两阶段 PO 算法
│   │   │   ├── bvh.rs           # AABB BVH（SAH 分割，Möller-Trumbore 求交）
│   │   │   ├── ray.rs           # Ray / RayHit 数据结构
│   │   │   ├── fresnel.rs       # Fresnel 系数 + PEC 镜面反射 + PO 感应电流
│   │   │   ├── excitation.rs    # 平面波激励 + 孔径射线发射
│   │   │   ├── po_integral.rs   # 远场 PO 积分 → RCS
│   │   │   └── output.rs        # RCS CSV + 感应电流 VTK
│   │   ├── src/bin/
│   │   │   └── gen_sbr_meshes.rs # 球面/平板测试网格生成工具
│   │   └── tests/
│   │       └── mie_validation.rs # 集成测试：SBR+ vs Mie（ka≈10.5，误差 < 0.1 dB）
│   │
│   └── wasm/                    # WASM 绑定（支持 SBR+）
│       └── src/
│           ├── lib.rs
│           └── api.rs           # wasm-bindgen JS API
│
├── crates/yew-app/              # Yew 前端（纯 Rust WASM，trunk 构建）
│   ├── Cargo.toml
│   ├── Trunk.toml
│   ├── index.html
│   └── src/
│       ├── main.rs              # Yew App 组件
│       ├── examples.rs          # 8 个 Palace 示例配置 + 网格数据
│       ├── solver.rs            # 直接调用 rem-* 求解器
│       └── style.css            # 样式
│
├── examples/
│   ├── parallel_plate.json      # Palace 格式示例配置
│   ├── coaxial.json
│   ├── resonator.json
│   └── meshes/
│
└── tests/
    ├── integration/
    │   ├── test_electrostatic.rs
    │   ├── test_magnetostatic.rs
    │   └── test_config_compat.rs
    └── fixtures/
        └── palace_configs/      # Palace 官方示例配置（用于兼容性测试）
```

---

## 5. 关键 Rust 接口设计

### 5.1 配置解析

```rust
// config/src/schema.rs

#[derive(Debug, Clone, Deserialize)]
pub struct PalaceConfig {
    #[serde(rename = "Problem")]
    pub problem: Problem,
    #[serde(rename = "Model")]
    pub model: Model,
    #[serde(rename = "Domains")]
    pub domains: Domains,
    #[serde(rename = "Boundaries")]
    pub boundaries: Boundaries,
    #[serde(rename = "Solver")]
    pub solver: SolverConfig,
}

#[derive(Debug, Clone, Deserialize)]
pub struct Problem {
    #[serde(rename = "Type")]
    pub problem_type: ProblemType,
    #[serde(rename = "Verbose", default = "default_verbose")]
    pub verbose: u8,
    #[serde(rename = "Output", default = "default_output")]
    pub output: String,
}

#[derive(Debug, Clone, Deserialize, PartialEq)]
#[serde(rename_all = "PascalCase")]
pub enum ProblemType {
    Electrostatic,
    Magnetostatic,
    Eigenmode,
    Driven,
    Transient,
    // REM extensions (ignored by Palace)
    MoM,
    BEM,
    SBR,
}

// config/src/lib.rs
pub fn load_config(path: &Path) -> Result<PalaceConfig, RemError> {
    let content = std::fs::read_to_string(path)?;
    match path.extension().and_then(|e| e.to_str()) {
        Some("json") => load_json(&content),
        Some("yaml") | Some("yml") => load_yaml(&content),
        _ => Err(RemError::UnknownFormat),
    }
}

fn load_json(s: &str) -> Result<PalaceConfig, RemError> {
    // 剥除 C++ 风格注释再解析
    let stripped = strip_comments(s);
    // 展开属性范围 "1,3-5,6" → [1,3,4,5,6]
    let expanded = expand_attribute_ranges(&stripped);
    serde_json::from_str(&expanded).map_err(Into::into)
}
```

### 5.2 主求解器派发

```rust
// rem-cli/src/main.rs
fn main() -> Result<(), RemError> {
    let args = Args::parse();
    let config = rem_config::load_config(&args.config)?;
    
    match config.problem.problem_type {
        ProblemType::Electrostatic  => rem_electrostatic::run(&config),
        ProblemType::Magnetostatic  => rem_magnetostatic::run(&config),
        ProblemType::Eigenmode      => rem_eigenmode::run(&config),
        ProblemType::Driven         => rem_driven::run(&config),
        ProblemType::Transient      => rem_transient::run(&config),
    }
}
```

### 5.3 静电场求解器核心

```rust
// electrostatic/src/solver.rs
pub fn run(config: &PalaceConfig) -> Result<(), RemError> {
    // 1. 加载网格
    let mesh = fem_io::read_msh_file(&config.model.mesh)?;
    let rem_mesh = RemMesh::from_fem_mesh(mesh, &config.domains.materials)?;
    
    // 2. 建立 FEM 空间
    let space = H1Space::new(&rem_mesh.inner, ElementOrder::P1);
    
    // 3. 组装刚度矩阵（变系数 ε）
    let epsilon_fn = rem_mesh.permittivity_fn();
    let integrator = DiffusionIntegrator::new(move |x, y, z| {
        EPS0 * epsilon_fn(x, y, z)
    });
    let mut K = Assembler::assemble_bilinear(&space, &[&integrator]);
    let mut f = Assembler::assemble_linear(&space, &[/* 体电荷 */]);
    
    // 4. 应用边界条件
    apply_pec_bc(&mut K, &mut f, &config.boundaries.pec, &space);
    apply_ground_bc(&mut K, &mut f, &config.boundaries.ground, &space);
    apply_neumann_bc(&mut f, &config.boundaries.zero_charge, &space);
    
    // 5. 求解
    let solver_cfg = SolverConfig::from(&config.solver.linear);
    let phi = solve_linear(K, f, solver_cfg)?;
    
    // 6. 后处理
    let e_field = gradient_recovery(&phi, &space);
    let capacitance = compute_capacitance_matrix(&phi, &space, &config)?;
    
    // 7. 输出
    write_vtk(&config.problem.output, &phi, &e_field, &rem_mesh)?;
    write_csv_energy(&config.problem.output, &phi, &space)?;
    
    Ok(())
}
```

### 5.4 边界元扩展接口（预留）

```rust
// core/src/bem.rs (v1.0 扩展点)
pub trait BemKernel: Send + Sync {
    /// Green 函数 G(r, r')
    fn green(&self, r: &[f64; 3], r_prime: &[f64; 3]) -> f64;
    /// ∂G/∂n'
    fn green_normal_deriv(&self, r: &[f64; 3], r_prime: &[f64; 3], n_prime: &[f64; 3]) -> f64;
}

pub struct ElectrostaticKernel;
impl BemKernel for ElectrostaticKernel {
    fn green(&self, r: &[f64; 3], r_prime: &[f64; 3]) -> f64 {
        let dist = euclidean_dist(r, r_prime);
        1.0 / (4.0 * PI * EPS0 * dist)
    }
    // ...
}
```

---

## 6. 依赖清单

```toml
# workspace Cargo.toml
[workspace]
resolver = "2"
members = [
    "crates/core",
    "crates/config",
    "crates/mesh",
    "crates/materials",
    "crates/bc",
    "crates/electrostatic",
    "crates/magnetostatic",
    "crates/eigenmode",
    "crates/driven",
    "crates/result",
    "crates/cli",
    "crates/wasm",
]

[workspace.dependencies]
# FEM 基础库
fem-core    = { git = "https://github.com/javagg/fem-rs.git" }
fem-mesh    = { git = "https://github.com/javagg/fem-rs.git" }
fem-element = { git = "https://github.com/javagg/fem-rs.git" }
fem-linalg  = { git = "https://github.com/javagg/fem-rs.git" }
fem-io      = { git = "https://github.com/javagg/fem-rs.git" }
fem-space   = { git = "https://github.com/javagg/fem-rs.git" }
fem-assembly= { git = "https://github.com/javagg/fem-rs.git" }
fem-solver  = { git = "https://github.com/javagg/fem-rs.git" }
fem-amg     = { git = "https://github.com/javagg/fem-rs.git" }

# 配置解析
serde       = { version = "1", features = ["derive"] }
serde_json  = "1"
serde_yaml  = "0.9"

# CLI
clap        = { version = "4", features = ["derive"] }

# 数值
nalgebra    = "0.33"
num-complex = "0.4"

# WASM
wasm-bindgen = { version = "0.2", optional = true }
js-sys       = { version = "0.3", optional = true }
web-sys      = { version = "0.3", optional = true, features = [
    "Worker", "MessageChannel", "MessagePort",
    "SharedArrayBuffer", "Atomics",              # WASM 并行同步原语
    "DedicatedWorkerGlobalScope",
] }

# MPI（native 目标，feature = "mpi"）
rsmpi        = { version = "0.8", optional = true }

# 错误处理
thiserror = "1"
anyhow    = "1"

# 日志
log       = "0.4"
env_logger = "0.11"
```

---

## 7. Palace 配置兼容性要求

> **原则**: REM 是 Palace 的功能超集。任何合法的 Palace 配置文件，不加修改即可被 REM 解析并运行，结果误差在数值精度允许范围内与 Palace 一致。MoM/BEM 新增能力通过 REM 专有的扩展字段提供，这些字段在 Palace 中被忽略，因此双向无破坏性影响。

### 7.1 必须支持的特性

| 特性 | 说明 | 实现位置 |
|------|------|---------|
| **属性范围展开** | `"Attributes": "1,3-5,6"` → `[1,3,4,5,6]` | `config/src/preprocess.rs` |
| **C++ 注释剥除** | `// 单行` 和 `/* 多行 */` 注释不影响解析 | `config/src/preprocess.rs` |
| **PascalCase 键名** | 完全匹配 Palace 文档中的大小写（如 `MinFreq` 非 `min_freq`） | `config/src/schema.rs` |
| **默认值填充** | 未指定字段使用 Palace 文档规定的默认值，而非 Rust 类型默认值 | `config/src/schema.rs` |
| **长度单位换算** | `Model.L0` 作用于网格坐标缩放，不影响频率/时间等其他量纲 | `mesh/src/mesh_data.rs` |
| **混合类型属性** | `"Attributes"` 同时接受整数数组和范围字符串 | `config/src/schema.rs` |
| **可选节宽松解析** | `Boundaries`、`Domains` 等顶层节完全缺失时视为空，不报错 | `#[serde(default)]` |

### 7.2 MoM/BEM 扩展字段（Palace 无损兼容）

REM v1.0 在 `Problem.Type` 中新增以下值，Palace 工作流完全不受影响：

```json
{
  "Problem": {
    "Type": "MoM",
    "Verbose": 1,
    "Output": "./output"
  },
  "Model": {
    "Mesh": "sphere.msh",
    "L0": 1.0e-3
  },
  "Boundaries": {
    "PEC": { "Attributes": [1] }
  },
  "Solver": {
    "MoM": {
      "Equation":    "CFIE",
      "Basis":       "RWG",
      "FreqMin":     1.0e9,
      "FreqMax":     3.0e9,
      "FreqStep":    0.5e9,
      "Alpha":       0.5,
      "SingularTol": 1.0e-6
    }
  },
  "Postprocessing": {
    "RCS": {
      "PhiDeg":   [0, 90],
      "ThetaDeg": "0:5:180"
    },
    "NearField": {
      "Attributes": [10]
    }
  }
}
```

**新增 `Problem.Type` 枚举值**:

| 值 | 说明 | Palace 兼容性 |
|----|------|-------------|
| `"MoM"` | 矩量法（RWG 基函数，EFIE/CFIE/PMCHWT；ACA 加速）✅ 已实现 | Palace 不识别，不影响其工作流 |
| `"BEM"` | 边界元法（Laplace P0 pulse）✅ 已实现 | 同上 |
| `"SBR"` | SBR+ 高频射线追踪 + 物理光学（PO）✅ 已实现 | 同上 |

**新增 `Solver.MoM` 节**（Palace 忽略未知节，无破坏性）:

| 字段 | 类型 | 默认 | 说明 |
|------|------|------|------|
| `Equation` | string | `"CFIE"` | `"EFIE"` \| `"MFIE"` \| `"CFIE"` \| `"PMCHWT"` |
| `Basis` | string | `"RWG"` | `"RWG"` \| `"Pulse"` |
| `FreqMin` | float | 必填 | 起始频率 [Hz] |
| `FreqMax` | float | 必填 | 终止频率 [Hz] |
| `FreqStep` | float | 必填 | 频率步进 [Hz] |
| `Alpha` | float | `0.5` | CFIE 混合系数 α（0=EFIE，1=MFIE） |
| `SingularTol` | float | `1e-6` | 奇异积分收敛容限 |
| `FastSolver` | string | `"Direct"` | `"Direct"` \| `"ACA"` \| `"FMM"` |

**新增 `Solver.SBR` 节**（SBR+ 高频射线追踪，✅ 已实现）:

| 字段 | 类型 | 默认 | 说明 |
|------|------|------|------|
| `FreqMin` | float | 必填 | 起始频率 [Hz] |
| `FreqMax` | float | 必填 | 终止频率 [Hz] |
| `FreqStep` | float | `0` | 频率步进 [Hz]（0 = 单频） |
| `RayDensity` | float | `1e4` | 射线密度 [rays/m²]，用于多次弹射 |
| `MaxBounces` | int | `5` | 最大弹射次数（1 = 纯一次弹射 PO） |
| `WeightThresh` | float | `1e-4` | 射线能量截断阈值 |
| `TargetType` | string | `"PEC"` | `"PEC"` \| `"Dielectric"` |
| `ThetaInc` | float | `0.0` | 入射仰角 [°] |
| `PhiInc` | float | `0.0` | 入射方位角 [°] |
| `Polarization` | string | `"theta"` | `"theta"` \| `"phi"` \| `"x"` \| `"y"` \| `"z"` |

**SBR+ 配置示例**（[examples/sbr_sphere/sbr_sphere.json](examples/sbr_sphere/sbr_sphere.json)）:

```json
{
  "Problem": { "Type": "SBR", "Output": "output/sbr_sphere" },
  "Model":   { "Mesh": "examples/sbr_sphere/mesh/sphere.msh", "L0": 1.0 },
  "Boundaries": { "PEC": { "Attributes": [1] } },
  "Solver": {
    "SBR": {
      "FreqMin": 3.0e9, "FreqMax": 3.0e9,
      "RayDensity": 5000.0, "MaxBounces": 3,
      "ThetaInc": 0.0, "PhiInc": 0.0, "Polarization": "theta"
    }
  },
  "Postprocessing": {
    "RCS": { "ThetaDeg": "0:5:180", "PhiDeg": [0.0] }
  }
}
```

**新增 `Postprocessing` 节**:

| 字段 | 说明 |
|------|------|
| `RCS.PhiDeg` | 远场 φ 角列表（度） |
| `RCS.ThetaDeg` | 远场 θ 角范围，支持 `"0:5:180"` 格式 |
| `NearField.Attributes` | 近场计算面的物理组 ID |

### 7.3 配置兼容性保证矩阵

| 场景 | 行为 |
|------|------|
| Palace 配置 → REM 解析 | ✅ 完全兼容，结果等价 |
| REM MoM 配置 → Palace 解析 | ⚠️ Palace 仅解析已知字段，未知字段被忽略；Palace 不执行 MoM 求解（预期行为） |
| 混合配置（同时含 FEM + MoM 字段） | ✅ REM 按 `Problem.Type` 派发到对应求解器 |
| `Boundaries.PEC` 在 MoM 中的含义 | 与 FEM 相同：被 PEC 属性标记的面用于 MoM 阻抗矩阵装配 |

### 7.4 兼容性测试基准

**Level 1 — Palace 原生示例（必须通过）**:
- `examples/rings/rings.json` — 静电场，电容矩阵误差 < 1%
- `examples/coaxial/coaxial.json` — 同轴线静电，电容/长度误差 < 0.5%
- `examples/cavity/cavity.json` — 谐振腔特征模，f₀ 误差 < 0.1%
- `examples/cpw/cpw.json` — 共面波导频域驱动，S₁₁ 误差 < 0.5 dB

**Level 2 — MoM/BEM/SBR+ 新增验证（已通过）**:
- PEC 球体 MoM 散射 — RCS vs Mie 级数（kα≈3，误差 < 0.5 dB）✅
- Laplace P0 BEM — 平行板电容与 FEM 结果误差 < 1% ✅
- SBR+ PEC 球体单站 RCS @ 1 GHz（ka≈10.5）— vs Mie 误差 < 0.1 dB ✅

---

## 8. 输出格式规范

### 8.1 文件结构（Palace 兼容）

```
{output_dir}/
├── palace.json              # 回写的解析后配置（调试用）
├── postpro/
│   ├── domain-E.csv         # 域能量（逐步）
│   ├── port-V.csv           # 端口电压（驱动问题）
│   ├── port-I.csv           # 端口电流（驱动问题）
│   ├── S.csv                # S 参数（驱动问题）
│   └── eig.csv              # 特征频率（特征模问题）
└── paraview/
    ├── mesh.vtu             # 网格
    ├── solution-0.vtu       # 第 0 步场量
    └── solution.pvd         # ParaView 时间序列索引
```

### 8.2 CSV 格式

`domain-E.csv`:
```
Freq (GHz),E_field (J),H_field (J),Total_E (J)
1.0,1.23e-12,1.23e-12,2.46e-12
```

`eig.csv`:
```
Mode,Freq (GHz),Q Factor
1,5.123456,1234.5
2,7.891234,987.6
```

---

## 9. WASM API

```typescript
// JavaScript / TypeScript 接口

interface RemSolverOptions {
  configJson: string;    // Palace JSON 配置（字符串）
  meshData: Uint8Array;  // GMSH .msh 文件内容
  onProgress?: (msg: string) => void;
}

interface RemSolverResult {
  success: boolean;
  nodeCoords: Float64Array;    // [x0,y0,z0, x1,y1,z1, ...]
  potential: Float64Array;     // 节点标量场
  eField: Float64Array;        // 单元矢量场 [ex0,ey0,ez0, ...]
  csvOutputs: Record<string, string>;  // 文件名 → CSV 内容
  errorMessage?: string;
}

// Rust 侧 wasm-bindgen 暴露
#[wasm_bindgen]
pub struct RemWasmSolver { ... }

#[wasm_bindgen]
impl RemWasmSolver {
    pub fn new(config_json: &str) -> Result<RemWasmSolver, JsValue>;
    pub fn load_mesh(&mut self, msh_bytes: &[u8]) -> Result<(), JsValue>;
    pub fn solve(&mut self) -> Result<JsValue, JsValue>;  // 返回 JSON
    pub fn get_vtk(&self) -> Result<Vec<u8>, JsValue>;
}
```

---

## 10. 测试与验证策略

### 10.1 单元测试

| 测试 | 目标 | 验证方法 |
|------|------|---------|
| 平行板电容器 | εᵣ 均匀 | C = ε₀εᵣA/d 解析解 |
| 同轴线电容 | 圆柱坐标 | C = 2πεL/ln(b/a) 解析解 |
| 方形截面线圈 | 2D 静磁 | L(解析近似) vs. FEM |
| Palace 配置解析 | 格式兼容 | 解析已知 JSON 后字段值断言 |

### 10.2 收敛测试

对所有静态求解器要求 H¹ 误差 O(h²)（P1 元）：
- 4×4、8×8、16×16、32×32、64×64 网格序列
- L2 误差相对于解析解下降斜率 ≥ 1.9

### 10.3 兼容性测试

解析 Palace 官方示例 JSON/YAML，确保无解析错误，配置树字段值与预期一致。

---

## 11. Web Demo 规范（crates/yew-app/）

### 11.1 技术栈

| 层 | 技术 | 说明 |
|----|------|------|
| 框架 | Yew 0.21 (Rust) | 纯 Rust 前端框架，函数组件 |
| 构建 | Trunk | WASM 构建 + 开发热重载 |
| 求解器 | 直接链接 rem-* crates | 同进程调用，无 JS/WASM 边界开销 |
| 代码展示 | `<pre>` | 配置 JSON + Rust 源码展示 |
| 部署 | 纯静态文件 | `trunk build` 产物可直接托管 |

### 11.2 页面布局

```
┌─────────────────────────────────────────────────────────────┐
│  REM EM Solver Demo (Yew + WASM)                            │
├──────────────────────┬──────────────────────────────────────┤
│                      │  [Palace Config] [Test Source]       │
│  Example: [▼ 选择]   │                                      │
│                      │  {                                   │
│  [▶ Run Simulation]  │    "Problem": {                      │
│                      │      "Type": "Electrostatic"         │
│  Summary Result:     │    },                                │
│  Energy: 0.123 pJ    │    "Model": { ... },                │
│  Nodes: 456          │    "Boundaries": { ... },            │
│  Max |E|: 1.23 V/m   │    "Solver": { ... }                │
│                      │  }                                   │
│  Logs:               │                                      │
│  ┌────────────────┐  │                                      │
│  │Starting sim... │  │                                      │
│  │Completed.      │  │                                      │
│  └────────────────┘  │                                      │
└──────────────────────┴──────────────────────────────────────┘
```

### 11.3 架构

Yew 应用直接链接 rem-* Rust crates，在 WASM 主线程中调用求解器：

```rust
// crates/yew-app/src/solver.rs
pub fn run_example(key: &str) -> Result<SimResult, String> {
    let config = load_config_from_str(config_json, ConfigFormat::Json)?;
    let mesh = load_mesh_from_bytes(&config, &mesh_bytes, &NoComm)?;
    mesh.partition(&comm);
    // 根据 problem type 调用 solve_es / solve_ms ...
}
```

优势：
- 无 Worker/jsmpi mock — 单进程直接调用
- 编译期类型安全 — 配置/求解器接口在编译时校验
- 代码量更少 — 无 JS/TS、无 package.json、无 node_modules

### 11.4 预置示例

8 个 Palace 示例（定义在 `examples.rs`）：
- **Spheres** (Electrostatic) ✅ — 生成网格 `annular_msh`
- **Rings** (Magnetostatic) ✅ — 生成网格 `rect_bimaterial_msh`
- **Coaxial** (Electrostatic) ✅ — 内嵌 .msh 文件
- **Cylinder** (Magnetostatic) ✅ — 内嵌 .msh 文件
- **Adapter** (Driven) 🔲 v0.2
- **Antenna** (Driven) 🔲 v0.2
- **CPW** (Driven) 🔲 v0.2
- **Transmon** (Eigenmode) 🔲 v0.2

### 11.5 构建与部署

```bash
# 开发模式
cd crates/yew-app && trunk serve      # http://localhost:8080

# 构建静态产物
cd crates/yew-app && trunk build      # 输出到 dist/

# 部署：dist/ 内容直接托管至 GitHub Pages / Nginx 等
```

---

## 12. 并行计算架构

### 12.1 统一 Comm 抽象

`crates/parallel` 提供一个 `Comm` trait，同时支持三种后端，编译时通过 feature flag 选择：

| 后端 | 条件 | 说明 |
|------|------|------|
| `SerialComm` | 始终可用（fallback） | 单进程，`rank=0 size=1` |
| `MpiComm` | `features = ["mpi"]`，native 目标 | 包装 `rsmpi` |
| `WorkerComm` | `target_arch = "wasm32"` | 用 jsmpi + Web Worker 模拟 |

```rust
// crates/parallel/src/lib.rs

/// 统一并行通信抽象
pub trait Comm: Send + Sync {
    fn rank(&self) -> usize;
    fn size(&self) -> usize;

    // --- 同步 ---
    fn barrier(&self);

    // --- 集合操作 ---
    fn allreduce_sum_f64(&self, local: &[f64], global: &mut [f64]);
    fn allreduce_sum_usize(&self, local: &[usize], global: &mut [usize]);
    fn broadcast_bytes(&self, root: usize, data: &mut Vec<u8>);

    // --- 散播/收集 ---
    fn scatter_f64(&self, root: usize, send: Option<&[f64]>, recv: &mut [f64]);
    fn gather_f64(&self, root: usize, send: &[f64]) -> Option<Vec<f64>>;

    // --- 点对点 ---
    fn send_f64(&self, dest: usize, tag: u32, data: &[f64]);
    fn recv_f64(&self, src: usize, tag: u32, buf: &mut [f64]);
}

/// 运行时构建：根据编译目标和 feature 自动选择后端
pub fn build_comm() -> Box<dyn Comm> {
    #[cfg(feature = "mpi")]
    { return Box::new(mpi_comm::MpiComm::init()); }

    #[cfg(all(target_arch = "wasm32", feature = "wasm-parallel"))]
    { return Box::new(worker_comm::WorkerComm::init()); }

    Box::new(serial::SerialComm)
}
```

### 12.2 Native MPI 后端（rsmpi）

```rust
// crates/parallel/src/mpi_comm.rs
#[cfg(feature = "mpi")]
use mpi::traits::*;

#[cfg(feature = "mpi")]
pub struct MpiComm {
    universe: mpi::environment::Universe,
}

#[cfg(feature = "mpi")]
impl MpiComm {
    pub fn init() -> Self {
        let universe = mpi::initialize().expect("MPI init failed");
        MpiComm { universe }
    }
    fn world(&self) -> mpi::topology::SimpleCommunicator {
        self.universe.world()
    }
}

#[cfg(feature = "mpi")]
impl Comm for MpiComm {
    fn rank(&self) -> usize { self.world().rank() as usize }
    fn size(&self) -> usize { self.world().size() as usize }

    fn barrier(&self) { self.world().barrier(); }

    fn allreduce_sum_f64(&self, local: &[f64], global: &mut [f64]) {
        self.world().all_reduce_into(local, global, &mpi::collective::SystemOperation::sum());
    }

    fn broadcast_bytes(&self, root: usize, data: &mut Vec<u8>) {
        let root_proc = self.world().process_at_rank(root as i32);
        // 先广播长度，再广播内容
        let mut len = [data.len() as u64];
        root_proc.broadcast_into(&mut len);
        if self.rank() != root { data.resize(len[0] as usize, 0u8); }
        root_proc.broadcast_into(data.as_mut_slice());
    }
    // ... scatter/gather/send/recv 类似
}
```

**Cargo feature 配置**：

```toml
# crates/parallel/Cargo.toml
[features]
default = []
mpi = ["dep:rsmpi"]

[dependencies]
rsmpi = { version = "0.8", optional = true }
```

**用法**：
```bash
# Native 单进程（默认）
cargo run -p cli -- config.json

# Native MPI 4 进程
cargo build --features parallel/mpi -p cli
mpirun -np 4 ./target/release/rem config.json
```

### 12.3 WASM Web Worker 模拟 MPI（WorkerComm）

#### 核心原理

每个 Web Worker 对应一个 MPI rank，通信通过 `SharedArrayBuffer` + `Atomics` 实现：
- **barrier**：所有 Worker 对共享计数器执行 `Atomics.add`，然后 `Atomics.wait` 等待全部到达
- **allreduce**：每个 Worker 将局部值写入共享缓冲区对应槽位，barrier 后求和
- **点对点**：每对 (src, dst) 维护一个环形缓冲区，`Atomics.store`/`Atomics.load` 交换数据

> **前提**：服务器需设置 `Cross-Origin-Opener-Policy: same-origin` 和
> `Cross-Origin-Embedder-Policy: require-corp`（已在 `vite.config.ts` 中配置）

#### 共享内存布局

```
SharedArrayBuffer (64 MB)
┌──────────────────────────────────────────────────────────────┐
│  0..63       │ 控制区（i32×16）：barrier计数、锁等           │
│  64..16447   │ 数据交换区（f64×2048 per rank pair）          │
│  16448..end  │ 全局 allreduce 临时区（f64×MAX_DOFS）         │
└──────────────────────────────────────────────────────────────┘
```

#### JS 协调器（主线程）

```typescript
// web/src/worker/mpi-coordinator.ts

export interface MpiOptions {
  nWorkers?: number;       // 默认 navigator.hardwareConcurrency
  wasmUrl: string;         // rem_wasm.wasm 路径
  configJson: string;
  meshBytes: ArrayBuffer;
}

export async function runParallel(opts: MpiOptions): Promise<SolveResult> {
  const size = opts.nWorkers ?? Math.min(navigator.hardwareConcurrency, 8);
  const sharedBuf = new SharedArrayBuffer(64 * 1024 * 1024);

  const workers = Array.from({ length: size }, (_, rank) => {
    const w = new Worker(
      new URL('./solver-mpi.worker.ts', import.meta.url),
      { type: 'module' }
    );
    // 注入 rank/size 和共享内存
    w.postMessage({
      type: 'init',
      rank,
      size,
      sharedBuffer: sharedBuf,
      wasmUrl: opts.wasmUrl,
    });
    return w;
  });

  // 广播网格和配置（transferable ArrayBuffer 只能给一个 worker，其余复制）
  workers[0].postMessage({
    type: 'solve',
    configJson: opts.configJson,
    meshBytes: opts.meshBytes,   // rank 0 拥有原始数据
  }, [opts.meshBytes]);
  for (let r = 1; r < size; r++) {
    workers[r].postMessage({ type: 'solve', configJson: opts.configJson });
  }

  return new Promise((resolve, reject) => {
    workers[0].onmessage = (e) => {
      if (e.data.type === 'result') resolve(e.data.payload);
      else if (e.data.type === 'error') reject(new Error(e.data.payload));
    };
  });
}
```

#### Rust WorkerComm 实现

```rust
// crates/parallel/src/worker_comm.rs
#[cfg(all(target_arch = "wasm32", feature = "wasm-parallel"))]
use wasm_bindgen::prelude::*;

#[wasm_bindgen(module = "jsmpi")]
extern "C" {
    #[wasm_bindgen(js_name = MPI_Comm_rank)]
    fn mpi_comm_rank(comm: i32, rank: &mut i32) -> i32;
    #[wasm_bindgen(js_name = MPI_Comm_size)]
    fn mpi_comm_size(comm: i32, size: &mut i32) -> i32;
    #[wasm_bindgen(js_name = MPI_Barrier)]
    fn mpi_barrier(comm: i32) -> i32;
    #[wasm_bindgen(js_name = MPI_Allreduce)]
    fn mpi_allreduce(sendbuf: &f64, recvbuf: &mut f64, count: i32, datatype: i32, op: i32, comm: i32) -> i32;
    // ... 其他 jsmpi 绑定
}

#[cfg(all(target_arch = "wasm32", feature = "wasm-parallel"))]
pub struct WorkerComm;

impl Comm for WorkerComm {
    fn rank(&self) -> usize {
        let mut r = 0;
        unsafe { mpi_comm_rank(MPI_COMM_WORLD, &mut r); }
        r as usize
    }
    fn size(&self) -> usize {
        let mut s = 0;
        unsafe { mpi_comm_size(MPI_COMM_WORLD, &mut s); }
        s as usize
    }
    fn barrier(&self) {
        unsafe { mpi_barrier(MPI_COMM_WORLD); }
    }
    fn allreduce_sum_f64(&self, local: &[f64], global: &mut [f64]) {
        unsafe {
            mpi_allreduce(
                local.as_ptr(),
                global.as_mut_ptr(),
                local.len() as i32,
                MPI_DOUBLE,
                MPI_SUM,
                MPI_COMM_WORLD
            );
        }
    }
    // ... 包装 jsmpi 的其他接口
}
```

### 12.4 并行 FEM 装配模式

并行装配使用 **METIS 分区**（负载均衡，默认）或 **几何分区**（fallback）：

```
整体网格（rank 0 读取）
    ↓  RemMesh::partition()（rem-mesh，feature = "metis"）
       → METIS k-way 对偶图分区（rmetis，纯 Rust / WASM 兼容）
       → 失败时自动 fallback 到 X 轴几何分区
各 rank 持有本地单元子集 + 幽灵节点层（ghost layer）
    ↓  本地装配（fem-assembly）→ K_local, f_local
    ↓  allreduce_sum（重叠 DOF 的贡献求和）
    ↓  全局 K_global, f_global（仅 rank 0 组装完整矩阵）
         或分布式矩阵（直接用分布式求解器）
    ↓  求解（PCG+AMG，分布式版本）
    ↓  gather（场量收集到 rank 0）
    ↓  输出
```

```rust
// crates/mesh/src/mesh_data.rs
impl RemMesh {
    /// 统一分区入口：feature="metis" 时使用 METIS k-way，否则几何分区
    pub fn partition(&mut self) { ... }

    /// METIS 对偶图分区（仅 feature="metis" 编译）
    #[cfg(feature = "metis")]
    pub fn partition_metis(&mut self) -> Result<(), rmetis::MetisError> {
        // 1. 枚举每个体单元的面（3-D: 三角面；2-D: 边），按排序节点集 hash
        // 2. 共享面 → 对偶图邻接边，构造 CSR xadj/adjncy
        // 3. rmetis::part_graph_kway(graph, nparts, None, None, &Options::for_kway())
        // 4. 按 result.part 赋值 element.rank
        // 5. 边界单元按最近体单元质心分配 rank
    }
}
```

### 12.5 rmetis 子模块

rmetis（`vendor/rmetis`，git submodule）是纯 Rust 实现的 METIS 5.1.x 兼容图分区库：
- 无 C/系统依赖，可编译至 `wasm32-unknown-unknown`
- 提供 `part_graph_kway`、`part_graph_recursive`、`node_nd` 三个主要 API
- 通过 `rem-mesh` 的 `metis` feature 可选启用

### 12.6 feature flag 汇总

```toml
# crates/mesh/Cargo.toml
[features]
metis = ["dep:rmetis"]   # METIS 网格分区（纯 Rust，无系统依赖）

# crates/parallel
default  = []            # 串行 fallback，无依赖
mpi      = ["dep:rsmpi"] # native MPI（需系统安装 OpenMPI/MPICH）

# crates/wasm（自动激活 WorkerComm）
wasm            = ["dep:wasm-bindgen", "dep:js-sys"]
wasm-parallel   = ["wasm", "dep:web-sys", "dep:js-sys"]
```

**编译目标与后端对应关系**：

| 编译命令 | Comm 后端 | 并行度 |
|----------|-----------|--------|
| `cargo build` | `SerialComm` | 1 进程 |
| `cargo build --features parallel/mpi` | `MpiComm` | N 进程（mpirun） |
| `cargo build --target wasm32-unknown-unknown` | `SerialComm` | 1 Worker |
| `cargo build --target wasm32-unknown-unknown --features parallel/wasm-parallel` | `WorkerComm` | N Workers (via jsmpi stubs) |

---

## 13. 非功能需求

| 需求 | 指标 |
|------|------|
| 纯 Rust | 无 unsafe（除 wasm-bindgen 必需处）、无 C FFI |
| WASM 编译 | `cargo build --target wasm32-unknown-unknown` 无错 |
| 编译时间 | 全量编译 < 90 秒（i7/Ryzen 桌面级） |
| 内存占用 | 100K DOF 问题峰值 < 2 GB |
| 求解速度 | 100K DOF 静电 PCG < 30 秒（单线程） |
| 错误信息 | 所有用户可见错误含文件行号和修正提示 |
| 日志级别 | `RUST_LOG=rem=debug` 可复现完整求解过程 |
