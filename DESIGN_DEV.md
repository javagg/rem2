# REM — 设计与开发文档
## AI Agent 工作指南 v0.2

> 本文档为 AI Agent 提供分阶段开发 rem 电磁仿真工具的具体任务指令。
> 每个阶段均含：目标、关键文件、实现步骤、验收标准。
> 阅读本文档前请先阅读 `TECHNICAL_SPEC.md`。

---

## 0. 准备工作

### 0.1 工具链要求

```bash
# 必须安装
rustup target add wasm32-unknown-unknown
cargo install wasm-pack          # 可选，WASM 打包用
cargo install cargo-nextest      # 更快的测试运行器

# 验证
rustc --version   # >= 1.75
cargo --version
```

### 0.2 fem-rs 依赖确认

在开始编码前，克隆并编译 fem-rs 以确认依赖可用：

```bash
git clone https://github.com/javagg/fem-rs.git /tmp/fem-rs
cd /tmp/fem-rs
cargo test --workspace 2>&1 | tail -20
cargo build --target wasm32-unknown-unknown -p fem-wasm --no-default-features
```

若编译失败，记录错误并在 rem Cargo.toml 中使用 git rev 锁定可用 commit。

### 0.3 工作区初始化

```bash
cd c:/Users/lilu/works/rem2
# 若目录为空，初始化 workspace Cargo.toml
# 注意：不要运行 cargo new，直接创建文件
```

---

## 阶段 1: 工作区骨架 + 配置解析（优先级：最高）

**目标**: 建立 Cargo workspace，实现 Palace 配置文件的完整解析，通过所有兼容性测试。

### 1.1 创建 Cargo.toml（workspace 根）

创建文件 `c:/Users/lilu/works/rem2/Cargo.toml`：

```toml
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
    "crates/parallel",
    "crates/cli",
    "crates/wasm",
]

[workspace.dependencies]
# FEM 库（使用 path 依赖调试，之后换 git）
serde       = { version = "1", features = ["derive"] }
serde_json  = "1"
serde_yaml  = "0.9"
clap        = { version = "4", features = ["derive"] }
nalgebra    = "0.33"
num-complex = "0.4"
thiserror   = "1"
anyhow      = "1"
log         = "0.4"
env_logger  = "0.11"
wasm-bindgen = { version = "0.2", optional = true }

[profile.wasm]
inherits = "release"
opt-level = "s"
lto = true
```

### 1.2 创建 rem-core

**文件**: `crates/core/src/lib.rs`

```rust
pub mod constants {
    pub const EPS0: f64 = 8.854_187_817e-12;  // [F/m]
    pub const MU0: f64 = 1.256_637_061e-6;    // [H/m]
    pub const ETA0: f64 = 376.730_313;        // [Ω] 自由空间阻抗
    pub const C0: f64 = 2.997_924_58e8;       // [m/s] 光速
}

pub mod error {
    use thiserror::Error;
    
    #[derive(Debug, Error)]
    pub enum RemError {
        #[error("IO error: {0}")]
        Io(#[from] std::io::Error),
        #[error("Config parse error at {file}:{line}: {msg}")]
        Config { file: String, line: u32, msg: String },
        #[error("Mesh error: {0}")]
        Mesh(String),
        #[error("Solver did not converge after {max_iter} iterations (residual={residual:.2e})")]
        SolverNotConverged { max_iter: usize, residual: f64 },
        #[error("Unknown problem type: {0}")]
        UnknownProblemType(String),
        #[error("Unknown file format")]
        UnknownFormat,
    }
    
    pub type RemResult<T> = Result<T, RemError>;
}
```

### 1.3 创建 rem-config（关键模块）

**文件结构**:
```
crates/config/
├── Cargo.toml
└── src/
    ├── lib.rs
    ├── schema.rs       ← Palace 完整 schema
    ├── preprocess.rs   ← 注释剥除 + 属性范围展开
    └── defaults.rs     ← 所有默认值函数
```

**schema.rs 关键要点**：

所有结构体必须完全匹配 Palace 文档中的字段名（PascalCase），使用 `#[serde(rename_all = "PascalCase")]` 或逐字段 `#[serde(rename = "...")]`。

Palace 特殊规则：
1. `Attributes` 字段可以是整数数组 `[1,2,3]` 或字符串 `"1,3-5,6"`，需自定义 Deserializer
2. 所有浮点数支持科学记数法 `1.0e9`
3. 未设置的可选节（如 `Eigenmode`）应为 `Option<...>`

**preprocess.rs 实现要点**：

```rust
/// 剥除 // 和 /* */ 注释，保留字符串内容不变
pub fn strip_comments(s: &str) -> String { ... }

/// 展开 "1,3-5,6" 格式为 [1,3,4,5,6]
pub fn expand_ranges(s: &str) -> Vec<u32> { ... }

/// 自定义反序列化：接受 Vec<u32> 或字符串形式的 Attributes
pub fn deserialize_attributes<'de, D>(d: D) -> Result<Vec<u32>, D::Error>
where D: Deserializer<'de> { ... }
```

**lib.rs 公开 API**：

```rust
pub fn load_config(path: &std::path::Path) -> RemResult<PalaceConfig>;
pub fn load_config_from_str(content: &str, format: ConfigFormat) -> RemResult<PalaceConfig>;
pub enum ConfigFormat { Json, Yaml }
```

### 1.4 配置兼容性测试

创建 `tests/fixtures/` 目录，放入以下最小 Palace 兼容配置：

**tests/fixtures/parallel_plate.json**:
```json
{
  // Parallel plate capacitor test
  "Problem": {
    "Type": "Electrostatic",
    "Verbose": 1,
    "Output": "./output/parallel_plate"
  },
  "Model": {
    "Mesh": "tests/meshes/unit_square.msh",
    "L0": 1.0e-3
  },
  "Domains": {
    "Materials": [
      {
        "Attributes": [1],
        "Permittivity": 4.5,
        "LossTan": 0.02
      }
    ]
  },
  "Boundaries": {
    "Ground": { "Attributes": [1, 2] },
    "PEC":    { "Attributes": [3, 4] }
  },
  "Solver": {
    "Order": 1,
    "Linear": {
      "Type": "GMRES",
      "Tol": 1.0e-8,
      "MaxIter": 500
    }
  }
}
```

**验收标准**:
- [ ] `cargo test -p rem-config` 全部通过
- [ ] 解析上述 JSON 后，`config.problem.problem_type == ProblemType::Electrostatic`
- [ ] `config.domains.materials[0].permittivity == 4.5`
- [ ] 属性范围 `"1,3-5"` 被正确展开为 `[1,3,4,5]`
- [ ] C++ 注释被正确剥除（不影响字段值）
- [ ] YAML 格式与 JSON 格式解析结果完全等价

---

## 阶段 2: 网格与材料适配层

**目标**: 包装 fem-rs 的网格/IO，添加物理组到材料的映射。

### 2.1 rem-mesh 实现

依赖 fem-rs 的 `fem-io` crate。核心数据结构：

```rust
// crates/mesh/src/mesh_data.rs
pub struct RemMesh {
    /// fem-rs 原生网格
    pub inner: fem_mesh::SimplexMesh<3>,
    /// GMSH 物理组 ID → 材料索引
    pub domain_tags: std::collections::HashMap<u32, usize>,
    /// GMSH 物理组 ID → 边界条件类型
    pub boundary_tags: std::collections::HashMap<u32, BoundaryTag>,
}

pub enum BoundaryTag {
    Pec,
    Pmc,
    Impedance { rs: f64, ls: f64, cs: f64 },
    Ground,
    ZeroCharge,
    Absorbing { order: u8 },
    LumpedPort { index: u32, r: f64 },
}
```

**关键函数**：
```rust
pub fn load_mesh(config: &PalaceConfig) -> RemResult<RemMesh>;
```

此函数需：
1. 调用 `fem_io::read_msh_file()` 读取网格
2. 遍历 `config.domains.materials` 建立 `domain_tags` 映射
3. 遍历 `config.boundaries.*` 建立 `boundary_tags` 映射
4. 应用 `config.model.l0` 缩放所有坐标

### 2.2 rem-materials 实现

```rust
// crates/materials/src/material.rs
#[derive(Debug, Clone)]
pub struct Material {
    pub permittivity: f64,   // εᵣ
    pub permeability: f64,   // μᵣ
    pub conductivity: f64,   // σ [S/m]
    pub loss_tangent: f64,   // tan δ
}

impl Material {
    pub fn epsilon_eff(&self, freq: f64) -> num_complex::Complex64 {
        // εᵣ(1 - j·tan_δ) - jσ/(ωε₀)
        let eps_r = num_complex::Complex64::new(self.permittivity, -self.permittivity * self.loss_tangent);
        if freq > 0.0 {
            let sigma_term = self.conductivity / (2.0 * std::f64::consts::PI * freq * EPS0);
            eps_r - num_complex::Complex64::new(0.0, sigma_term)
        } else {
            eps_r
        }
    }
    
    /// 每个积分点的 ε 系数（静电场用）
    pub fn epsilon_scalar(&self) -> f64 {
        self.permittivity
    }
    
    /// 每个积分点的 ν = 1/μ 系数（静磁场用）
    pub fn reluctivity(&self) -> f64 {
        1.0 / (MU0 * self.permeability)
    }
}
```

**验收标准**:
- [ ] 成功加载 GMSH .msh 文件，物理组正确映射到材料
- [ ] 材料属性查询函数单元测试通过
- [ ] 坐标缩放（L0 = 1e-3 时坐标乘以 0.001）正确

---

## 阶段 3: 静电场求解器（v0.1 核心功能）

**目标**: 实现完整静电场求解流程，通过平行板和同轴线解析解验证。

### 3.1 实现流程（electrostatic/src/solver.rs）

```
load_mesh(config)
    ↓
H1Space::new(mesh, P1)        ← fem-space
    ↓
assemble_stiffness(space, ε)  ← fem-assembly DiffusionIntegrator
assemble_load(space, ρ)       ← fem-assembly DomainSourceIntegrator（通常 ρ=0）
    ↓
apply_dirichlet_bc(K, f, pec_dofs, 0.0)   ← PEC: φ=0
apply_dirichlet_bc(K, f, gnd_dofs, 0.0)   ← Ground: φ=0
apply_dirichlet_bc(K, f, hot_dofs, V)     ← 激励端口: φ=V
apply_neumann_bc(f, zero_charge_faces, 0.0) ← 齐次 Neumann（默认）
    ↓
PCG+AMG(K, f, tol, max_iter)  ← fem-solver
    ↓
phi: Vec<f64>                  ← 节点电位
    ↓
e_field = gradient_recovery(phi, space)  ← 后处理
capacitance_matrix = compute_C(phi, space, ports)
    ↓
write_vtk(output_dir, phi, e_field, mesh)
write_csv_energy(output_dir, ...)
```

### 3.2 DiffusionIntegrator 变系数使用

fem-rs 的 `DiffusionIntegrator` 接受 `Fn(f64, f64, f64) -> f64` 系数函数。需要从材料图构建此函数：

```rust
fn build_epsilon_fn(
    mesh: &RemMesh, 
    materials: &[Material]
) -> impl Fn(f64, f64, f64) -> f64 + '_ {
    move |x, y, z| {
        let elem_id = mesh.find_element_at(x, y, z);
        let mat_idx = mesh.domain_tags[&elem_id];
        materials[mat_idx].epsilon_scalar()
    }
}
```

**注意**: fem-rs 可能尚未提供 `find_element_at` 的空间查询。替代方案：
- 在组装循环中，已知当前单元 ID，直接用单元 ID 查材料（需改造 Assembler 接口或使用 per-element 系数）
- 构建 `elem_id → epsilon` 的 `Vec<f64>` 查找表，传入 per-element 系数版本的 DiffusionIntegrator

### 3.3 梯度恢复（E 场后处理）

P1 元的梯度在每个单元内为常数，需做节点平均（ZZ 恢复）：

```rust
pub fn gradient_recovery(phi: &[f64], space: &H1Space) -> Vec<[f64; 3]> {
    let mesh = space.mesh();
    let mut e_field = vec![[0.0f64; 3]; mesh.n_nodes()];
    let mut counts  = vec![0usize; mesh.n_nodes()];
    
    for elem in 0..mesh.n_elements() {
        let nodes = mesh.element_nodes(elem);
        let grad = compute_p1_gradient(phi, nodes, mesh);  // 常数梯度
        for &node in nodes {
            e_field[node][0] -= grad[0];  // E = -∇φ
            e_field[node][1] -= grad[1];
            e_field[node][2] -= grad[2];
            counts[node] += 1;
        }
    }
    
    for (i, (e, c)) in e_field.iter_mut().zip(counts.iter()).enumerate() {
        if *c > 0 {
            e[0] /= *c as f64;
            e[1] /= *c as f64;
            e[2] /= *c as f64;
        }
    }
    e_field
}
```

### 3.4 验证测试

**测试 1: 平行板电容器**

- 网格: 单位正方形，顶边 tag=1, 底边 tag=2, 左右边 tag=3
- 配置: 顶边 PEC φ=1V, 底边 Ground φ=0, 左右 ZeroCharge
- 解析解: φ(x,y) = y, E_y = -1 V/m
- 验收: L2 误差 < 1e-10（P1 在线性解上精确）

**测试 2: 同轴线（GMSH 网格）**

- 网格: 同轴截面，内圆 r_i=0.5mm, 外圆 r_o=2mm
- 解析解: φ(r) = ln(r/r_o) / ln(r_i/r_o)
- 验收: L2 误差 < 0.5% (P1 in GMSH mesh)

**验收标准**:
- [ ] 平行板解析解 L2 误差 < 1e-8
- [ ] 同轴线 L2 误差 O(h²) 收敛（4 级网格细化）
- [ ] VTK 输出文件可用 ParaView 打开
- [ ] domain-E.csv 格式与 Palace 兼容（列名完全一致）

---

## 阶段 4: 静磁场求解器

**目标**: 实现 2D (A_z) 和 3D (A 矢量位) 静磁场求解。

### 4.1 2D 静磁（优先实现）

方程: `-∇·(ν ∇A_z) = J_z`

实现与静电场完全类似，仅：
- 系数从 ε 换为 ν = 1/(μ₀μᵣ)
- 右端项为电流密度 J_z 而非电荷密度 ρ
- 后处理: B_x = ∂A_z/∂y, B_y = -∂A_z/∂x（P1 梯度恢复）
- 电感矩阵: L = Φ/I = (∫A_z J_z dΩ) / I²

### 4.2 3D 静磁（需 Nedelec 元，依赖 fem-rs Phase 5）

**如果 fem-rs 已实现 Nedelec 元**：直接使用 `NedelecSpace`（H(curl)）

**如果未实现**（当前状态）：
- 暂用 H1 节点元近似（仅适用于简单几何）
- 在代码注释中标注 `// TODO: upgrade to Nedelec when fem-rs Phase 5 ships`
- 提供 `--2d` 标志仅使用 2D A_z 公式

### 4.3 验证测试

**测试: 方形截面直导线**
- J_z = 1 MA/m² 在中心 0.2×0.2 mm²
- 边界: A_z = 0 on Γ_D
- 验证: B 场围绕导线的正确方向和量级

**验收标准**:
- [ ] 2D 静磁正确收敛，B 场方向与分析预期一致
- [ ] 变磁导率（μᵣ=1000 铁芯）材料正确处理
- [ ] 3D 情形有明确错误消息或 TODO 标注

---

## 阶段 5: CLI 入口与输出模块

### 5.1 rem-cli 实现

```rust
// crates/cli/src/main.rs
use clap::Parser;

#[derive(Parser)]
#[command(name = "rem", about = "Rust Electromagnetic Solver")]
struct Args {
    /// Palace 格式配置文件 (.json 或 .yaml)
    config: std::path::PathBuf,
    
    /// 覆盖输出目录
    #[arg(short, long)]
    output: Option<std::path::PathBuf>,
    
    /// 详细日志（重复使用增加详细度：-v -vv -vvv）
    #[arg(short, long, action = clap::ArgAction::Count)]
    verbose: u8,
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();
    
    // 设置日志级别
    let log_level = match args.verbose {
        0 => "warn",
        1 => "info",
        2 => "debug",
        _ => "trace",
    };
    env_logger::Builder::from_env(
        env_logger::Env::default().default_filter_or(log_level)
    ).init();
    
    let mut config = rem_config::load_config(&args.config)?;
    if let Some(output) = args.output {
        config.problem.output = output.to_string_lossy().into_owned();
    }
    
    match config.problem.problem_type {
        ProblemType::Electrostatic  => rem_electrostatic::run(&config, comm.as_ref())?,
        ProblemType::Magnetostatic  => rem_magnetostatic::run(&config, comm.as_ref())?,
        ProblemType::Eigenmode      => rem_eigenmode::run(&config, comm.as_ref())?,
        ProblemType::Driven         => rem_driven::run(&config, comm.as_ref())?,
        ProblemType::Transient      => {
            eprintln!("Transient solver not yet implemented (v1.0)");
            std::process::exit(1);
        }
    }
    
    Ok(())
}
```

### 5.2 result 实现

**VTK 输出（包装 fem-io）**：

```rust
pub struct SolutionData<'a> {
    pub mesh: &'a RemMesh,
    pub scalar_fields: Vec<(&'static str, &'a [f64])>,    // (name, node_values)
    pub vector_fields: Vec<(&'static str, &'a [[f64; 3]])>, // (name, element_vectors)
}

pub fn write_vtk(output_dir: &Path, data: &SolutionData) -> RemResult<()> {
    std::fs::create_dir_all(output_dir)?;
    let vtu_path = output_dir.join("paraview/solution.vtu");
    // 调用 fem_io::VtkWriter
    let mut writer = fem_io::VtkWriter::new(&vtu_path)?;
    writer.write_mesh(&data.mesh.inner)?;
    for (name, values) in &data.scalar_fields {
        writer.add_nodal_scalar(name, values)?;
    }
    for (name, vectors) in &data.vector_fields {
        writer.add_cell_vector(name, vectors)?;
    }
    writer.finish()?;
    Ok(())
}
```

**CSV 输出（Palace 格式）**：

```rust
pub fn write_domain_energy(
    output_dir: &Path,
    freq_or_step: f64,
    e_energy: f64,
    h_energy: f64,
) -> RemResult<()> {
    let path = output_dir.join("postpro/domain-E.csv");
    let header = "Freq (GHz),E_field (J),H_field (J),Total_E (J)";
    let line = format!("{:.6e},{:.6e},{:.6e},{:.6e}",
        freq_or_step, e_energy, h_energy, e_energy + h_energy);
    append_csv(&path, header, &line)
}
```

**验收标准**:
- [ ] `rem path/to/config.json` 命令正常运行
- [ ] 未实现的问题类型给出有意义的错误消息
- [ ] VTU 文件可在 ParaView 中打开并显示正确场量
- [ ] CSV 文件列名与 Palace 输出完全一致

---

## 阶段 6: WASM 绑定（crates/wasm）

**目标**: 将静电/静磁求解器封装为 WASM 模块，可在浏览器调用。

### 6.1 实现要点

```rust
// crates/wasm/src/lib.rs
#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
pub struct RemSolver {
    config: Option<rem_config::PalaceConfig>,
    mesh: Option<rem_mesh::RemMesh>,
}

#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
impl RemSolver {
    #[cfg_attr(target_arch = "wasm32", wasm_bindgen(constructor))]
    pub fn new() -> Self {
        // 初始化 panic hook（WASM 调试）
        #[cfg(target_arch = "wasm32")]
        console_error_panic_hook::set_once();
        RemSolver { config: None, mesh: None }
    }
    
    /// 从 JSON 字符串加载配置
    pub fn load_config(&mut self, json: &str) -> Result<(), String> {
        self.config = Some(
            rem_config::load_config_from_str(json, rem_config::ConfigFormat::Json)
                .map_err(|e| e.to_string())?
        );
        Ok(())
    }
    
    /// 从字节数组加载 GMSH .msh
    pub fn load_mesh(&mut self, data: &[u8]) -> Result<(), String> {
        let config = self.config.as_ref().ok_or("Config not loaded")?;
        self.mesh = Some(
            rem_mesh::load_mesh_from_bytes(config, data)
                .map_err(|e| e.to_string())?
        );
        Ok(())
    }
    
    /// 执行求解，返回 JSON 结果
    pub fn solve(&self) -> Result<String, String> {
        let config = self.config.as_ref().ok_or("Config not loaded")?;
        let mesh = self.mesh.as_ref().ok_or("Mesh not loaded")?;
        
        let result = match config.problem.problem_type {
            ProblemType::Electrostatic => {
                rem_electrostatic::solve_to_json(config, mesh)
            }
            ProblemType::Magnetostatic => {
                rem_magnetostatic::solve_to_json(config, mesh)
            }
            _ => Err(rem_core::error::RemError::UnknownProblemType(
                format!("{:?} not yet supported in WASM", config.problem.problem_type)
            ))
        };
        
        result.map_err(|e| e.to_string())
    }
}
```

### 6.2 WASM 编译限制处理

**必须避免的**:
- `std::fs::*`（文件系统 API）
- `std::process::*`（进程 API）
- 任何依赖系统时钟的计时（`std::time::Instant` 在 wasm32 可用）
- MPI（编译时用 `#[cfg(not(target_arch = "wasm32"))]` 隔离）

**crates/mesh 的 WASM 适配**:

```rust
// 需要添加从字节加载的版本
pub fn load_mesh_from_bytes(config: &PalaceConfig, data: &[u8]) -> RemResult<RemMesh> {
    // fem_io 的 GMSH 解析器需要支持从 &[u8] 读取
    // 若 fem-io 只有 read_msh_file，先写入临时内存缓冲区
    let content = std::str::from_utf8(data).map_err(|_| RemError::Mesh("Invalid UTF-8".into()))?;
    fem_io::read_msh_str(content)  // 如果 fem-io 没有此函数，需向 fem-rs 提 PR
        .map_err(|e| RemError::Mesh(e.to_string()))
        .and_then(|inner| RemMesh::from_fem_mesh(inner, &config.domains.materials))
}
```

### 6.3 构建与测试

```bash
# 构建 WASM
cargo build --target wasm32-unknown-unknown -p rem-wasm --no-default-features --features wasm

# 使用 trunk 构建 Yew 前端（输出到 crates/yew-app/dist/）
cd crates/yew-app && trunk build

# 测试 WASM 功能（在 node 中运行）
wasm-pack test crates/wasm --node
```

**验收标准**:
- [ ] `cargo build --target wasm32-unknown-unknown -p rem-wasm` 无错
- [ ] 从 JavaScript 调用 `RemSolver.new()`, `load_config()`, `load_mesh()`, `solve()` 返回正确结果
- [ ] WASM bundle 大小 < 5 MB（--release + opt-level=s）

---

## 阶段 7: 特征模与频域求解器（v0.2）[COMPLETED]

> **前置条件**: 阶段 1-5 完成
> **状态**: 已实现基础版本，所有示例可运行

### 7.1 特征模求解器 (rem-eigenmode)

**方程**: `K x = λ M x`（广义特征值问题）

**实际实现**:
1. P1 FEM 组装刚度矩阵 K 和一致质量矩阵 M（支持 Tri3/Tet4/Tet10）
2. Lanczos shift-invert 迭代：`(K - σM)^{-1} M v` 使用 PCG 内层求解
3. nalgebra `SymmetricEigen` 分解 m×m 三对角矩阵得到 Ritz 特征值
4. σ = (2πf_target/c)²，从 Palace config `Solver.Eigenmode.Target` 读取

**输出**:
- `eigenfrequencies.csv`: m, f (Hz)
- `mode_N.vtk`: 每个模态的标量场 VTK 输出

**已知限制**:
- 标量 P1 基函数（非 Nedelec 矢量元），适用于 TEM/quasi-static 问题
- Tet10 使用 corner-only P1 近似
- Ritz 向量恢复使用 Lanczos 基向量近似（非完整 Ritz 向量重建）

### 7.2 频域驱动求解器 (rem-driven)

**方程**: `[K - k₀² M] φ = f`（标量实数波动方程）

**实际实现**:
- 频率扫描 [MinFreq, MaxFreq] 步进 FreqStep
- 每频率点组装 A = K - k²M，PCG 求解
- LumpedPort 边界条件：激励端口 φ=1V（Dirichlet）
- 端口 V/I 计算：V = mean(φ_port)，I = Σ (K·φ)_port
- S₁₁ = (Z - Z₀)/(Z + Z₀)

**输出**:
- `port-S.csv`: f (Hz), Re(S11), Im(S11), |S11| (dB)
- `driven_NNNN.vtk`: 场分布快照（按 SaveStep）

**已知限制**:
- 实数算法（无复数导纳/损耗支持，Im(S11)=0）
- 无 PML（Absorbing BC 处理为 Dirichlet φ=0）
- 当 k² > λ_min 时系统失去 SPD，PCG 可能不收敛

**验收标准** (v0.2):
- [x] transmon 示例：特征模求解输出 5 个模态频率（GHz 范围）
- [x] cpw 示例：41 频率点驱动求解完成，S 参数 CSV 输出
- [x] adapter / antenna 示例：驱动求解完成（adapter 有 PCG 收敛警告，属预期）
- [x] S 参数 CSV 输出与 Palace 格式兼容

---

## 阶段 6.5: 并行计算层（crates/parallel） [COMPLETED]

**目标**: 实现 `Comm` trait 的三种后端，并将并行装配接入静电/静磁求解器。

> **前置条件**: 阶段 1-5 完成（串行求解器可用）

### 6.5.1 crates/parallel 实现顺序

**Step 1 — 先实现 SerialComm（无依赖，立即可测试）**

```rust
// crates/parallel/src/serial.rs
pub struct SerialComm;

impl Comm for SerialComm {
    fn rank(&self) -> usize { 0 }
    fn size(&self) -> usize { 1 }
    fn barrier(&self) {}
    fn allreduce_sum_f64(&self, local: &[f64], global: &mut [f64]) {
        global.copy_from_slice(local);   // 单进程：直接复制
    }
    fn broadcast_bytes(&self, _root: usize, _data: &mut Vec<u8>) {}
    fn scatter_f64(&self, _root: usize, send: Option<&[f64]>, recv: &mut [f64]) {
        recv.copy_from_slice(send.unwrap());
    }
    fn gather_f64(&self, _root: usize, send: &[f64]) -> Option<Vec<f64>> {
        Some(send.to_vec())
    }
    fn send_f64(&self, _dest: usize, _tag: u32, _data: &[f64]) {
        panic!("SerialComm: send to other rank impossible")
    }
    fn recv_f64(&self, _src: usize, _tag: u32, _buf: &mut [f64]) {
        panic!("SerialComm: recv from other rank impossible")
    }
}
```

**Step 2 — MpiComm（feature = "mpi"，native 专用）**

参照 TECHNICAL_SPEC.md §12.2，包装 rsmpi API。构建前确认系统 MPI：
```bash
which mpicc && mpicc --version   # OpenMPI 或 MPICH 均可
cargo add rsmpi --optional -p parallel
```

验证编译：
```bash
cargo build -p parallel --features mpi
cargo test  -p parallel --features mpi -- --test-threads=1
```

测试用例（单机 4 进程）：
```bash
mpirun -np 4 cargo test -p parallel --features mpi --test test_allreduce
```

**Step 3 — WorkerComm（target = wasm32）**

参照 TECHNICAL_SPEC.md §12.3。关键实现细节：

1. **`from_global()` 中的 JS 全局变量注入**：JS coordinator 在 `postMessage({ type: 'init' })` 前需先在 Worker 全局设置：
   ```typescript
   // solver-mpi.worker.ts 顶部
   self.addEventListener('message', async (e) => {
     if (e.data.type === 'init') {
       (self as any)._MPI_RANK = e.data.rank;
       (self as any)._MPI_SIZE = e.data.size;
       (self as any)._MPI_SAB  = e.data.sharedBuffer;
       // 然后初始化 WASM
       await init(e.data.wasmUrl);
       self.postMessage({ type: 'ready' });
     }
   });
   ```

2. **`Atomics.wait` 在主线程不可用**：所有 WorkerComm 操作必须在 Worker 线程内执行，禁止在主线程调用 `barrier()`。

3. **SharedArrayBuffer size 计算**（保守估计，可运行时动态调整）：
   ```rust
   // CTRL 区：64 × i32（256 bytes）
   // DATA 区：size × max_local_dofs × f64
   // 建议 64 MB 上限，超过则退化为串行
   const SAB_SIZE: usize = 64 * 1024 * 1024;
   ```

4. **Atomics 的 float 问题**：WebAssembly Atomics 仅支持 `i32/i64`，不支持直接原子操作 f64。
   解决方案：使用 `i64::from_bits(f64::to_bits(v))` 转换后通过 `Atomics.store`（`BigInt64Array`）写入。

验证编译：
```bash
cargo build -p parallel --target wasm32-unknown-unknown --features wasm-parallel
```

### 6.5.2 网格分区（mesh_data.rs）

METIS 分区已通过 `rmetis`（`vendor/rmetis` git submodule，纯 Rust）实现并默认启用：

```rust
// crates/mesh/src/mesh_data.rs
impl RemMesh {
    /// 统一入口：feature="metis" → METIS k-way；否则几何分区
    pub fn partition(&mut self) { ... }

    /// METIS 对偶图 k-way 分区
    #[cfg(feature = "metis")]
    pub fn partition_metis(&mut self) -> Result<(), rmetis::MetisError> {
        // 构造元素对偶图（共享面 → 邻接边），调用 rmetis::part_graph_kway
        // 边界单元按最近体单元质心分配 rank
    }

    /// 几何分区（X 轴均分），作为 fallback
    pub fn partition_geometric(&mut self) { ... }
}
```

启用方式（rem-cli/rem-wasm 的 Cargo.toml）：

```toml
rem-mesh = { workspace = true, features = ["metis"] }
```

### 6.5.3 接入静电场求解器

修改 `crates/electrostatic/src/solver.rs` 的 `run()` 签名：

```rust
// 原来
pub fn run(config: &PalaceConfig) -> RemResult<()>

// 改为
pub fn run(config: &PalaceConfig, comm: &dyn Comm) -> RemResult<()> {
    let mesh = load_mesh(config)?;
    let partition = partition_mesh(&mesh, comm);

    // 每个 rank 只装配本地单元
    let local_K = assemble_local_stiffness(&partition, &mesh, config)?;
    let local_f = assemble_local_rhs(&partition, &mesh, config)?;

    // allreduce 合并全局矩阵（仅适用于小问题；大问题用分布式求解器）
    let (K, f) = if comm.size() == 1 {
        (local_K, local_f)
    } else {
        allreduce_csr(local_K, local_f, comm)?
    };

    // rank 0 求解并输出
    if comm.rank() == 0 {
        let phi = solve_linear(K, f, &config.solver)?;
        write_results(config, &phi, &mesh)?;
    }
    Ok(())
}
```

在 `crates/cli/src/main.rs` 中：

```rust
let comm = parallel::build_comm();
match config.problem.problem_type {
    ProblemType::Electrostatic => electrostatic::run(&config, comm.as_ref())?,
    // ...
}
```

### 6.5.4 验收标准

- [ ] `SerialComm` 通过所有现有串行测试（不破坏已有功能）
- [ ] `MpiComm`（2 进程）对 64×64 网格静电场的结果与串行解 L2 误差 < 1e-10
- [ ] `WorkerComm`（2 Worker）WASM 构建无编译错误
- [ ] `build_comm()` 在三种编译条件下分别返回正确后端
- [ ] `mpirun -np 1/2/4 rem config.json` 结果一致（无数值差异）

---

## 阶段 7.5: Web Demo（crates/yew-app/）

**目标**: 构建纯 Rust Yew 前端，直接链接 rem-* 求解器，可零配置部署至 GitHub Pages。

> **前置条件**: 阶段 6（WASM 绑定）完成，`trunk` 已安装。

### 7.5.1 技术栈

- **Yew 0.21** — Rust 前端框架，函数组件 + hooks
- **Trunk** — WASM 构建工具，开发热重载
- **直接链接 rem-* crates** — 无 Worker/jsmpi，同进程调用求解器

### 7.5.2 项目结构

```
crates/yew-app/
├── Cargo.toml         # 依赖 yew, rem-*, wasm-bindgen 等
├── Trunk.toml         # Trunk 构建配置
├── index.html         # 入口 HTML
└── src/
    ├── main.rs        # Yew App 组件（示例选择、运行、结果、日志、代码查看）
    ├── examples.rs    # 8 个 Palace 示例配置 + 网格数据
    ├── solver.rs      # 封装 rem-* 求解器调用
    └── style.css      # 样式
```

### 7.5.3 求解器调用

Yew 组件通过 `wasm_bindgen_futures::spawn_local` 异步调用：

```rust
// crates/yew-app/src/solver.rs
pub fn run_example(key: &str) -> Result<SimResult, String> {
    let config = load_config_from_str(config_json, ConfigFormat::Json)?;
    let mesh = load_mesh_from_bytes(&config, &mesh_bytes, &NoComm)?;
    mesh.partition(&comm);
    match cfg.problem.problem_type {
        ProblemType::Electrostatic => { solve_es(...) }
        ProblemType::Magnetostatic => { solve_ms(...) }
        _ => Err("Not implemented"),
    }
}
```

### 7.5.4 构建与部署

```bash
# 开发模式
cd crates/yew-app && trunk serve      # http://localhost:8080

# 构建静态产物
cd crates/yew-app && trunk build      # 输出到 dist/
```

### 7.5.5 验收标准

- [x] `trunk serve` 启动无错误，浏览器显示完整 UI
- [x] 选择 Spheres 示例运行，显示能量/节点数/场强结果
- [x] 选择 Rings 示例运行磁静场
- [x] Config/Source tab 切换正常
- [x] 未实现的示例（Driven/Eigenmode）按钮禁用并显示提示

---

## 阶段 8: 全面 Palace 兼容测试

**目标**: 与 Palace 官方示例的数值结果对比。

### 8.1 测试矩阵

| 算例 | 问题类型 | 指标 | 目标误差 |
|------|----------|------|---------|
| `rings` | Electrostatic | 电容矩阵 C₁₁ | < 1% |
| `coaxial` | Electrostatic | 电容/长度 | < 0.5% |
| `cavity` | Eigenmode | f₀ (TM₀₁₀) | < 0.1% |
| `cpw` | Driven | S₁₁, S₂₁ @5GHz | < 0.5 dB |
| `inductance` | Magnetostatic | L 矩阵 | < 2% |

### 8.2 自动回归测试脚本

```bash
#!/bin/bash
# tests/regression/run_all.sh
for case in parallel_plate coaxial rings; do
  echo "=== $case ==="
  ./target/release/rem tests/fixtures/${case}.json --output /tmp/${case}_out
  python3 tests/compare_output.py /tmp/${case}_out tests/expected/${case}
done
```

---

## 开发注意事项

### A. fem-rs 集成陷阱

1. **Nedelec 元未实现**: 频域/特征模问题需要 H(curl) 空间，fem-rs Phase 5 尚未完成。
   - 临时方案: 用节点元（精度较差）+ 明确文档说明限制
   - 长期方案: 向 fem-rs 贡献 Nedelec 元实现，或在 rem 中自行实现

2. **GMSH 文件仅从路径读取**: `fem_io::read_msh_file` 只接受文件路径，WASM 中无文件系统。
   - 需要在 rem-mesh 中实现 `read_msh_bytes(&[u8])` 包装器
   - 若 fem-io 不支持，考虑内嵌简化 GMSH 解析器

3. **MPI 隔离**: `crates/parallel` 的 `MpiComm` 通过 `features = ["mpi"]` 条件编译，`WorkerComm` 仅在 `target_arch = "wasm32"` 下编译。
   - WASM 构建默认使用 `SerialComm`（单 Worker），需要多 Worker 并行时加 `--features parallel/wasm-parallel`
   - WASM 目标上 `rsmpi` 绝对不能进入依赖树（其依赖 C MPI 库）；`crates/parallel/Cargo.toml` 中 `rsmpi` 必须设 `optional = true`

4. **复数支持**: fem-rs 当前全部使用 `f64`，频域问题需要 `Complex64`。
   - 解决方案 A: 使用实部/虚部分离系统（倍维度实系统）
   - 解决方案 B: 自行实现复数版本的 CsrMatrix 和 Assembler

### B. Palace 格式边缘情况

1. **整数数组 vs 字符串属性**: `"Attributes": [1,2,3]` 和 `"Attributes": "1,3-5"` 必须都支持
2. **缺失节的处理**: `Boundaries` 可能完全缺少（空边界情况），需要 `#[serde(default)]`
3. **嵌套 vs 扁平**: Palace 某些版本的配置格式在不同示例中有细微差异，需宽松解析

### C. 数值稳定性

1. **病态矩阵**: 高对比度材料（εᵣ=1 vs εᵣ=10000）可能导致条件数极大
   - 应用 AMG 预条件器
   - 输出条件数估计（`log::warn` 当条件数 > 1e10）

2. **精度一致性**: 与 Palace 数值对比时，注意网格和单位的一致性
   - Palace 默认长度单位可能与用户网格不一致，确保 `L0` 正确应用

### D. 代码质量要求

- 所有 `pub` 函数必须有 rustdoc 注释
- `clippy --deny warnings` 编译必须无警告
- 所有错误路径必须通过 `RemError` 类型传播，禁止 `.unwrap()` 在库代码中
- 测试覆盖率目标：crates/config ≥ 90%，crates/electrostatic ≥ 80%

---

## 版本里程碑

| 版本 | 内容 | 对应阶段 | 状态 |
|------|------|---------|------|
| v0.1.0 | 工作区 + 配置解析 + 静电/静磁 + CLI + VTK 输出 | 1-5 | ✅ |
| v0.1.1 | WASM 绑定 + Yew Web Demo（可部署） | 6, 7.5 | ✅ |
| v0.2.0 | 并行层 (jsmpi/Comm trait) + 分布式组装 | 6.5, 7 | ✅ |
| v0.2.1 | rmetis 子模块 + METIS k-way 对偶图分区 | — | ✅ |
| v0.3.0 | Palace 官方示例完整兼容测试 | 8 | ✅ |
| v0.4.0 | MoM 基础设施：密集复数矩阵 + 表面网格提取 + 高斯求积 | 9 | ✅ |
| v0.5.0 | MoM EFIE 求解器：Green 函数 + 脉冲基函数 + 验证 | 10 | ✅ |
| v0.6.0 | MoM CFIE 求解器：RWG 基函数 + 奇异积分 + RCS 输出 | 11 | ✅ |
| v0.7.0 | BEM 静态求解器：Laplace BEM + 与 FEM 交叉验证 | 12 | ✅ |
| **v0.8.0** | **SBR+ 高频 PO 求解器：BVH + 两阶段 PO + Mie 验证** | **—** | **✅** |
| **v0.8.1** | **警告清理 + fem-rs submodule 更新（NCMesh/AMR/DenseTensor）** | **—** | **✅** |
| v1.0.0 | 时域瞬态 (TD-FEM) + 生产就绪 | FDTD_PLAN.md | 🔲 |

---

## 快速参考：关键 fem-rs API

```rust
// 网格加载
let mesh = fem_io::read_msh_file(Path::new("mesh.msh"))?;

// FEM 空间
let space = fem_space::H1Space::new(&mesh, fem_element::ElementOrder::P1);

// 组装
let integrator = fem_assembly::DiffusionIntegrator::constant(1.0);
let K = fem_assembly::Assembler::assemble_bilinear(&space, &[&integrator]);
let f = fem_assembly::Assembler::assemble_linear(&space, &[&source_integrator]);

// 边界条件
fem_space::apply_dirichlet(&mut K, &mut f, &dof_ids, 0.0);

// 求解
let config = fem_solver::SolverConfig { rtol: 1e-8, max_iter: 1000, ..Default::default() };
let solution = fem_solver::solve_pcg_amg(K, f, &config)?;

// 输出
let writer = fem_io::VtkWriter::new(Path::new("output.vtu"))?;
writer.write_mesh(&mesh)?;
writer.add_nodal_scalar("phi", &solution)?;
writer.finish()?;
```

---

### 常见问题与排查 (Web/WASM)

- **TypeError: jsmpi.Init is not a function**: 确保 `worker.js` 中 `self.jsmpi` 完整实现了 `Init`, `Finalize`, `Comm_size`, `Comm_rank` 等桩函数。由于 `wasm-bindgen` 的命名空间绑定，这些函数必须存在于 JavaScript 对象中。
- **文件路径无法访问**: WASM 运行在浏览器沙箱中，无法直接读取 `Model.Mesh` 指定的本地路径。建议通过 JavaScript 层 `fetch` 或 `FileReader` 获取字节流后传入 WASM 接口。

*最后更新: 2026-04-05 | 作者: Claude (AI Agent 开发指南)*

---

## 阶段 9: MoM/BEM 基础设施（v0.4.0）[COMPLETED]

**目标**: 建立 `crates/mom` crate 骨架，实现密集复数矩阵、表面网格提取、高斯求积规则，不含任何 Green 函数逻辑。

**前置条件**: 阶段 1-8 完成（FEM 求解器可用，Palace 兼容测试通过）。

### 9.1 创建 crates/mom

在根 `Cargo.toml` 的 `workspace.members` 中追加 `"crates/mom"`，并创建如下结构：

```
crates/mom/
├── Cargo.toml
└── src/
    ├── lib.rs
    ├── surface_mesh.rs   # 表面网格提取
    ├── quadrature.rs     # 高斯求积规则
    ├── matrix.rs         # 密集复数矩阵（Z 矩阵）
    └── config.rs         # MoM 配置节解析
```

**Cargo.toml 核心依赖**:

```toml
[dependencies]
rem-core     = { workspace = true }
rem-config   = { workspace = true }
rem-mesh     = { workspace = true }
num-complex  = { workspace = true }
faer         = { version = "0.19", features = ["complex"] }   # 密集矩阵 + LU
rayon        = "1"                                             # 并行装配行
```

### 9.2 配置解析（config.rs）

扩展 `crates/config/src/schema.rs`：

```rust
// 在 ProblemType 中新增
#[derive(Debug, Clone, Deserialize, PartialEq, Eq)]
pub enum ProblemType {
    Electrostatic,
    Magnetostatic,
    Eigenmode,
    Driven,
    Transient,
    MoM,   // 新增
    BEM,   // 新增
}

// 在 SolverConfig 中新增
#[derive(Debug, Clone, Deserialize, Default)]
pub struct SolverConfig {
    // ...existing fields...
    #[serde(rename = "MoM", default)]
    pub mom: Option<MomSolverConfig>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct MomSolverConfig {
    #[serde(rename = "Equation", default = "default_mom_equation")]
    pub equation: String,       // "EFIE" | "MFIE" | "CFIE" | "PMCHWT"

    #[serde(rename = "Basis", default = "default_mom_basis")]
    pub basis: String,          // "RWG" | "Pulse"

    #[serde(rename = "FreqMin")]
    pub freq_min: f64,

    #[serde(rename = "FreqMax")]
    pub freq_max: f64,

    #[serde(rename = "FreqStep")]
    pub freq_step: f64,

    #[serde(rename = "Alpha", default = "default_cfie_alpha")]
    pub alpha: f64,             // CFIE 混合系数，0=EFIE，1=MFIE

    #[serde(rename = "SingularTol", default = "default_singular_tol")]
    pub singular_tol: f64,

    #[serde(rename = "FastSolver", default = "default_fast_solver")]
    pub fast_solver: String,    // "Direct" | "ACA" | "FMM"
}

fn default_mom_equation() -> String { "CFIE".to_string() }
fn default_mom_basis()     -> String { "RWG".to_string()  }
fn default_cfie_alpha()    -> f64    { 0.5 }
fn default_singular_tol()  -> f64    { 1e-6 }
fn default_fast_solver()   -> String { "Direct".to_string() }
```

同时在 `PalaceConfig` 顶层新增 `Postprocessing` 节：

```rust
#[derive(Debug, Clone, Deserialize, Default)]
pub struct PalaceConfig {
    // ...existing fields...
    #[serde(rename = "Postprocessing", default)]
    pub postprocessing: Postprocessing,
}

#[derive(Debug, Clone, Deserialize, Default)]
pub struct Postprocessing {
    #[serde(rename = "RCS", default)]
    pub rcs: Option<RcsConfig>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct RcsConfig {
    #[serde(rename = "PhiDeg", default)]
    pub phi_deg: Vec<f64>,

    /// 支持 "0:5:180" 范围字符串
    #[serde(rename = "ThetaDeg", deserialize_with = "deserialize_angle_range")]
    pub theta_deg: Vec<f64>,
}
```

### 9.3 表面网格提取（surface_mesh.rs）

从 `RemMesh.boundary_elements` 提取三角面元，并构建 RWG 所需的边-面拓扑：

```rust
/// 单个三角面元
#[derive(Debug, Clone)]
pub struct TriFace {
    pub nodes: [usize; 3],          // 全局节点索引
    pub centroid: [f64; 3],
    pub normal: [f64; 3],           // 单位外法向
    pub area: f64,
}

/// 共享边（RWG 基函数载体）
#[derive(Debug, Clone)]
pub struct SharedEdge {
    pub nodes: [usize; 2],          // 边的两个节点（升序）
    pub plus_face: usize,           // T+ 面索引
    pub minus_face: usize,          // T- 面索引
    pub length: f64,
}

pub struct SurfaceMesh {
    pub nodes: Vec<[f64; 3]>,       // 节点坐标
    pub faces: Vec<TriFace>,
    pub edges: Vec<SharedEdge>,     // 内部共享边（RWG 基函数）
    pub boundary_edges: Vec<[usize; 2]>, // 边界边（仅属于一个面）
}

impl SurfaceMesh {
    /// 从 RemMesh 提取具有指定物理组属性的三角面
    pub fn extract(rem_mesh: &RemMesh, pec_attrs: &[u32]) -> Self { ... }

    /// 构建边-面拓扑（O(E log E)，排序后 hash）
    fn build_edge_topology(faces: &[TriFace]) -> (Vec<SharedEdge>, Vec<[usize; 2]>) { ... }
}
```

### 9.4 高斯求积规则（quadrature.rs）

提供三角形面元上的数值积分规则：

```rust
/// 三角形上的 Dunavant 求积点（n 阶精确）
pub struct TriQuad {
    pub points: Vec<[f64; 3]>,   // 重心坐标 (ξ₁, ξ₂, ξ₃)
    pub weights: Vec<f64>,
}

impl TriQuad {
    pub fn new(order: usize) -> Self { ... }

    /// 将重心坐标映射到全局坐标
    pub fn global_point(tri: &TriFace, nodes: &[[f64; 3]], bary: &[f64; 3]) -> [f64; 3] { ... }
}

/// 在面元上对标量函数积分：∫_T f(x) dS ≈ Σ wᵢ f(xᵢ) |J|
pub fn integrate_scalar<F: Fn(&[f64; 3]) -> f64>(
    face: &TriFace,
    nodes: &[[f64; 3]],
    quad: &TriQuad,
    f: F,
) -> f64 { ... }
```

支持的阶次与求积点数：

| 阶次 | 求积点数 | 精确多项式阶 | 用途 |
|------|---------|------------|------|
| 1 | 1 | 1 | 粗略估计 |
| 3 | 4 | 3 | 远场常规积分 |
| 5 | 7 | 5 | 默认标准积分 |
| 7 | 13 | 7 | 高精度积分 |
| 9 | 19 | 9 | 近奇异积分外层 |

**验收标准**:
- [ ] `crates/mom` 加入 workspace，`cargo build -p rem-mom` 编译通过
- [ ] 新增 `ProblemType::MoM`/`BEM`，配置解析测试通过（`cargo test -p rem-config`）
- [ ] 从 sphere.msh 提取表面网格，面元数、边数、节点数正确
- [ ] 7 点高斯规则对 2 次多项式积分误差 < 1e-14

---

## 阶段 10: EFIE 求解器（脉冲基函数）（v0.5.0）[COMPLETED]

**目标**: 实现最简单的 EFIE 求解器（脉冲/常数基函数，标量版），验证端到端流程，建立可对比的参考解。

### 10.1 Green 函数（green.rs）

```rust
use num_complex::Complex64;
use std::f64::consts::PI;

/// 3D 自由空间标量 Green 函数
/// G(r, r') = exp(-jkR) / (4πR)，R = |r - r'|
pub fn green3d(r: &[f64; 3], r_prime: &[f64; 3], k: f64) -> Complex64 {
    let rx = r[0] - r_prime[0];
    let ry = r[1] - r_prime[1];
    let rz = r[2] - r_prime[2];
    let dist = (rx*rx + ry*ry + rz*rz).sqrt();
    if dist < 1e-14 { return Complex64::ZERO; }  // 奇异点由专用积分处理
    let phase = Complex64::new(0.0, -k * dist);
    phase.exp() / (4.0 * PI * dist)
}

/// ∂G/∂n' = G(jkR + 1)/R² * (r-r')·n'
pub fn green3d_normal_deriv(
    r: &[f64; 3], r_prime: &[f64; 3], n_prime: &[f64; 3], k: f64
) -> Complex64 { ... }
```

### 10.2 阻抗矩阵装配（pulse 基函数版）

脉冲基函数最简单：每个面元一个未知量，基函数 = 1 on Tₙ, 0 otherwise。

```rust
/// 装配 N×N 阻抗矩阵（脉冲基函数 EFIE）
/// Z[m,n] = -jωμ₀ ∫_Tm ∫_Tn G(r,r') dS' dS  （远场块）
///         + 奇异修正项                         （对角块 m==n）
pub fn assemble_efie_pulse(
    surf: &SurfaceMesh,
    freq: f64,
    quad: &TriQuad,
    singular_tol: f64,
) -> faer::Mat<Complex64> {
    let n = surf.faces.len();
    let mut z = faer::Mat::<Complex64>::zeros(n, n);
    let k = 2.0 * PI * freq / C0;
    let omega = 2.0 * PI * freq;

    // 外层循环可并行（rayon）
    z.par_col_chunks_mut(1).enumerate().for_each(|(n_idx, mut col)| {
        for m_idx in 0..n {
            col[m_idx] = zmn_pulse(&surf.faces[m_idx], &surf.faces[n_idx], &surf.nodes, k, omega, quad);
        }
    });
    z
}

fn zmn_pulse(...) -> Complex64 {
    if m_idx == n_idx {
        zmn_singular(...)   // Duffy 变换
    } else {
        zmn_regular(...)    // 标准高斯积分
    }
}
```

### 10.3 Duffy 变换（奇异自积分）

```
对角块 Z[m,m]：将三角形 T 分成 3 个子三角形（以源点为顶点），
Duffy 变换消去 1/R 奇异性，转化为规则积分。

参考：Rao, Wilton, Glisson (1982)，附录 B
```

```rust
/// 计算自积分 ∫_T ∫_T G(r,r') dS' dS（Duffy 变换）
fn zmn_singular(face: &TriFace, nodes: &[[f64; 3]], k: f64, omega: f64) -> Complex64 {
    // 1. 将 T 分成 3 个以 r₁, r₂, r₃ 为源点极的子三角形
    // 2. 每个子三角形做 Duffy 变换: u = ρ cos θ, v = ρ sin θ（极坐标）
    // 3. 分母 R 被 Jacobian ρ 约消，得到光滑被积函数
    // 4. 对 (ρ, θ) 做标准高斯积分
    ...
}
```

### 10.4 主求解流程

```
SurfaceMesh::extract(rem_mesh, pec_attrs)
  └─ assemble_efie_pulse(surf, freq, quad)   → Z (N×N 复数密集)
       └─ 激励向量 V (入射平面波在面元中心)
            └─ faer::LU::factorize(Z)
                 └─ LU.solve(V)              → I (面电流密度)
                      └─ compute_rcs(I, surf, freq, theta, phi)
                           └─ write_rcs_csv(output_dir, ...)
```

### 10.5 Palace 配置集成

在 `crates/cli/src/main.rs` 的 `match` 分支中追加：

```rust
ProblemType::MoM => rem_mom::run(&config)?,
ProblemType::BEM => rem_bem::run(&config)?,
```

**验收标准**:
- [ ] PEC 球体（半径 0.1λ）EFIE 脉冲基函数 RCS 与 Mie 解析解误差 < 15%（脉冲基函数精度有限，预期误差较大）
- [ ] Palace 配置格式解析：`Problem.Type = "MoM"` 正确路由到 MoM 求解器
- [ ] `PEC.Attributes` 在 MoM 中正确识别为导体面
- [ ] RCS CSV 输出格式与 `Postprocessing.RCS` 配置一致
- [ ] `cargo test -p rem-mom` 单元测试全部通过（Green 函数、求积规则、Duffy 变换）

---

## 阶段 11: RWG 基函数 + CFIE（v0.6.0）[COMPLETED]

**目标**: 实现 RWG 矢量基函数和 CFIE 方程，达到生产精度（PEC 球体 RCS 误差 < 5%）。

### 11.1 RWG 基函数实现

RWG 基函数定义在共享边 eₙ 上，载体为 T⁺ 和 T⁻ 两个三角形：

```
f_n(r) = {  lₙ/(2A⁺) * (r - r⁺_free)   if r ∈ T⁺
          { -lₙ/(2A⁻) * (r - r⁻_free)   if r ∈ T⁻
          {  0                            otherwise
```

其中 lₙ 为边长，A± 为面积，r±_free 为 T± 中不属于该边的顶点。

```rust
/// RWG 基函数（定义在 SharedEdge 上）
pub struct RwgBasis {
    pub edge_idx: usize,
    pub plus_face: usize,
    pub minus_face: usize,
    pub free_node_plus: usize,   // T+ 中的自由顶点
    pub free_node_minus: usize,  // T- 中的自由顶点
    pub length: f64,             // 边长 lₙ
}

impl RwgBasis {
    /// 在 T± 上对 f_n 求值
    pub fn eval(&self, r: &[f64; 3], surf: &SurfaceMesh, in_plus: bool) -> [f64; 3] { ... }

    /// ∇·f_n = ±lₙ/A± （常数）
    pub fn divergence(&self, surf: &SurfaceMesh, in_plus: bool) -> f64 { ... }
}
```

### 11.2 CFIE 阻抗矩阵

CFIE = α·EFIE + (1-α)·η₀·MFIE，α ∈ [0,1]。

**EFIE 项**（矢量 RWG）:
```
Z_EFIE[m,n] = -jωμ₀ ∫_Tm f_m(r)·∫_Tn G(r,r') f_n(r') dS' dS
             + j/(ωε₀) ∫_Tm ∇_s·f_m(r) ∫_Tn G(r,r') ∇'_s·f_n(r') dS' dS
```

**MFIE 项**（仅适用外问题）:
```
Z_MFIE[m,n] = δ_{mn}/2 * ∫_Tm f_m dS
             + ∫_Tm f_m(r)·n̂×∫_Tn ∇G×f_n dS' dS
```

### 11.3 Sauter-Schwab 奇异积分

需处理四种几何接触情形（比脉冲基函数更复杂）：

| 情形 | 条件 | 处理方法 |
|------|------|---------|
| 完全重合 (4D) | T = T' | Duffy 4D 变换（6 个子块） |
| 共享边 (3D) | T ∩ T' = 一条边 | Sauter-Schwab Rule 2 |
| 共享顶点 (2D) | T ∩ T' = 一个顶点 | Sauter-Schwab Rule 3 |
| 非接触 (0D) | T ∩ T' = ∅ | 标准 7 点高斯 |

**实现顺序**：先实现非接触 + 完全重合（最简单的两种），再逐步加入共享边和共享顶点。

### 11.4 RCS 后处理

```rust
/// 计算双站 RCS（雷达截面积）
/// σ(θ,φ) = 4π|F(θ,φ)|²/|E_inc|²
/// F = ∫_S J_s(r') × exp(jk r̂·r') dS'
pub fn compute_rcs(
    currents: &[Complex64],  // RWG 展开系数
    surf: &SurfaceMesh,
    freq: f64,
    theta_deg: &[f64],
    phi_deg: &[f64],
    output_dir: &Path,
) -> RemResult<()> { ... }
```

**输出文件**（Palace 扩展格式）：
```
{output_dir}/postpro/rcs.csv
Freq (GHz),Theta (deg),Phi (deg),RCS (dBsm)
1.0,0,0,-15.3
1.0,10,0,-14.8
...
```

**验收标准**:
- [ ] PEC 球体（半径 0.1λ）CFIE RCS 与 Mie 解析解误差 < 5%
- [ ] α=1（纯 EFIE）和 α=0（纯 MFIE）分别可独立运行
- [ ] Sauter-Schwab 4 种情形单元测试通过
- [ ] RCS CSV 输出可与 `postprocessing.RCS` 配置正确对应

---

## 阶段 12: Laplace BEM 静电求解器（v0.7.0）[COMPLETED]

**目标**: 实现 Laplace BEM，与现有 FEM 静电求解器交叉验证，同时扩展 `Problem.Type = "BEM"` 路由。

### 12.1 BEM 方程

对 Laplace 方程 ∇²φ = 0，外 Dirichlet 问题（PEC 表面）：

```
双层势表示:  φ(r) = ∫_S [φ(r') ∂G_L/∂n' - G_L(r,r') q(r')] dS'
           G_L(r,r') = 1/(4π|r-r'|)
```

离散化后为 2×2 块系统（H·φ = G·q），P1 节点基函数足以达到二阶精度。

### 12.2 与 Palace 配置的映射

- `Boundaries.PEC.Attributes` → 导体表面（φ = 已知值）
- `Boundaries.Terminal` → 激励端口（φ = Vₙ）
- `Boundaries.Ground` → 接地（φ = 0）
- 输出电容矩阵格式与 FEM 静电求解器完全相同（`capacitance.csv`）

### 12.3 验证基准

| 算例 | FEM 参考值 | BEM 目标误差 |
|------|-----------|------------|
| 平行板电容 | ε₀A/d | < 1% |
| 球体电容 | 4πε₀R | < 0.5% |
| 同轴线电容/长度 | 2πε₀/ln(b/a) | < 0.5% |

**验收标准**:
- [ ] 球体电容与解析解误差 < 0.5%
- [x] 与 FEM 静电结果（同网格）交叉误差 < 1%
- [x] `Problem.Type = "BEM"` 完整路由，Palace 其余配置字段正确传递

---

## MoM/BEM 开发注意事项

### E. Palace 兼容性守护规则

1. **不修改现有 Palace 字段语义**: `ProblemType` 已有值的处理逻辑绝对不能因 MoM/BEM 实现而改变。新增枚举值只做追加，不做替换。
2. **新增字段必须有默认值**: `Solver.MoM`、`Postprocessing` 等新节必须加 `#[serde(default)]`，保证 Palace 配置在 REM 中无需修改即可运行。
3. **兼容性回归测试不可删除**: `tests/integration/test_config_compat.rs` 中的 Palace 示例测试在每次 PR 中必须通过。
4. **MoM 输出目录与 FEM 并列**: 不覆盖 FEM 输出，`{output_dir}/postpro/rcs.csv` vs `{output_dir}/postpro/domain-E.csv`。

### F. MoM 数值稳定性

1. **奇异积分是最高风险点**: 优先实现 Duffy 自积分，其次是 Sauter-Schwab，最后是近场近奇异。每种情形单独单元测试，对 Green 函数的已知积分解析解比对。
2. **CFIE α 参数选择**: 默认 α=0.5（经典 CFIE）。纯 EFIE（α=1）在内谐振频率附近会发散——在文档和日志中明确标注。
3. **条件数监控**: 装配完 Z 矩阵后输出估计条件数（`log::info!`），若 > 1e12 发出 `log::warn!`。
4. **WASM 限制**: `faer` 支持 WASM，但 Z 矩阵 LU 分解在 WASM 下仅单线程（`rayon` WASM 支持有限）。建议 WASM 模式限制 N < 1000，超出则提示用 native 模式。

---

## 阶段 13: SBR+ 高频射线追踪 + PO（v0.8.0）[COMPLETED]

**目标**: 实现 SBR+（Shooting and Bouncing Rays Plus）高频渐近散射求解器，与 Mie 级数解析解对比验证。

**参考文档**: [SBR_PLUS_PLAN.md](SBR_PLUS_PLAN.md)

### 13.1 核心实现

| 模块 | 文件 | 说明 |
|------|------|------|
| BVH 加速结构 | `crates/sbr/src/bvh.rs` | AABB BVH，SAH 分割，Möller-Trumbore 三角形求交 |
| 射线数据结构 | `crates/sbr/src/ray.rs` | Ray / RayHit / RayPath |
| 平面波激励 | `crates/sbr/src/excitation.rs` | 孔径铺设 + 平面波 E/H 场 |
| Fresnel 系数 | `crates/sbr/src/fresnel.rs` | PEC 镜面反射 + PO 感应电流 `J = 2n̂×H` |
| 远场 PO 积分 | `crates/sbr/src/po_integral.rs` | N(r̂) = Σ J_m A_m exp(jkr̂·r_m) → RCS |
| 输出 | `crates/sbr/src/output.rs` | RCS CSV + 感应电流 VTK |
| 主流程 | `crates/sbr/src/lib.rs` | 两阶段算法：first_bounce_po + multibounce_rays |

### 13.2 关键设计决策

**两阶段 PO 算法**：
- **一次弹射**（`first_bounce_po`）：逐面迭代（per-face），与射线密度无关，阴影测试从 `face.centroid + ε·face.normal` 发出
- **多次弹射**（`multibounce_rays`）：射线管追踪，`J` 贡献乘以 `A_ray/A_face` 比例系数

**mesh 分辨率约束**：PO 相位积分要求面片尺寸 < λ/4。验证测试用 1 GHz（ka=10.5），24 纬度环充分。

### 13.3 验证结果

| 测试 | 参考解 | 结果 |
|------|--------|------|
| PEC 球（r=0.5m）@ 1 GHz，ka=10.5，单站 RCS | Mie 级数 | 误差 0.05 dB（< 3 dB 限值）✅ |

---

## 阶段 13.5: 警告清理 + fem-rs 更新（v0.8.1）[COMPLETED]

### 工作内容

1. **fem-rs submodule 更新**（`vendor/fem-rs`）：bc55578 → bf0ab3a
   - 新增：NCMesh 悬挂节点支持
   - 新增：DenseTensor / LU 分解（`dense.rs`）
   - 新增：ZZ/Kelly/Dörfler AMR 误差估计器
   - 新增：Quad4/Hex8 等参元（完整积分公式）

2. **Rust 编译警告清理**：
   - `crates/bem/src/kernel.rs`：`#[allow(non_snake_case)]` 保留数学命名（`laplace_G`, `laplace_dG_dn`）
   - `crates/bem/src/assemble.rs`：移除死代码 `type C64 = faer::c64`
   - `crates/driven/src/lib.rs` + `output.rs`：`FreqResult` / `write_s_params` → `pub(crate)`
   - `crates/eigenmode/src/lib.rs`：移除 `BoundaryTag` 无用导入；`tol` → `_tol`
   - `crates/electrostatic/src/lib.rs` + `crates/magnetostatic/src/lib.rs`：`NoComm` 移入 `#[cfg(test)]` 模块
   - `crates/mesh/src/gmsh.rs`：`phys_names` → `_phys_names`（仅用于数值 tag，不用字符串名）

---

*最后更新: 2026-04-05 | v0.8.1 完成*

---

## 阶段 14: Palace v0.16 差距补全 P11–P13（v0.15.0）[COMPLETED]

**目标**: 补全 Palace v0.16 相对 REM v0.14 的关键差距，并新增 Driven 求解器扩展能力。

### 完成项清单

| 任务 | 关键文件 | 说明 |
|------|---------|------|
| **P11-1** 完整 N×N S 参数矩阵 | `crates/driven/src/lib.rs` | `FreqResult.s_matrix: Vec<Vec<Complex64>>`；多端口 Z→S 转换；WASM flat 格式 |
| **P11-2** 材料各向异性 ε/μ 张量装配 | `crates/materials/src/material.rs`, `crates/electrostatic/src/assemble.rs` | `epsilon_tensor: [[f64;3];3]`；`MaterialAxes` 旋转矩阵；`assemble_stiffness_aniso()`；静电/特征模/驱动三路均已接入 |
| **P11-3** 导体 Q 因子（R_s 表面欧姆损耗） | `crates/eigenmode/src/lib.rs` | R_s = √(ωμ₀/2σ)；微扰法 1/Q_c；与介质 tan δ Q_d 合并；`Boundaries.Conductivity` 触发 |
| **P11-4** 电流偶极子点源激励 | `crates/driven/src/lib.rs` | Palace v0.16 `Domains.CurrentDipole`；Hertz 偶极子 jω μ₀ Il 注入最近自由节点 |
| **P12-1** Floquet 周期边界条件（Γ 点） | `crates/core/src/sparse.rs`, `crates/electrostatic/src/bc.rs`, `crates/eigenmode/src/lib.rs` | `TripletMatrix::remap_periodic_nodes()`；`collect_periodic_node_pairs()` 几何平移匹配；非零 FloquetWaveVector → 警告跳过 |
| **P12-2** Drude-Lorentz 频变材料 | `crates/materials/src/material.rs`, `crates/materials/src/domain_map.rs`, `crates/driven/src/lib.rs` | ε(ω) = ε∞ + Σ ωp²/(ω0²−ω²+jγω)；每频点修正刚度矩阵；扣除静态损耗重叠 |
| **P12-3** JSON Schema 运行时校验 | `crates/config/src/validate.rs`, `crates/config/src/lib.rs` | 两阶段：结构预校验 + 语义校验；5 个单元测试 |
| **P12-4** 内存峰值报告 | `crates/core/src/memory.rs` | Linux VmPeak；Windows GetProcessMemoryInfo (raw psapi)；WASM → None；各求解器完成时 log::info 输出 |
| **P12-5** 近远场变换（Kirchhoff 辐射方向图） | `crates/driven/src/far_field.rs`, `crates/wasm/src/lib.rs`, `crates/yew-app/src/solver.rs` | F(r̂)=∫E e^{jkr̂·r'}dS'；E=-∇φ 梯度恢复；dBi 归一化；`Solver.FarField` 配置；`far_field.csv` artifact |
| **P13-1** 快照 ROM 频率扫描加速 | `crates/driven/src/rom.rs`, `crates/config/src/schema.rs` | 修正 Gram-Schmidt 正交化快照基；`A_r(ω)=V†A(ω)V` r×r LU 求解；展开频率均匀/对数分布；`DrivenSolver.RomOrder` 配置；仅单端口可用；3 个单元测试 |

### 验收（构建验证）

```bash
cargo build -p rem-yew    # 0 错误
cargo test -p rem-driven rom    # 3 通过
cargo test -p rem-config        # 全部通过
```

*最后更新: 2026-04-09 | v0.15.0 完成*

---

## 阶段 15 — ROM 电路综合（P14-1）

| 编号 | 功能 | 状态 |
|------|------|------|
| P14-1 | Vector Fitting 电路综合 | ✅ |

### P14-1 Vector Fitting 电路综合

- 新文件 `crates/driven/src/vf.rs`：Gustavsen-Semlyen VF 算法（实数化 LS + Schur 特征值）
- Config: `DrivenSolver.CircuitSynthesis: bool`
- 输出：`s_params.s1p`（Touchstone）、`circuit_model.csv`（极点-留数）、`equivalent_circuit.cir`（SPICE）
- 11 项测试全通过，`cargo build -p rem-yew` 0 错误

验证：
```bash
cargo test -p rem-driven vf   # 4 tests ok
cargo build -p rem-yew         # 0 errors
```

*最后更新: 2026-04-09 | v0.16.0 完成*

---

## 阶段 16 — MoM 端口激励 + S 参数提取（v0.17.0）🔲

> 详细规格：[docs/MOM_Sonnet19_Alignment_Plan.md](docs/MOM_Sonnet19_Alignment_Plan.md) § 阶段 16

**目标**：在 MoM 框架内引入集总端口激励，输出 Touchstone `.sNp`，
使 REM 能仿真微带/CPW 无源器件的 S 参数（Sonnet 核心场景）。

### 关键任务

| 编号 | 功能 | 文件 | 状态 |
|------|------|------|------|
| P16-1 | `MomPort` 配置解析 + `RefImpedance` | `crates/config/src/schema.rs` | 🔲 |
| P16-2 | `crates/mom/src/port.rs`（集总端口模型） | 新文件 | 🔲 |
| P16-3 | `crates/mom/src/sparams.rs`（S 矩阵 + Touchstone） | 新文件 | 🔲 |
| P16-4 | 主流程端口分支（`lib.rs`） | `crates/mom/src/lib.rs` | 🔲 |
| P16-5 | 单/双端口验证测试（偶极子 + 传输线） | `crates/mom/tests/` | 🔲 |

### 验收标准

- [ ] `Ports` 配置解析正确，S 参数路由激活
- [ ] 双端口：2×2 S 矩阵写出 `.s2p`，ADS/Qucs 可读
- [ ] 半波偶极子 S11 vs 解析阻抗误差 < 5%
- [ ] 无端口时原有 RCS 路径零回归
- [ ] `cargo test -p rem-mom` 全部通过

---

## 阶段 17 — 分层介质 Green 函数（v0.18.0）🔲

> 详细规格：[docs/MOM_Sonnet19_Alignment_Plan.md](docs/MOM_Sonnet19_Alignment_Plan.md) § 阶段 17

**目标**：引入多层介质 Sommerfeld 积分 Green 函数（DCIM 离散复像法），
使 REM MoM 能仿真嵌入 PCB/MMIC 基板中的导体结构（Sonnet 最核心的物理能力）。

### 关键任务

| 编号 | 功能 | 文件 | 状态 |
|------|------|------|------|
| P17-1 | 新 crate `crates/layered_green`（传输矩阵 + DCIM） | 新 crate | 🔲 |
| P17-2 | `GreenFunction` trait 抽象（自由空间 / 分层介质） | `crates/mom/src/green_trait.rs` | 🔲 |
| P17-3 | `assemble_efie_rwg` / `assemble_mfie_rwg` 接受 trait | `crates/mom/src/assemble.rs` | 🔲 |
| P17-4 | `SubstrateConfig` / `LayerConfig` 配置解析 | `crates/config/src/schema.rs` | 🔲 |
| P17-5 | 贴片天线验证（FR4 基板） | `examples/patch_antenna/` | 🔲 |

### 验收标准

- [ ] DCIM 残差 vs 数值 Sommerfeld < 1e-3
- [ ] 贴片天线谐振频率误差 < 2%
- [ ] 原有自由空间 RCS 测试零回归

---

## 阶段 18 — 有损导体 SIBC + FFT 加速（v0.19.0）🔲

> 详细规格：[docs/MOM_Sonnet19_Alignment_Plan.md](docs/MOM_Sonnet19_Alignment_Plan.md) § 阶段 18

**目标**：引入表面阻抗边界条件（SIBC）建模铜导体损耗；
在严格平面结构上实现 FFT 加速矩阵-向量积（O(N log N)）。

### 关键任务

| 编号 | 功能 | 文件 | 状态 |
|------|------|------|------|
| P18-1 | `crates/mom/src/sibc.rs`（SIBC 修正阻抗矩阵） | 新文件 | 🔲 |
| P18-2 | `ConductivityBc` 扩展（MoM 路径） | `crates/config/src/schema.rs` | 🔲 |
| P18-3 | `crates/mom/src/fft_accel.rs`（FFT 加速求解器） | 新文件 | 🔲 |
| P18-4 | 铜微带传输线有损验证 | `examples/microstrip_lossy/` | 🔲 |

### 验收标准

- [ ] 趋肤深度公式验证 @ 1-30 GHz
- [ ] 铜微带 S21 vs Sonnet 参考 ΔS21 < 0.1 dB
- [ ] FFT 加速：N=2000 时速度提升 > 5×，精度 < 0.1 dB

---

## 阶段 19 — MoM AMR + 频率 ROM + Touchstone 完整（v0.20.0）🔲

> 详细规格：[docs/MOM_Sonnet19_Alignment_Plan.md](docs/MOM_Sonnet19_Alignment_Plan.md) § 阶段 19

**目标**：补齐剩余工程可用性差距：自适应网格细化、频率扫描 ROM 加速、
完整 Touchstone 2.0 格式兼容。

### 关键任务

| 编号 | 功能 | 文件 | 状态 |
|------|------|------|------|
| P19-1 | `crates/mom/src/amr.rs`（表面电流梯度误差指示 + Dörfler 标记 + 面元细化） | 新文件 | 🔲 |
| P19-2 | MoM 频率扫描快照 ROM（`RomOrder` 配置） | `crates/mom/src/rom.rs` | 🔲 |
| P19-3 | Touchstone 2.0 完整（MA/RI/DB 格式、注释行、多端口规范） | `crates/mom/src/sparams.rs` | 🔲 |
| P19-4 | `FEATURE_COMPARISON.md` 全面更新 | `FEATURE_COMPARISON.md` | 🔲 |

### 验收标准

- [ ] AMR：PEC 球 3 次迭代收敛（RCS 变化 < 0.1 dB）
- [ ] ROM：100 频率点 = 10 锚点精度，S 参数误差 < 0.05 dB
- [ ] `.s2p` 可被 ADS、Qucs、scikit-rf 读取
- [ ] `cargo test --workspace` 全部通过

---

## 版本里程碑（更新）

| 版本 | 内容 | 对应阶段 | 状态 |
|------|------|---------|------|
| v0.1.0 | 工作区 + 配置解析 + 静电/静磁 + CLI + VTK 输出 | 1-5 | ✅ |
| v0.1.1 | WASM 绑定 + Yew Web Demo（可部署） | 6, 7.5 | ✅ |
| v0.2.0 | 并行层 (jsmpi/Comm trait) + 分布式组装 | 6.5, 7 | ✅ |
| v0.2.1 | rmetis 子模块 + METIS k-way 对偶图分区 | — | ✅ |
| v0.3.0 | Palace 官方示例完整兼容测试 | 8 | ✅ |
| v0.4.0 | MoM 基础设施：密集复数矩阵 + 表面网格提取 + 高斯求积 | 9 | ✅ |
| v0.5.0 | MoM EFIE 求解器：Green 函数 + 脉冲基函数 + 验证 | 10 | ✅ |
| v0.6.0 | MoM CFIE 求解器：RWG 基函数 + 奇异积分 + RCS 输出 | 11 | ✅ |
| v0.7.0 | BEM 静态求解器：Laplace BEM + 与 FEM 交叉验证 | 12 | ✅ |
| v0.8.0 | SBR+ 高频 PO 求解器：BVH + 两阶段 PO + Mie 验证 | — | ✅ |
| v0.8.1 | 警告清理 + fem-rs submodule 更新 | — | ✅ |
| v0.14.0–v0.16.0 | Palace v0.16 差距补全（材料、端口、ROM、电路综合） | 14-15 | ✅ |
| **v0.17.0** | **MoM 集总端口 + S 参数 + Touchstone 基础** | **16** | **🔲** |
| **v0.18.0** | **分层介质 Green 函数（DCIM）+ 基板配置** | **17** | **🔲** |
| **v0.19.0** | **SIBC 有损导体 + FFT 加速平面 MoM** | **18** | **🔲** |
| **v0.20.0** | **MoM AMR + 频率扫描 ROM + Touchstone 完整** | **19** | **🔲** |
| v1.0.0 | 时域瞬态 (TD-FEM) + 生产就绪 | FDTD_PLAN.md | 🔲 |

*最后更新: 2026-04-10 | 新增 Sonnet 19 对齐路线图（v0.17–v0.20）*
