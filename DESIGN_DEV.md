# REM — 设计与开发文档
## AI Agent 工作指南 v0.1

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
        ProblemType::Electrostatic  => rem_electrostatic::run(&config)?,
        ProblemType::Magnetostatic  => rem_magnetostatic::run(&config)?,
        ProblemType::Eigenmode      => {
            eprintln!("Eigenmode solver not yet implemented (v0.2)");
            std::process::exit(1);
        }
        ProblemType::Driven         => {
            eprintln!("Driven solver not yet implemented (v0.2)");
            std::process::exit(1);
        }
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

# 使用 wasm-pack（输出到 web/src/pkg，供 Vue3 直接 import）
wasm-pack build crates/wasm --target web --out-dir ../../web/src/pkg

# 测试 WASM 功能（在 node 中运行）
wasm-pack test crates/wasm --node
```

**验收标准**:
- [ ] `cargo build --target wasm32-unknown-unknown -p rem-wasm` 无错
- [ ] 从 JavaScript 调用 `RemSolver.new()`, `load_config()`, `load_mesh()`, `solve()` 返回正确结果
- [ ] WASM bundle 大小 < 5 MB（--release + opt-level=s）

---

## 阶段 7: 特征模与频域求解器（v0.2）

> **前置条件**: 阶段 1-5 完成，fem-rs Phase 5（Nedelec 元）已合并

### 7.1 特征模求解器

**方程**: `curl(μᵣ⁻¹ curl **E**) = k₀² εᵣ **E**`

**离散化**:
1. 构建 curl-curl 矩阵 A（刚度）和质量矩阵 B（质量）
2. 解广义特征值问题: `A x = λ B x`，λ = k₀² = (ω/c)²

**特征值求解器**:
- 纯 Rust 实现: LOBPCG（局部最优块预条件共轭梯度）
- 外部接口: 预留 `slepc-rs` / `arpack` FFI 接口
- 默认返回最小 N 个正实特征值

**输出**:
```
postpro/eig.csv:
Mode,Freq (GHz),Q Factor,Error
1,5.123456,1234.5,1.23e-6
```

### 7.2 频域驱动求解器

**方程**: `[A - k₀² B + jωC] E = f`

其中 C 包含导体损耗和集总端口阻抗。

**频率扫描**:
```rust
for freq in linspace(config.solver.driven.min_freq, max_freq, n_steps) {
    let k0 = 2.0 * PI * freq / C0;
    let system = assemble_driven_system(k0, &mesh, &config);
    let solution = solve_complex_system(system)?;
    let s_params = compute_s_params(&solution, &config.boundaries.lumped_port);
    write_s_params(freq, &s_params, output_dir)?;
}
```

**S 参数计算**（集总端口）:
```
S_ij = 2√(R_j/R_i) · V_i / V_j^inc - δ_ij
```

**验收标准** (v0.2):
- [ ] 圆柱谐振腔特征频率与解析解误差 < 0.1%（TM₀₁₀ 模式）
- [ ] 简单传输线 S₂₁ 在通带内 > -3 dB
- [ ] S 参数 CSV 输出与 Palace 格式兼容

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

## 阶段 7.5: Web Demo（web/）

**目标**: 构建 Vue3 静态 Demo 页面，集成 WASM 求解器，可零配置部署至 GitHub Pages。

> **前置条件**: 阶段 6（WASM 绑定）完成，`web/src/pkg/` 下有有效的 wasm-pack 产物。

### 7.5.1 初始化项目

```bash
cd web
npm create vue@latest . -- --typescript --router=false --pinia --no-jsx
npm install naive-ui @vicons/ionicons5
npm install three @types/three
npm install @monaco-editor/vue
npm install pinia
```

**`vite.config.ts`** 关键配置（启用 WASM 支持）：

```typescript
import { defineConfig } from 'vite'
import vue from '@vitejs/plugin-vue'

export default defineConfig({
  plugins: [vue()],
  optimizeDeps: {
    exclude: ['rem-wasm'],  // 不预构建 WASM 包
  },
  server: {
    headers: {
      // Web Worker + SharedArrayBuffer 需要这两个头
      'Cross-Origin-Opener-Policy': 'same-origin',
      'Cross-Origin-Embedder-Policy': 'require-corp',
    },
  },
  build: {
    target: 'esnext',  // 支持顶层 await（WASM init 用）
  },
})
```

### 7.5.2 组件实现顺序

**Step 1 — ExampleSelector.vue**

侧边栏，从 `src/examples/*/meta.json` 动态加载示例列表。点击示例后：
1. 读取对应 `config.json` 内容 → 更新 Pinia `configStore.json`
2. 读取对应 `mesh.msh` 字节 → 更新 `meshStore.bytes`
3. 触发 MeshViewer 重新渲染

```typescript
// src/stores/example.ts
import { defineStore } from 'pinia'
import { ref } from 'vue'

export const useExampleStore = defineStore('example', () => {
  const configJson = ref('')
  const meshBytes = ref<Uint8Array | null>(null)
  const currentExample = ref('')

  async function loadExample(id: string) {
    const config = await fetch(`./examples/${id}/config.json`).then(r => r.text())
    const meshBuf = await fetch(`./examples/${id}/mesh.msh`).then(r => r.arrayBuffer())
    configJson.value = config
    meshBytes.value = new Uint8Array(meshBuf)
    currentExample.value = id
  }

  return { configJson, meshBytes, currentExample, loadExample }
})
```

**Step 2 — ConfigEditor.vue**

Monaco Editor 包装组件，绑定 `exampleStore.configJson`，提供 JSON 语法高亮。

```vue
<template>
  <MonacoEditor
    v-model:value="store.configJson"
    language="json"
    theme="vs-dark"
    :options="{ minimap: { enabled: false }, fontSize: 13 }"
    style="height: 100%; width: 100%"
  />
</template>
```

**Step 3 — SolverPanel.vue + Web Worker**

`SolverPanel` 支持两种模式，通过 UI 开关切换：
- **单 Worker 模式**（默认）：调用 `solver.worker.ts`，使用 `SerialComm`
- **多 Worker 并行模式**：调用 `mpi-coordinator.ts`，使用 `WorkerComm`（需浏览器支持 SharedArrayBuffer）

```typescript
// src/composables/useSolver.ts
export function useSolver() {
  const worker = new Worker(
    new URL('../worker/solver.worker.ts', import.meta.url),
    { type: 'module' }
  )
  const status = ref<'idle' | 'running' | 'done' | 'error'>('idle')
  const result = ref<SolveResult | null>(null)
  const log = ref<string[]>([])

  worker.onmessage = (e) => {
    if (e.data.type === 'ready') {
      status.value = 'idle'
    } else if (e.data.type === 'result') {
      result.value = e.data.payload
      status.value = 'done'
    } else if (e.data.type === 'error') {
      log.value.push(`ERROR: ${e.data.payload}`)
      status.value = 'error'
    }
  }

  async function solve(configJson: string, meshBytes: Uint8Array) {
    status.value = 'running'
    worker.postMessage({ type: 'solve', payload: { configJson, meshBytes: meshBytes.buffer } }, [meshBytes.buffer])
  }

  worker.postMessage({ type: 'init' })
  return { status, result, log, solve }
}
```

**Step 4 — MeshViewer.vue**

Three.js 渲染 GMSH 网格轮廓（2D 截面），从 `SolveResult.nodeCoords` + `connectivity` 构建 `THREE.BufferGeometry`。

**Step 5 — FieldViewer.vue**

渲染求解结果场量：
- 标量场（电位 φ）：节点颜色映射（jet colormap）
- 矢量场（E/B 场）：单元箭头（`THREE.ArrowHelper`）
- 色标条：CSS 渐变 + min/max 数值显示

**Step 6 — ResultTable.vue**

将 `SolveResult.csvOutputs` 中的 CSV 字符串解析后渲染为 Naive UI 的 `NDataTable`。

### 7.5.3 示例文件制作

每个预置示例网格必须满足：
- GMSH v4.1 ASCII 格式（`-format msh4`）
- 2D 示例文件 < 200 KB，3D 示例 < 500 KB
- 物理组 ID 与对应 `config.json` 中的 `Attributes` 完全一致

生成命令（需本地安装 gmsh）：
```bash
gmsh examples/meshes/parallel_plate.geo -2 -o web/src/examples/parallel_plate/mesh.msh -format msh4
gmsh examples/meshes/coaxial.geo -2 -o web/src/examples/coaxial/mesh.msh -format msh4
```

### 7.5.4 验收标准

- [ ] `npm run dev` 启动无控制台错误
- [ ] 点击"平行板"示例后，Monaco 编辑器显示对应 JSON，Three.js 渲染网格轮廓
- [ ] 点击"运行求解"后，Web Worker 执行，进度反馈显示，最终渲染电位场热图
- [ ] `npm run build` 产物 `dist/` 可用 `npx serve dist` 本地验证
- [ ] 所有预置示例（平行板、同轴线、方形导线）均可成功求解并显示结果
- [ ] 页面在 Chrome/Firefox/Edge 最新版均正常运行

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

| 版本 | 内容 | 对应阶段 |
|------|------|---------|
| v0.1.0 | 工作区 + 配置解析 + 静电/静磁 + CLI + VTK 输出 | 1-5 |
| v0.1.1 | WASM 绑定 + Vue3 Web Demo（可部署） | 6, 7.5 |
| v0.2.0 | 并行层 (jsmpi/Comm trait) + 分布式组装 | 6.5, 7 |
| v0.2.1 | rmetis 子模块 + METIS k-way 对偶图分区 | — |
| v0.3.0 | Palace 官方示例完整兼容测试 | 8 |
| v1.0.0 | 时域瞬态 + AMR | 未规划 |

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

*最后更新: 2026-04-01 | 作者: Claude (AI Agent 开发指南)*
