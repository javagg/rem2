# REM2 FDTD 求解器技术方案

> 版本：v1.0 草案  
> 日期：2026-04-05  
> 定位：时域全波电磁求解器，作为现有频域/静态求解器的时域补充

---

## 1. 背景与现状

### 1.1 项目现状

REM2 目前已具备以下求解能力：

| 求解器 | 物理方程 | 状态 |
|--------|---------|------|
| Electrostatic | −∇·(ε₀εᵣ ∇φ) = ρ | 完成 |
| Magnetostatic | −∇·(ν ∇A_z) = J_z | 完成 |
| Eigenmode | K x = λ M x | 完成 |
| Driven (频域) | −∇·(ε ∇φ) − k² ε φ = J | 完成 |
| Method of Moments | EFIE/MFIE/CFIE | v0.6 开发中 |
| **Transient (FDTD)** | **∇×∇×E + με ∂²E/∂t² = −μ ∂J/∂t** | **❌ 待实现** |

### 1.2 现有基础设施（可复用）

- `core`：稀疏矩阵 (CSR)、PCG 求解器、物理常数 (ε₀, μ₀, c₀)
- `mesh`：RemMesh 数据结构（支持 Tet4、Hex8 等）
- `materials`：材料属性映射（εᵣ, μᵣ, σ, tan δ）
- `result`：VTK/CSV 输出基础设施
- `parallel`：Comm trait（支持 NoComm 和 WASM/jsmpi）
- `config`：`ProblemType::Transient` 和 `TransientConfig` 骨架已定义

### 1.3 CLI 占位符

[cli/src/main.rs:62](crates/cli/src/main.rs#L62) 当前返回 `"Transient solver not yet implemented (v1.0)"`，待替换为实际调用。

---

## 2. 技术路线选型

### 2.1 时域有限元法（TD-FEM）vs. 传统 FDTD

本项目选择 **时域有限元法（TD-FEM / FETD）** 而非传统 Yee 格式 FDTD，理由如下：

| 特性 | 传统 FDTD (Yee) | TD-FEM (本方案) |
|------|----------------|----------------|
| 网格类型 | 结构化正交网格 | 任意非结构化网格（与现有 mesh 模块一致）|
| 曲面几何 | 阶梯近似，精度低 | 精确贴合曲面 |
| 与现有代码集成 | 需要重新建网格 | 直接复用 RemMesh |
| 时间步约束 | CFL 条件（显式）| 可用隐式方案（无条件稳定）|
| 吸收边界 | Mur/PML 成熟 | Mur/PML 均可实现 |
| 实现复杂度 | 简单，但与现有架构割裂 | 适中，与现有架构高度一致 |

> **结论**：采用 TD-FEM，基于 Nedelec H(curl) 矢量基函数，时间推进使用 Newmark-β 隐式方案。

### 2.2 控制方程

Maxwell 旋度方程（二阶矢量波动方程）：

```
μ ∂²E/∂t² + σ ∂E/∂t + ∇×(μ⁻¹ ∇×E) = −∂J_src/∂t
```

弱形式（乘以测试函数 v，在 Ω 上积分）：

```
∫ με ∂²E/∂t² · v dΩ
+ ∫ σ ∂E/∂t · v dΩ
+ ∫ μ⁻¹ (∇×E)·(∇×v) dΩ
= −∫ ∂J_src/∂t · v dΩ
```

对应矩阵方程：

```
M_ε ë + C_σ ė + K E = F(t)
```

其中：
- `M_ε` = ε 质量矩阵 `∫ ε E·v dΩ`
- `C_σ` = σ 阻尼矩阵 `∫ σ E·v dΩ`
- `K` = 刚度矩阵 `∫ μ⁻¹ (∇×E)·(∇×v) dΩ`
- `F(t)` = 激励源向量

### 2.3 时间推进方案：Newmark-β

采用无条件稳定的 Newmark-β 方案（β=0.25, γ=0.5，即梯形法则）：

```
E_{n+1} = E_n + Δt Ė_n + Δt²[(0.5−β)Ë_n + β Ë_{n+1}]
Ė_{n+1} = Ė_n + Δt[(1−γ)Ë_n + γ Ë_{n+1}]
```

每时间步求解的线性系统：

```
A_eff · E_{n+1} = F_eff
A_eff = M_ε/(βΔt²) + γC_σ/(βΔt) + K
```

> 优势：无条件稳定，时间步长仅受精度约束，不受 CFL 条件限制。  
> 替代方案：如追求效率可使用显式中心差分（γ=0.5, β=0），但需满足 CFL。

---

## 3. 核心组件设计

### 3.1 Nedelec 矢量基函数（H(curl)，P1）

这是 FDTD 实现的**关键路径**，现有代码中不存在。

**DOF 分配**：每条棱边（edge）分配 1 个自由度（P1 Nedelec）

对于四面体单元 $K$ 的第 $i$ 条棱 $e_i = (a_i, b_i)$，基函数为：

```
φ_i = l_i (λ_{a_i} ∇λ_{b_i} − λ_{b_i} ∇λ_{a_i})
```

其中 `λ_j` 为第 `j` 个节点的重心坐标，`l_i` 为棱长。

**局部旋度**：

```
∇×φ_i = 2 l_i ∇λ_{a_i} × ∇λ_{b_i}
```

**文件位置**：新建 `crates/transient/src/nedelec.rs`

```rust
pub struct NedelecBasis {
    pub dof_map: Vec<(usize, usize)>,  // edge → (node_a, node_b)
    pub edge_to_dof: HashMap<[usize;2], usize>,
    pub n_dofs: usize,
}

impl NedelecBasis {
    pub fn from_mesh(mesh: &RemMesh) -> Self { ... }
    pub fn phi(edge_idx: usize, lambda: &[f64;4], grad_lambda: &[[f64;3];4]) -> [f64;3] { ... }
    pub fn curl_phi(edge_idx: usize, grad_lambda: &[[f64;3];4]) -> [f64;3] { ... }
}
```

### 3.2 矩阵装配

**文件**：`crates/transient/src/assembly.rs`

```rust
/// 装配三个系统矩阵
pub fn assemble_matrices(
    mesh: &RemMesh,
    domain_map: &DomainMap,
    basis: &NedelecBasis,
) -> (CsrMatrix, CsrMatrix, CsrMatrix) {
    // 返回 (K, M_eps, C_sigma)
    // 遍历体单元 → 计算局部矩阵 → 累加到 triplet → 转 CSR
}
```

局部刚度矩阵（逐元素）：

```rust
// 对四面体 e 的 6 条棱 i, j:
K_local[i][j] = |J_e|/6 * mu_inv_e * dot(curl_phi_i, curl_phi_j)
```

局部质量矩阵：

```rust
M_local[i][j] = |J_e| * eps_e * ∫ phi_i · phi_j dΩ  // 数值积分
```

### 3.3 时间推进器

**文件**：`crates/transient/src/time_stepper.rs`

```rust
pub struct NewmarkStepper {
    pub beta: f64,   // 0.25
    pub gamma: f64,  // 0.5
    pub dt: f64,
    pub n_steps: usize,
    pub a_eff: CsrMatrix,  // 预组装的有效系统矩阵
}

impl NewmarkStepper {
    pub fn step(&self, state: &mut TimeDomainState, f_n: &[f64], f_n1: &[f64]) {
        // 计算 F_eff，调用 PCG 求解 A_eff * E_{n+1} = F_eff
    }
}

pub struct TimeDomainState {
    pub e: Vec<f64>,    // E 场（边自由度）
    pub e_dot: Vec<f64>,
    pub e_ddot: Vec<f64>,
}
```

### 3.4 激励源

**文件**：`crates/transient/src/excitation.rs`

支持三种激励类型：

1. **高斯脉冲**（宽带激励，适合 S 参数提取）
   ```
   J(t) = J₀ · exp(−(t − t₀)² / (2σ²))
   ```

2. **正弦波**（单频稳态响应）
   ```
   J(t) = J₀ · sin(2πf₀t)
   ```

3. **平面波入射**（散射问题）
   ```
   E_inc(r, t) = E₀ · f(t − k̂·r/c)
   ```

### 3.5 吸收边界条件

**文件**：`crates/transient/src/absorbing_bc.rs`

**v1.0：Mur 一阶 ABC**

在外边界 ∂Ω 施加：

```
∂E/∂t + c ∂E/∂n = 0
```

弱形式贡献到边界积分：

```
∫_{∂Ω} (1/c) ∂E/∂t · v dS
```

这为质量矩阵增加一个边界贡献项：

```
M_mur[i][j] = (1/c) ∫_{∂Ω} phi_i · phi_j dS
```

**v1.1（后续）：PML（完美匹配层）**

通过坐标拉伸实现，无反射吸收。

---

## 4. 新建 Crate 结构

```
crates/transient/
├── Cargo.toml
└── src/
    ├── lib.rs              # 入口，pub fn run(config, comm) -> RemResult<()>
    ├── nedelec.rs          # Nedelec P1 基函数、DOF 编号
    ├── assembly.rs         # 矩阵装配 (K, M_ε, C_σ, M_mur)
    ├── time_stepper.rs     # Newmark-β 推进器
    ├── excitation.rs       # 激励源（高斯脉冲、正弦、平面波）
    ├── absorbing_bc.rs     # Mur ABC 边界条件
    └── output.rs           # 逐时间步 VTK 输出、能量 CSV
```

**Cargo.toml 依赖**：

```toml
[dependencies]
rem-core    = { path = "../core" }
rem-config  = { path = "../config" }
rem-mesh    = { path = "../mesh" }
rem-materials = { path = "../materials" }
rem-result  = { path = "../result" }
rem-parallel = { path = "../parallel" }
```

---

## 5. 配置扩展

扩展 [config/src/schema.rs](crates/config/src/schema.rs) 中的 `TransientConfig`：

```rust
#[derive(Deserialize, Debug, Clone)]
pub struct TransientConfig {
    pub dt: f64,                     // 时间步长 [s]
    pub t_end: f64,                  // 终止时间 [s]
    pub output_interval: usize,      // 每 N 步输出一次 VTK
    pub time_scheme: TimeScheme,     // "newmark" | "leapfrog"
    pub excitation: ExcitationType,  // "gaussian" | "sinusoidal" | "planewave"
    pub absorbing_bc: AbcType,       // "mur1" | "pml"
    pub newmark_beta: Option<f64>,   // 默认 0.25
    pub newmark_gamma: Option<f64>,  // 默认 0.5
}
```

Palace JSON 配置示例：

```json
{
  "Problem": { "Type": "Transient" },
  "Model": { "Mesh": "cavity.msh" },
  "Transient": {
    "dt": 1e-11,
    "t_end": 5e-9,
    "output_interval": 10,
    "time_scheme": "newmark",
    "excitation": "gaussian",
    "absorbing_bc": "mur1"
  }
}
```

---

## 6. CLI 接入

修改 [cli/src/main.rs:62](crates/cli/src/main.rs#L62)：

```rust
ProblemType::Transient => {
    rem_transient::run(&config, &comm)
}
```

---

## 7. WASM 兼容性

本方案与现有 WASM 架构天然兼容：

- 纯 Rust 实现，无 C FFI 依赖
- 复用 `parallel::NoComm`（单线程 WASM 环境）
- PCG 求解器已有 WASM 版本
- 输出通过 Blob URL 而非文件系统

WASM 绑定扩展（`crates/wasm/src/lib.rs`）：

```rust
#[wasm_bindgen]
pub fn run_transient(config_json: &str, mesh_bytes: &[u8]) -> JsValue {
    // 调用 rem_transient::run(...)
}
```

---

## 8. 验证方案

| 测试用例 | 分析解/参考解 | 验证量 |
|---------|-------------|--------|
| 平面波在均匀介质中传播 | E(x,t) = E₀sin(kx − ωt) | 波形、相速度误差 < 1% |
| 矩形腔谐振频率 | f_mnp = c/(2)√(m/a)²+(n/b)²+(p/d)² | 频率误差 < 0.5% |
| Mur ABC 反射系数 | R < −40 dB | 边界吸收效果 |
| PEC 球散射（与 MoM 对比） | Mie 级数解析解 | RCS 误差 < 2 dB |

---

## 9. 实施计划

| 阶段 | 内容 | 工作量 | 依赖 |
|------|------|--------|------|
| P1 | 配置扩展 + Crate 骨架 | 1 天 | 无 |
| P2 | **Nedelec P1 基函数**（关键路径）| 5–7 天 | 无 |
| P3 | 矩阵装配（K, M, C） | 3–4 天 | P2 完成 |
| P4 | Newmark-β 时间推进 | 2–3 天 | P3 完成 |
| P5 | 激励源 | 2 天 | P4 完成 |
| P6 | Mur ABC | 2 天 | P4 完成 |
| P7 | 输出 + WASM 绑定 | 2 天 | P4 完成 |
| P8 | 验证与调试 | 3–5 天 | P5-P7 完成 |
| **合计** | | **20–26 天** | |

> Nedelec 基函数（P2）是唯一不可绕过的技术壁垒，建议优先突破。

---

## 10. 与现有求解器的对比

| 维度 | 频域求解器 (Driven) | 时域求解器 (FDTD/TD-FEM) |
|------|-------------------|------------------------|
| 适用场景 | 窄带、稳态 | 宽带、瞬态 |
| 计算量 | 每频点一次求解 | N_steps 次求解 |
| 激励形式 | 单频端口激励 | 脉冲/时变激励 |
| 色散建模 | 频域直接处理 | 需辅助差分方程 (ADE) |
| 非线性支持 | 困难 | 天然支持（逐步更新）|
| S 参数提取 | 直接 | 需傅里叶变换后处理 |

---

## 参考文献

1. Jin, J.-M. *The Finite Element Method in Electromagnetics*, 3rd ed. Wiley-IEEE Press, 2014.
2. Nédélec, J.-C. "A new family of mixed finite elements in ℝ³." *Numerische Mathematik* 50, 57–81 (1986).
3. Newmark, N.M. "A method of computation for structural dynamics." *J. Eng. Mech. Div.* 85, 67–94 (1959).
4. Mur, G. "Absorbing boundary conditions for the finite-difference approximation of the time-domain electromagnetic-field equations." *IEEE Trans. EMC* 23(4), 377–382 (1981).
