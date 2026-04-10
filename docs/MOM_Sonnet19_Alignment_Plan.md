# REM MoM → Sonnet 19 功能对齐开发计划

> 版本：2026-04-10  
> 基准：REM v0.16.0（阶段15完成）  
> 目标：使 REM MoM 具备 Sonnet Suite 19 的平面电路分析核心能力  
> 来源分析：[docs/MOM_vs_Sonnet19.md](MOM_vs_Sonnet19.md)

---

## 执行摘要

Sonnet 19 的核心竞争力集中在**四个维度**，REM 当前均缺失：

| 差距维度 | Sonnet 核心能力 | REM 现状 | 优先级 |
|---------|---------------|---------|--------|
| **Green 函数** | Sommerfeld 分层介质积分 | 仅自由空间 | P0（物理根基） |
| **端口 + S 参数** | 集总/波导端口，Touchstone 输出 | MoM 无端口 | P0（用户接口） |
| **快速算法** | FFT O(N log N) 加速 | ACA O(N·r)（三维散射） | P1（规模扩展） |
| **导体损耗** | 有限电导率 SIBC | 仅 PEC | P1（材料完整性） |

**策略**：按依赖顺序分 4 个版本迭代，每个版本产出可独立验证的里程碑。

---

## 版本路线图

```
v0.16.0（当前）
    │
    ▼
v0.17.0 ── MoM 端口激励 + S 参数（集总端口）
    │
    ▼
v0.18.0 ── 分层介质 Green 函数（Sommerfeld / 离散复像法）
    │
    ▼
v0.19.0 ── 有损导体 SIBC + FFT 加速 MoM（平面结构）
    │
    ▼
v0.20.0 ── MoM AMR + 参数化扫描 + Touchstone 完整兼容
```

---

## 阶段 16：MoM 端口激励 + S 参数提取（v0.17.0）

> **目标**：在 MoM 框架内引入集总端口，输出 Touchstone `.s2p`，
> 使 REM 能仿真微带/CPW 无源器件的 S 参数。

### 16.1 动机与技术路径

Sonnet 通过端口电流激励替代平面波。端口 p 在 MoM 系统中表现为：
- 在端口边上施加已知切向电场（Dirichlet 激励）
- 后处理阶段从端口处的电流/电压提取 S 参数

对于三维 MoM（目前 REM 的形式），集总端口等效为：
1. 在端口面上增加激励 RHS 项（端口电压 V₀ → ∂φ/∂n = j·k·V₀）
2. 端口电流：I_p = ∮_{Γ_p} J·dl（沿端口边积分）
3. S 参数：S_{pq} = (V_p - Z₀ I_p)/(V_q)（扫描各激励端口）

### 16.2 配置扩展

扩展 `crates/config/src/schema.rs`，在 `MomSolverConfig` 中新增端口支持：

```rust
#[derive(Debug, Clone, Deserialize)]
pub struct MomSolverConfig {
    // 原有字段不变 ...

    /// 端口列表（可选；有端口时输出 S 参数，无端口时仍输出 RCS）
    #[serde(rename = "Ports", default)]
    pub ports: Vec<MomPort>,

    /// S 参数参考阻抗 [Ω]，默认 50.0
    #[serde(rename = "RefImpedance", default = "default_ref_impedance")]
    pub ref_impedance: f64,
}

#[derive(Debug, Clone, Deserialize)]
pub struct MomPort {
    /// 端口编号（1-based，与 Touchstone 端口顺序对应）
    #[serde(rename = "Index")]
    pub index: u32,

    /// 端口所在边界的物理组属性 ID
    #[serde(rename = "Attributes")]
    pub attributes: Vec<u32>,

    /// 端口朝向（激励场方向），默认 "x"
    #[serde(rename = "Direction", default = "default_port_direction")]
    pub direction: String,

    /// 端口参考阻抗 [Ω]（覆盖全局 RefImpedance）
    #[serde(rename = "Impedance", default)]
    pub impedance: Option<f64>,
}

fn default_ref_impedance() -> f64 { 50.0 }
fn default_port_direction() -> String { "x".to_string() }
```

### 16.3 核心文件新增

**`crates/mom/src/port.rs`**：

```rust
/// MoM 集总端口模型
pub struct MomLumpedPort {
    pub index: u32,
    /// 端口边上的 RWG 基函数索引列表
    pub rwg_indices: Vec<usize>,
    /// 端口激励向量（仅在该端口激励时非零）
    pub excitation: Vec<Complex64>,
    /// 参考阻抗
    pub z0: f64,
}

impl MomLumpedPort {
    /// 从表面网格和端口属性列表中识别端口 RWG 基函数
    pub fn from_surface(
        surf: &SurfaceMesh,
        port_cfg: &MomPort,
        ref_z0: f64,
    ) -> RemResult<Self>;

    /// 端口激励 RHS 向量（入射场设为端口电压 V₀=1V）
    pub fn excitation_rhs(&self, k: f64, v0: Complex64) -> Vec<Complex64>;

    /// 从电流系数提取端口电压、电流
    /// V = -∫_port E_tan · dl，I = ∮_port J · dl
    pub fn extract_vi(&self, currents: &[Complex64], surf: &SurfaceMesh) -> (Complex64, Complex64);
}
```

**`crates/mom/src/sparams.rs`**：

```rust
/// 从 N 端口激励结果中计算 S 矩阵
///
/// 算法：
///   对每个激励端口 p 运行一次完整 MoM 求解 → 得电流 I_p
///   计算所有端口 q 的 V_q, I_q
///   S_{qp} = (V_q - Z₀_q I_q) / (V_p_incident)
pub fn compute_s_matrix(
    surf: &SurfaceMesh,
    ports: &[MomLumpedPort],
    z_mat: &DMatrix<Complex64>,
    freq: f64,
) -> RemResult<SMatrix>;

/// S 参数矩阵（N×N 复数）
pub struct SMatrix {
    pub n_ports: usize,
    pub freq_hz: f64,
    pub data: Vec<Vec<Complex64>>,  // data[i][j] = S_{i+1, j+1}
}

impl SMatrix {
    /// 写入 Touchstone 格式 .s{N}p 文件
    pub fn write_touchstone(&self, path: &Path) -> RemResult<()>;
    /// 追加写入 Palace 兼容 port-S.csv
    pub fn append_palace_csv(&self, path: &Path) -> RemResult<()>;
}
```

### 16.4 主流程修改（`crates/mom/src/lib.rs`）

```
当 mom_cfg.ports 非空时：
  1. 构建所有端口：Vec<MomLumpedPort>
  2. 对每个激励端口 p:
       rhs = port_p.excitation_rhs(k, V0=1.0)
       currents_p = solve(z_mat, rhs)
  3. compute_s_matrix(surf, ports, currents_p_list) → SMatrix
  4. SMatrix::write_touchstone(output/s_params.s{N}p)
  5. SMatrix::append_palace_csv(output/postpro/port-S.csv)

当 ports 为空时（原有路径）：
  平面波激励 → RCS 输出（原路径不变）
```

### 16.5 配置示例

```json
{
  "Problem": { "Type": "MoM", "Output": "output/microstrip_filter" },
  "Model": { "Mesh": "filter.msh" },
  "Boundaries": { "PEC": { "Attributes": [1, 2, 3] } },
  "Solver": {
    "MoM": {
      "Equation": "CFIE",
      "FreqMin": 1.0e9, "FreqMax": 10.0e9, "FreqStep": 0.5e9,
      "FastSolver": "ACA",
      "RefImpedance": 50.0,
      "Ports": [
        { "Index": 1, "Attributes": [10], "Direction": "x" },
        { "Index": 2, "Attributes": [11], "Direction": "x" }
      ]
    }
  }
}
```

### 16.6 验证基准

| 测试 | 目标 | 方法 |
|------|------|------|
| 单端口 PEC 偶极子 | S11（输入阻抗）vs 解析解 | 半波偶极子，Z_in ≈ 73+j42 Ω @ f_res |
| 双端口传输线段 | S21 < -0.1 dB（无损 PEC） | λ/2 微带线，前向传输接近 0 dB |
| Touchstone 格式 | `.s2p` 文件可被 ADS/Qucs 读取 | 标准格式验证 |

### 16.7 新增/修改文件

```
crates/mom/src/port.rs        ← 新增（端口模型）
crates/mom/src/sparams.rs     ← 新增（S 参数提取）
crates/mom/src/lib.rs         ← 修改（主流程分支）
crates/mom/src/excitation.rs  ← 修改（端口激励 RHS）
crates/config/src/schema.rs   ← 修改（MomPort, RefImpedance）
```

### 16.8 验收标准

- [ ] `Problem.Type = "MoM"` + `Ports` 配置正确解析，S 参数路由激活
- [ ] 单端口激励：一次 MoM 求解 + S11 提取，输出 `port-S.csv`
- [ ] 双端口：两次激励 → 完整 2×2 S 矩阵，写出 `s_params.s2p`
- [ ] 半波偶极子 S11 vs 解析阻抗误差 < 5%
- [ ] 无端口时原有 RCS 路径不受影响（零回归）
- [ ] `cargo test -p rem-mom` 全部通过

---

## 阶段 17：分层介质 Green 函数（v0.18.0）

> **目标**：引入多层介质 Sommerfeld 积分 Green 函数，使 REM MoM
> 能仿真嵌入 PCB/MMIC 基板中的导体结构，这是 Sonnet 最核心的物理能力。

### 17.1 物理模型

Sommerfeld 积分 Green 函数适用于平面分层介质（基板叠层）：

```
z
│  Air (ε₀, μ₀)       z > z_{N}
│─────────────────────
│  Layer N  (ε_N, μ_N) z_{N-1} < z < z_N
│  ...
│  Layer 1  (ε₁, μ₁)  0 < z < z_1
│─────────────────────
│  Ground plane (PEC)  z = 0  （可选）
```

分层介质 Green 函数：

```
G_A(r, r') = (1/4π) ∫₀^∞ g_A(k_ρ, z, z') J₀(k_ρ ρ) k_ρ dk_ρ

其中 ρ = √((x-x')²+(y-y')²)，g_A 通过传输矩阵法计算
```

### 17.2 新 Crate 规划：`crates/layered_green`

```
crates/layered_green/
├── Cargo.toml
└── src/
    ├── lib.rs               ← 公共接口
    ├── layer.rs             ← 层参数（ε, μ, 厚度）
    ├── transfer_matrix.rs   ← 传输矩阵递推（TMM）
    ├── sommerfeld.rs        ← Sommerfeld 积分（数值 + DCIM）
    └── discrete_image.rs    ← 离散复像法（DCIM）近似
```

**选用离散复像法（DCIM）而非纯数值积分的理由**：

| 方法 | 复杂度 | 精度 | 适用范围 |
|------|--------|------|---------|
| 直接数值 Sommerfeld 积分 | 每次 O(N_q)（N_q~1000 积分点） | 高 | 通用，但极慢 |
| DCIM（Chow 1994，Golub-Welch GPOF） | 一次拟合后 O(1) per pair | 中高 | 均匀介质层，快速 |
| GPOF + DCIM | 一次拟合后 O(M_pole) per pair | 高 | 多极点，宽带 |

**推荐**：实现 DCIM（多级离散复像），准确度满足工程需求且速度提升 100-1000×。

### 17.3 核心数据结构

```rust
/// 单一介质层参数
#[derive(Debug, Clone)]
pub struct DielectricLayer {
    /// 相对介电常数（复数，含损耗：εᵣ(1-j·tanδ)）
    pub eps_r: Complex64,
    /// 相对磁导率
    pub mu_r: Complex64,
    /// 层厚度 [m]（顶层可以为无穷大）
    pub thickness: f64,
}

/// 分层介质栈（从底层到顶层，底层 z=0）
pub struct LayeredMedium {
    pub layers: Vec<DielectricLayer>,
    /// 底部是否为 PEC 地板
    pub bottom_pec: bool,
}

/// 预计算的 DCIM Green 函数（在给定频率和介质栈下有效）
pub struct DcimGreen {
    /// GPOF 拟合极点 s_i
    poles: Vec<Complex64>,
    /// GPOF 拟合留数 a_i
    residues: Vec<Complex64>,
    /// 有效频率范围 [Hz]
    freq_hz: f64,
    /// 有效 ρ 范围 [m]（小 ρ 需奇异提取）
    rho_min: f64,
    rho_max: f64,
}

impl DcimGreen {
    /// 从层叠栈和频率预计算 DCIM 系数（一次性开销）
    pub fn compute(medium: &LayeredMedium, freq: f64) -> RemResult<Self>;

    /// 磁矢量位 Green 函数 G_A(ρ, z, z')
    pub fn g_a(&self, rho: f64, z: f64, z_prime: f64) -> Complex64;

    /// 标量位 Green 函数 G_φ(ρ, z, z')
    pub fn g_phi(&self, rho: f64, z: f64, z_prime: f64) -> Complex64;

    /// 法向分量 ∂G/∂z
    pub fn dg_a_dz(&self, rho: f64, z: f64, z_prime: f64) -> Complex64;
}
```

### 17.4 与 MoM 的集成接口

修改 `crates/mom/src/assemble.rs`，抽象 Green 函数依赖：

```rust
/// Green 函数 trait：自由空间和分层介质均实现此接口
pub trait GreenFunction: Send + Sync {
    /// 标量 Green 函数 G(r, r')
    fn g(&self, r: &[f64; 3], rp: &[f64; 3]) -> Complex64;
    /// ∇G（梯度，用于 MFIE）
    fn grad_g(&self, r: &[f64; 3], rp: &[f64; 3]) -> [Complex64; 3];
}

/// 自由空间（现有）
pub struct FreeSpaceGreen { pub k: f64 }

/// 分层介质（新增）
pub struct LayeredGreen {
    pub dcim: DcimGreen,
    pub k: f64,
}

impl GreenFunction for FreeSpaceGreen { ... }
impl GreenFunction for LayeredGreen { ... }
```

修改 `assemble_efie_rwg`、`assemble_mfie_rwg` 签名：

```rust
pub fn assemble_efie_rwg(
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    green: &dyn GreenFunction,   // ← 替换原有 k: f64 参数
    omega: f64,
    quad: &TriQuad,
) -> RemResult<DMatrix<Complex64>>;
```

### 17.5 配置扩展

在 `MomSolverConfig` 中新增介质层叠配置：

```rust
/// 分层介质基板描述（可选；为空时使用自由空间 Green 函数）
#[serde(rename = "Substrate", default)]
pub substrate: Option<SubstrateConfig>,
```

```rust
#[derive(Debug, Clone, Deserialize)]
pub struct SubstrateConfig {
    /// 从底部（地板）到顶部（空气）依次列出
    pub layers: Vec<LayerConfig>,
    /// 底部是否有 PEC 地板，默认 true
    #[serde(default = "default_true")]
    pub bottom_pec: bool,
}

#[derive(Debug, Clone, Deserialize)]
pub struct LayerConfig {
    pub permittivity: f64,
    pub loss_tangent: f64,     // 默认 0.0
    pub permeability: f64,     // 默认 1.0
    pub thickness: f64,        // [m]，最顶层设为极大值
}
```

**Palace JSON 配置示例**（REM 专有扩展，Palace 忽略 `Substrate` 字段）：

```json
{
  "Problem": { "Type": "MoM", "Output": "output/patch_antenna" },
  "Model": { "Mesh": "patch.msh" },
  "Boundaries": { "PEC": { "Attributes": [1] } },
  "Solver": {
    "MoM": {
      "Equation": "EFIE",
      "FreqMin": 2.4e9, "FreqMax": 2.4e9,
      "Substrate": {
        "BottomPec": true,
        "Layers": [
          { "Permittivity": 4.4, "LossTan": 0.02, "Thickness": 1.6e-3 },
          { "Permittivity": 1.0, "LossTan": 0.0,  "Thickness": 1.0    }
        ]
      },
      "Ports": [
        { "Index": 1, "Attributes": [2], "Direction": "y" }
      ]
    }
  }
}
```

### 17.6 验证基准

| 测试 | 参考解 | 目标精度 |
|------|--------|---------|
| FR4 基板（εᵣ=4.4）上的贴片天线谐振频率 | 解析估算 ± 公开测量数据 | 谐振频率误差 < 1% |
| 半波长微带传输线 S21（PEC 地板） | HFSS/CST 参考 | ΔS21 < 0.2 dB |
| DCIM 残差 vs 数值 Sommerfeld 直接积分 | 数值 Sommerfeld（精度参考） | 相对误差 < 1e-3 |

### 17.7 新增/修改文件

```
crates/layered_green/          ← 新 crate（Green 函数库）
  src/lib.rs
  src/layer.rs
  src/transfer_matrix.rs
  src/sommerfeld.rs
  src/discrete_image.rs
crates/mom/src/green_trait.rs  ← 新增（GreenFunction trait）
crates/mom/src/assemble.rs     ← 修改（使用 trait 替换裸 k 参数）
crates/config/src/schema.rs    ← 修改（SubstrateConfig, LayerConfig）
Cargo.toml（workspace）        ← 添加 crates/layered_green
```

### 17.8 验收标准

- [ ] `crates/layered_green` 编译并通过单元测试（DCIM vs 数值 Sommerfeld < 1e-3）
- [ ] 单层 FR4 基板 + PEC 地板：`DcimGreen::compute` 在 1-10 GHz 无 panic
- [ ] MoM 装配切换：`Substrate` 为空时用 `FreeSpaceGreen`，有 `Substrate` 时用 `LayeredGreen`
- [ ] 贴片天线谐振频率 vs 解析估算误差 < 2%
- [ ] 原有自由空间 PEC 球 RCS 测试不受影响（零回归）
- [ ] `cargo test --workspace` 全部通过

---

## 阶段 18：有损导体 SIBC + FFT 加速（v0.19.0）

> **目标**：引入表面阻抗边界条件（SIBC）和 FFT 加速矩阵-向量积，
> 使 REM MoM 能处理实际铜导体损耗并在大规模平面电路上保持效率。

### 18.1 有损导体（SIBC）

**物理模型**：

对有限电导率 σ 的导体，表面阻抗边界条件（Leontovich SIBC）为：

```
E_tan = Z_s · (n̂ × H_tan)

Z_s = (1+j) / (σ · δ_s)，δ_s = √(2/(ωμσ)) 为趋肤深度
```

**SIBC 修正 CFIE 阻抗矩阵**：

在原有 PEC CFIE 基础上，对角块添加阻抗修正：

```
Z_SIBC = Z_CFIE_PEC + Z_s · M_identity

M_identity[m,n] = ∫_Tm f_m · f_n dS（质量矩阵型积分）
```

### 18.2 新文件

**`crates/mom/src/sibc.rs`**：

```rust
/// 表面阻抗（复数，单位：Ω/sq）
#[derive(Debug, Clone, Copy)]
pub struct SurfaceImpedance {
    pub z_s: Complex64,
}

impl SurfaceImpedance {
    /// 从电导率和频率计算 Z_s
    pub fn from_conductivity(sigma: f64, freq: f64) -> Self {
        let omega = 2.0 * PI * freq;
        let delta_s = (2.0 / (omega * MU0 * sigma)).sqrt();
        let z_s = Complex64::new(1.0, 1.0) / (sigma * delta_s);
        Self { z_s }
    }

    /// 铜导体近似：sigma = 5.8e7 S/m
    pub fn copper(freq: f64) -> Self {
        Self::from_conductivity(5.8e7, freq)
    }
}

/// 将 SIBC 修正叠加到现有阻抗矩阵
pub fn apply_sibc(
    z_mat: &mut DMatrix<Complex64>,
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    z_s: SurfaceImpedance,
    quad: &TriQuad,
);
```

### 18.3 配置扩展

在 `Boundaries.Conductivity` 中扩展 MoM 路径：

```rust
#[derive(Debug, Clone, Deserialize, Default)]
pub struct ConductivityBc {
    pub attributes: Vec<u32>,
    /// 电导率 [S/m]（用于 FEM Q 因子 + MoM SIBC）
    pub sigma: f64,
    /// 导体厚度 [m]（用于厚金属建模，可选）
    #[serde(default)]
    pub thickness: Option<f64>,
}
```

### 18.4 FFT 加速 MoM（平面结构）

Sonnet 的 FFT 加速利用平面结构的**移位不变性**：当所有导体在同一平面 z=0 时，Green 函数仅依赖 (x-x', y-y')，构成卷积形式，可用 FFT 加速。

**适用条件**：
- 所有 RWG 基函数的源/观测点均在同一水平层（z≈常数）
- 网格为矩形或三角形均可，但矩形网格 FFT 更高效

**算法**：

```
Z[m,n] = Z(x_m - x_n, y_m - y_n)  ← 仅平面结构成立

矩阵-向量积：Z · I = IFFT(FFT(Z_kernel) × FFT(I))

复杂度：O(N log N)（vs 稠密 O(N²)）
内存：O(N)（仅存 FFT kernel，无需存完整 Z）
```

**实现规划**：

```
crates/mom/src/fft_accel.rs  ← 新增

pub struct FftMomSolver {
    /// 预计算的 FFT 化 Green 函数核
    z_kernel_fft: ndarray::Array2<Complex64>,
    /// 平面网格参数（dx, dy, nx, ny）
    grid: PlanarGrid,
}

impl FftMomSolver {
    /// 检测当前表面网格是否满足平面 FFT 条件
    pub fn is_applicable(surf: &SurfaceMesh, tol: f64) -> bool;

    /// 矩阵-向量积（O(N log N)）
    pub fn matvec(&self, x: &[Complex64]) -> Vec<Complex64>;
}
```

**激活条件**：`FastSolver: "FFT"` + 表面网格满足平面性检测。

### 18.5 验证基准

| 测试 | 参考解 | 目标精度 |
|------|--------|---------|
| 铜微带传输线 S21（有损） | Sonnet 19 参考值 | ΔS21 < 0.1 dB @ 10 GHz |
| 贴片天线 Q 因子（铜 σ=5.8e7） | Sonnet 19 参考 | ΔQ_c/Q_c < 5% |
| FFT 加速 vs 直接装配（N=2000 平面网格） | 直接装配（参考精度） | 结果误差 < 0.1 dB，速度提升 > 5× |

### 18.6 新增/修改文件

```
crates/mom/src/sibc.rs        ← 新增（SIBC 修正）
crates/mom/src/fft_accel.rs   ← 新增（FFT 加速求解器）
crates/mom/src/lib.rs         ← 修改（激活 SIBC/FFT 路径）
crates/config/src/schema.rs   ← 修改（ConductivityBc 扩展）
```

### 18.7 验收标准

- [ ] 铜导体 SIBC：趋肤深度在 1-30 GHz 范围计算正确（vs 解析式）
- [ ] SIBC 修正后阻抗矩阵对角块变化满足 Z_s 量级预期
- [ ] FFT 加速：N=2000 平面网格，速度比稠密至少 5×
- [ ] FFT 加速：与稠密装配的 S 参数误差 < 0.1 dB
- [ ] `cargo test -p rem-mom` 全部通过

---

## 阶段 19：MoM AMR + 参数化扫描 + 完整 Touchstone（v0.20.0）

> **目标**：补齐 Sonnet 在工程可用性上的剩余差距：
> 自适应网格、参数化频率扫描、完整 Touchstone 格式兼容。

### 19.1 MoM 自适应网格细化（AMR）

**误差指示器**：对每个三角面元，计算表面电流密度梯度作为误差指示：

```
η_m = ||∇ J_s||_{T_m} × h_m

其中 h_m = √(A_m) 为面片等效尺寸
```

**细化策略**：
- Dörfler 标记（标记误差贡献排名前 30% 的面元）
- 面元中线分割（Tri3 → 4×Tri3）
- RWG 基函数重映射（细化后重新生成 edges 拓扑）

```
crates/mom/src/amr.rs  ← 新增

pub fn mom_amr_loop(
    config: &PalaceConfig,
    mom_cfg: &MomSolverConfig,
    mesh: &RemMesh,
    max_iter: usize,
    tol: f64,           // RCS/S 参数收敛容限
) -> RemResult<MomResult>;
```

### 19.2 参数化频率扫描优化

引入 MoM 版本的快照 ROM 加速：

```rust
/// MoM 频率扫描 ROM（原理与 Driven FEM ROM 相同）
/// 在若干"锚点"频率求解完整系统，在其余频率用低维近似
pub struct MomRom {
    /// 锚点频率的 LU 分解（O(1) 复用）
    snapshots: Vec<(f64, Vec<Complex64>)>,  // (freq, currents)
    /// 压缩基矩阵 V（N × r）
    basis: DMatrix<Complex64>,
}
```

配置扩展：

```json
"Solver": {
  "MoM": {
    "RomOrder": 8,         // 0 = 禁用；ROM 锚点数
    "FreqMin": 1e9, "FreqMax": 10e9, "FreqStep": 0.1e9
  }
}
```

### 19.3 完整 Touchstone 格式

完善 `crates/mom/src/sparams.rs` 的 Touchstone 输出：

| 功能 | 说明 |
|------|------|
| `.s1p` / `.s2p` / `.sNp` | 自动选择端口数对应后缀 |
| 频率单位 | `GHz` / `MHz` / `Hz` 可配 |
| 数据格式 | `MA`（模+角度）/ `RI`（实+虚）/ `DB`（dB+角度）|
| 选项行 | 标准 `# GHz S MA R 50` 格式 |
| 多端口排列 | N×N 矩阵按 Touchstone 2.0 规范 |
| 注释行 | `!` 前缀，包含 REM 版本和仿真参数 |

### 19.4 验收标准

- [ ] AMR：PEC 球 RCS 在 3 次迭代内收敛（RCS 变化 < 0.1 dB）
- [ ] ROM：10 倍频率扫描点（100 频率点 vs 10 锚点），S 参数误差 < 0.05 dB
- [ ] Touchstone：生成的 `.s2p` 可被 ADS、Qucs、Python scikit-rf 正确读取
- [ ] `cargo test --workspace` 全部通过
- [ ] FEATURE_COMPARISON.md 中 MoM 对应条目全部更新为 ✅

---

## 工作量与时间估计

| 阶段 | 版本 | 核心工作量 | 估计工期 |
|------|------|-----------|---------|
| 阶段 16：端口 + S 参数 | v0.17.0 | 中等（建立在现有 MoM 上，新增 2 个文件） | 2-3 周 |
| 阶段 17：分层介质 Green 函数 | v0.18.0 | 高（新 crate，数学最复杂） | 4-6 周 |
| 阶段 18：SIBC + FFT 加速 | v0.19.0 | 中等（SIBC 简单；FFT 需 ndarray 依赖） | 3-4 周 |
| 阶段 19：AMR + ROM + Touchstone | v0.20.0 | 中等（复用现有 FEM AMR/ROM 模式） | 2-3 周 |
| **合计** | v0.17-0.20 | — | **11-16 周** |

---

## 依赖与风险

### 依赖关系

```
阶段 16（端口/S 参数）
    └─ 依赖：现有 MoM CFIE（阶段 11 ✅）

阶段 17（分层 Green 函数）
    ├─ 依赖：阶段 16（端口验证需要 S 参数输出）
    └─ 外部：ndarray 或 rustfft crate（FFT 库）

阶段 18（SIBC + FFT）
    ├─ 依赖：阶段 17（分层 Green 函数作为物理基础）
    └─ 外部：rustfft（FFT 加速）

阶段 19（AMR + ROM + Touchstone）
    └─ 依赖：阶段 16-18 全部完成
```

### 主要风险

| 风险 | 概率 | 影响 | 缓解措施 |
|------|------|------|---------|
| DCIM 数值稳定性（Sommerfeld 积分 GPOF 拟合失败） | 中 | 高 | 保留直接数值积分作后备；分段频率区间拟合 |
| FFT 加速仅对严格平面网格有效，三维散射体退化 | 低 | 低 | 运行时检测平面性，非平面时自动降级 ACA |
| 端口定义在三维任意网格上不直观 | 中 | 中 | 文档中说明端口仅推荐平面结构使用 |
| WASM 内存不足（分层 Green 函数预计算开销） | 中 | 中 | WASM 下禁用 DCIM，仅支持自由空间（现有路径） |

---

## 与现有文档的关联

本计划直接扩展 `DESIGN_DEV.md`：
- 阶段 16 = DESIGN_DEV.md 阶段 16（接续阶段 15 v0.16.0）
- 阶段 17 = DESIGN_DEV.md 阶段 17
- 阶段 18 = DESIGN_DEV.md 阶段 18
- 阶段 19 = DESIGN_DEV.md 阶段 19

版本里程碑更新：

| 版本 | 内容 | 对应阶段 |
|------|------|---------|
| v0.17.0 | MoM 集总端口 + S 参数 + Touchstone 基础 | 阶段 16 |
| v0.18.0 | 分层介质 Green 函数（DCIM）+ 基板配置 | 阶段 17 |
| v0.19.0 | SIBC 有损导体 + FFT 加速平面 MoM | 阶段 18 |
| v0.20.0 | MoM AMR + 频率扫描 ROM + Touchstone 完整 | 阶段 19 |

---

## 参考资料

- Harrington, *Field Computation by Moment Methods*, IEEE Press, 1993
- Fang, *Analytical and Numerical Methods in Electromagnetic Wave Theory*, 1993（DCIM 基础）
- Chow, Mosig, et al., "Discrete complex image representation," *IEEE TAP*, 1994
- Gustavsen & Semlyen, "Rational approximation of frequency domain responses," *IEEE TPWRD*, 1999
- Touchstone File Format Specification, Version 2.0, IPC-2141, 2009
- [Sonnet Software 白皮书：分层 Green 函数](https://www.sonnetsoftware.com)
- REM `crates/mom/src/` — 现有 MoM 实现基准
- REM `docs/MOM_vs_Sonnet19.md` — 差距分析原文
