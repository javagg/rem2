# MoM/BEM 增量推进计划

> 项目：rem2 电磁仿真工具  
> 文档版本：v1.0，2026-04-05  
> 前置状态：v0.2 已完成（FEM 静电/静磁/特征模/频域驱动，Palace 配置兼容）

---

## 一、目标与范围

在不破坏任何现有 Palace 兼容性的前提下，为 rem2 新增：

1. **矩量法（MoM）**：RWG 基函数 + CFIE，用于 PEC/介质体的全波散射和天线分析
2. **边界元法（BEM）**：Laplace/Helmholtz 积分方程，与 FEM 静电求解器交叉验证

**不在本计划内**：FMM/ACA 快速算法（可作 v1.0 后续项），介质 MoM（PMCHWT，作为 v0.7 选项）。

---

## 二、Palace 兼容性约束（硬性要求）

> 本计划所有阶段必须满足以下约束，任何 PR 合并前验证。

| 约束 | 验证方式 |
|------|---------|
| 现有 Palace 示例配置解析不变 | `cargo test -p rem-config` |
| FEM 静电/静磁/特征模/驱动结果不变 | `tests/integration/` 回归测试 |
| 新增字段全部有默认值（`#[serde(default)]`） | 编译器检查 + 单测 |
| `ProblemType` 已有枚举值处理逻辑不修改 | Code review |
| WASM 编译不引入新的 C FFI 依赖 | `cargo build --target wasm32-unknown-unknown` |

---

## 三、总体里程碑

```
v0.3  [已有] Palace 兼容测试
 │
v0.4  [Phase A] MoM 基础设施（4 周）
 │     ├─ crates/mom 骨架
 │     ├─ 配置扩展（MoM/BEM Problem.Type）
 │     ├─ 表面网格提取
 │     └─ 高斯求积规则
 │
v0.5  [Phase B] EFIE 脉冲基函数（3 周）
 │     ├─ Green 函数
 │     ├─ 密集复数矩阵装配
 │     ├─ Duffy 自积分
 │     └─ 端到端验证：PEC 球体 RCS
 │
v0.6  [Phase C] RWG + CFIE（5 周）
 │     ├─ RWG 基函数
 │     ├─ Sauter-Schwab 奇异积分（4 种情形）
 │     ├─ CFIE 阻抗矩阵
 │     └─ 生产级 RCS 精度（< 5% vs Mie）
 │
v0.7  [Phase D] Laplace BEM（3 周）
       ├─ BEM 积分算子
       ├─ P0/P1 基函数
       └─ 与 FEM 静电交叉验证（< 1% 误差）
```

**总计：约 15 周**（Phase A-D 顺序执行，不可并行，因为每阶段是下一阶段的前置）

---

## 四、Phase A：MoM 基础设施（v0.4.0）

**时长**：4 周  
**负责人**：主开发  
**风险**：低（纯基础设施，无数学难点）

### 任务列表

#### W1：crate 骨架 + 配置扩展

- [ ] **A1** 创建 `crates/mom/Cargo.toml`，加入 workspace
- [ ] **A2** 扩展 `crates/config/src/schema.rs`：
  - 新增 `ProblemType::MoM` 和 `ProblemType::BEM`
  - 新增 `SolverConfig.mom: Option<MomSolverConfig>`
  - 新增 `PalaceConfig.postprocessing: Postprocessing`（含 `RcsConfig`）
  - 所有新字段加 `#[serde(default)]`
- [ ] **A3** 更新 `crates/cli/src/main.rs` 的 `match` 分支，添加 `MoM`/`BEM` 路由（返回 `not yet implemented` 错误）
- [ ] **A4** 配置解析测试：解析 `examples/sphere_mom.json`（新建示例文件）

**A4 示例文件** (`examples/sphere_mom.json`):
```json
{
  "Problem": { "Type": "MoM", "Output": "./output/sphere" },
  "Model": { "Mesh": "examples/meshes/sphere.msh", "L0": 1.0e-3 },
  "Boundaries": { "PEC": { "Attributes": [1] } },
  "Solver": {
    "MoM": {
      "Equation": "CFIE", "Basis": "RWG",
      "FreqMin": 1.0e9, "FreqMax": 1.0e9, "FreqStep": 1.0e9,
      "Alpha": 0.5
    }
  },
  "Postprocessing": {
    "RCS": { "PhiDeg": [0], "ThetaDeg": "0:10:180" }
  }
}
```

**验收**：
- `cargo test -p rem-config` 全部通过（含新增配置测试）
- 现有 Palace 配置回归测试不变

#### W2：表面网格提取

- [ ] **A5** 实现 `crates/mom/src/surface_mesh.rs`：
  - `SurfaceMesh::extract(rem_mesh, pec_attrs)` — 从 `RemMesh.boundary_elements` 过滤 PEC 面
  - `TriFace`：节点索引、质心、单位外法向、面积
  - `build_edge_topology(faces)` — 构建共享边列表（排序节点对 hash，O(E log E)）
  - `SharedEdge`：两端节点、T±面索引、边长

**验收**：
- 单测：从 10×10 球面网格提取，面元数、内部边数（= 顶点数 + 面元数 - 2，欧拉公式）符合预期
- 单测：所有内部边的 T+/T- 法向指向一致（外向）

#### W3-W4：高斯求积规则

- [ ] **A6** 实现 `crates/mom/src/quadrature.rs`：
  - Dunavant 规则阶次 1/3/5/7/9（直接内嵌求积点和权重，不做运行时计算）
  - `TriQuad::new(order)` 工厂函数
  - `integrate_scalar(face, nodes, quad, f)` 辅助函数
  - 重心坐标 → 全局坐标映射

**验收**：
- 单测：对面积为 1 的参考三角形积分常数 f=1，结果 = 1.0，误差 < 1e-14
- 单测：7 点规则对完整 4 次多项式精确积分（误差 < 1e-12）
- 单测：高斯点坐标在三角形内（所有重心坐标 > 0 且和 = 1）

---

## 五、Phase B：EFIE 脉冲基函数（v0.5.0）

**时长**：3 周  
**负责人**：主开发  
**风险**：中（Duffy 变换实现有一定难度，但有大量参考文献）

### 任务列表

#### W5：Green 函数 + 密集矩阵

- [ ] **B1** 实现 `crates/mom/src/green.rs`：
  - `green3d(r, r_prime, k) -> Complex64`
  - `green3d_normal_deriv(r, r_prime, n_prime, k) -> Complex64`
  - 单测：|G| 正比于 1/r，相位线性增长（各自对比 r = 0.1, 1.0, 10.0 λ）
- [ ] **B2** 添加 `faer` 依赖（`crates/mom/Cargo.toml`）
- [ ] **B3** 实现 `crates/mom/src/matrix.rs`：
  - `DenseComplexMat`（薄包装 `faer::Mat<Complex64>`）
  - 从 Z 矩阵求解：`lu_solve(z: &DenseComplexMat, rhs: &[Complex64]) -> Vec<Complex64>`

#### W6：Duffy 自积分

- [ ] **B4** 实现 `crates/mom/src/singular.rs`：
  - `zmn_self_duffy(face, nodes, k, omega, n_quad) -> Complex64`
  - 实现：将 T 分 3 个子三角（以 r₁/r₂/r₃ 为极点），每个做极坐标变换（ρ, θ），标准高斯积分
  - 单测：对比已知解析值（均匀电流板的自阻抗）

- [ ] **B5** 实现 `crates/mom/src/assemble.rs`（脉冲基函数版）：
  - `assemble_efie_pulse(surf, freq, quad, singular_tol) -> faer::Mat<Complex64>`
  - `rayon` 并行外层循环
  - 对角块调用 `zmn_self_duffy`，非对角块调用标准高斯积分

#### W7：端到端验证

- [ ] **B6** 实现 `crates/mom/src/excitation.rs`：
  - `plane_wave_rhs(surf, k, pol, inc_dir) -> Vec<Complex64>`：平面波激励向量
- [ ] **B7** 实现 `crates/mom/src/postprocess.rs`：
  - `compute_rcs_pulse(currents, surf, freq, theta_deg, phi_deg) -> Vec<(f64, f64, f64)>`
  - `write_rcs_csv(output_dir, freq, data) -> RemResult<()>`（Palace 扩展 CSV 格式）
- [ ] **B8** 实现 `crates/mom/src/lib.rs`：
  - `pub fn run(config: &PalaceConfig) -> RemResult<()>` 主入口
- [ ] **B9** 准备 PEC 球面网格（`examples/meshes/sphere_100.msh`，约 100 面元）
- [ ] **B10** 集成测试：`tests/integration/test_mom_pulse.rs`
  - PEC 球体 (r=0.1λ) @ 1 GHz，单站 RCS 与 Mie 误差 < 15%（脉冲基函数的预期精度）

**验收**：
- `rem examples/sphere_mom.json` 运行完成，生成 `rcs.csv`
- 单元测试 `cargo test -p rem-mom` 全部通过
- Palace FEM 回归测试不变

---

## 六、Phase C：RWG 基函数 + CFIE（v0.6.0）

**时长**：5 周  
**负责人**：主开发（建议引入 BEM 专家顾问审查奇异积分）  
**风险**：**高**（Sauter-Schwab 4 种情形是本计划最复杂部分）

### 任务列表

#### W8-W9：RWG 基函数

- [ ] **C1** 实现 `crates/mom/src/basis/rwg.rs`：
  - `RwgBasis` 结构体（`edge_idx`, `plus/minus_face`, `free_node_plus/minus`, `length`）
  - `RwgBasis::eval(r, surf, in_plus) -> [f64; 3]` — 面内矢量值
  - `RwgBasis::divergence(surf, in_plus) -> f64` — 常数 ±lₙ/Aₙ
  - `generate_rwg_bases(surf) -> Vec<RwgBasis>` — 对所有内部共享边生成基函数
  - 单测：T+ 和 T- 上法向分量连续（共享边法向分量相等），不连续分量跳变正确

- [ ] **C2** 实现 `crates/mom/src/assemble.rs` RWG 版（替换脉冲版）：
  - `zmn_rwg_regular(bm, bn, surf, nodes, k, omega, quad) -> Complex64` — 非接触对
  - 奇异情形的分类：计算 T_m 与 T_n 的共享节点数（0/1/2/3）

#### W10-W11：Sauter-Schwab 奇异积分

- [ ] **C3** 实现 `crates/mom/src/singular.rs` 完整版：
  - **情形 0（完全重合，3 次共享节点）**：`zmn_self_rwg_duffy` — 6 个 2D 子块，每块 `n_quad²` 点
  - **情形 1（共享边，2 次共享节点）**：`zmn_edge_sauter_schwab` — Sauter-Schwab Rule 2
  - **情形 2（共享顶点，1 次共享节点）**：`zmn_vertex_sauter_schwab` — Sauter-Schwab Rule 3
  - **情形 3（无接触）**：标准高斯（已有）
  - 单测（每种情形）：对比 BEM++ 或 FEKO 参考值（预先计算存为 fixture）
  - 单测：所有情形连续性（接触程度递减时，积分值平滑过渡）

实现参考：Sauter & Schwab，*Boundary Element Methods*（2011），Chapter 5。

- [ ] **C4** 实现 MFIE 核：
  - `kmn_mfie(bm, bn, surf, k, quad) -> Complex64`（curl-Green 积分）
  - 单测：MFIE 对角块（delta 项）的解析值验证

#### W12：CFIE 装配 + 验证

- [ ] **C5** 实现 `assemble_cfie(surf, freq, bases, alpha, quad) -> faer::Mat<Complex64>`：
  - `Z_CFIE = alpha * Z_EFIE + (1 - alpha) * eta0 * Z_MFIE`
  - 并行装配（`rayon`，按行分块）
- [ ] **C6** 更新 `run()` 使用 RWG + CFIE
- [ ] **C7** 集成测试：`tests/integration/test_mom_cfie.rs`
  - PEC 球体 (r=0.1λ) @ 1 GHz，单站 RCS 与 Mie 解析解误差 < **5%**
  - PEC 球体 (r=0.5λ) @ 1 GHz，避开内谐振（CFIE 优势），与 EFIE 对比稳定性
  - 偶极子天线输入阻抗（与 NEC2 参考值误差 < 3%）

**验收**：
- 主要精度验收：PEC 球体 RCS 误差 < 5%（Mie 解析解）
- `cargo test -p rem-mom` 全部通过
- WASM 编译：`cargo build --target wasm32-unknown-unknown -p rem-wasm` 无错

---

## 七、Phase D：Laplace BEM 静电求解器（v0.7.0）

**时长**：3 周  
**负责人**：主开发  
**风险**：低（比 MoM 简单，Laplace 核无高频奇异性，且有 FEM 参考解可对比）

### 任务列表

#### W13-W14：BEM 算子实现

- [ ] **D1** 创建 `crates/bem/` crate（或在 `crates/mom/` 中新增 `src/laplace.rs`）
- [ ] **D2** 实现 Laplace Green 函数（已有，无 k 参数的 `green3d` 特例）
- [ ] **D3** 实现 P0（常数）基函数版 BEM 算子：
  - `assemble_H(surf, quad)` — 双层势算子 H（法向导数 Green）
  - `assemble_G(surf, quad)` — 单层势算子 G（Green 函数）
  - 对角块奇异处理：P0 自积分有解析公式 `G_self = area/(2π) * (1 - ln(2√(area/π)))`
- [ ] **D4** 施加 Dirichlet 边界条件（PEC/Terminal/Ground）
- [ ] **D5** 求解线性系统（`faer::LU`），输出节点电位和法向通量

#### W15：验证与集成

- [ ] **D6** 集成测试：`tests/integration/test_bem_electrostatic.rs`
  - 球体电容：`C = 4πε₀R`，误差 < 0.5%
  - 平行板电容（宽板近似）：误差 < 1%
  - 与 FEM 静电求解器交叉验证（同一问题，两种方法误差 < 1%）
- [ ] **D7** 输出电容矩阵（与 FEM 静电 `capacitance.csv` 格式相同）
- [ ] **D8** 扩展 Yew demo，支持 `BEM` 问题类型示例（可选，按时间决定）

**验收**：
- 球体电容误差 < 0.5%（Mie 解析解）
- 与 FEM 静电交叉误差 < 1%
- `Problem.Type = "BEM"` 完整路由并正确输出

---

## 八、关键风险与缓解措施

| 风险 | 概率 | 影响 | 缓解措施 |
|------|------|------|---------|
| **Sauter-Schwab 奇异积分数值不稳定** | 高 | 高 | 分情形逐步实现；每情形有参考值对比；可引入 BEM 专家审查 |
| **CFIE 在某些频率发散（内谐振）** | 中 | 中 | 文档明确标注；α=0.5 是标准取值；添加条件数监控日志 |
| **WASM faer LU 性能不足** | 中 | 低 | WASM 限制 N<1000；文档注明；native 模式无限制 |
| **faer API 变化（版本升级）** | 低 | 中 | 锁定版本；维护最小包装层隔离 faer API |
| **fem-rs 网格 API 不兼容 BEM 需求** | 低 | 高 | `SurfaceMesh` 独立于 fem-rs，直接从 `RemMesh` 构建 |
| **Sauter-Schwab 共享边/顶点情形**（最难） | 高 | 高 | Phase C W10-W11 集中处理；先写测试，再写实现（TDD） |

---

## 九、Palace 兼容性守护检查表（每个 PR 必查）

```bash
# 1. 配置兼容测试
cargo test -p rem-config -- --nocapture

# 2. FEM 回归测试（确保 MoM/BEM 代码不影响 FEM 路径）
cargo test -p rem-electrostatic -p rem-magnetostatic -p rem-eigenmode -p rem-driven

# 3. WASM 编译检查
cargo build --target wasm32-unknown-unknown -p rem-wasm --no-default-features --features wasm

# 4. 新增字段默认值检查（Palace 示例不含 MoM 字段）
cargo test -p rem-config -- palace_compat

# 5. clippy
cargo clippy --workspace -- -D warnings
```

---

## 十、输出文件格式（Palace 扩展）

### MoM 新增输出（不覆盖 FEM 输出）

```
{output_dir}/
├── postpro/
│   ├── domain-E.csv        # FEM 已有格式（不变）
│   ├── rcs.csv             # MoM 新增：RCS 方向图
│   └── port-Z.csv          # MoM 新增：端口阻抗（天线输入阻抗）
└── paraview/
    ├── solution.vtu         # FEM 已有（不变）
    └── surface_current.vtu  # MoM 新增：表面电流 VTK
```

**rcs.csv 格式**：
```
Freq (GHz),Theta (deg),Phi (deg),RCS (dBsm)
1.000000e+00,0.0,0.0,-15.30
1.000000e+00,10.0,0.0,-14.85
...
```

**port-Z.csv 格式**（天线输入阻抗）：
```
Freq (GHz),Port,Re(Z) (Ohm),Im(Z) (Ohm),|Z| (Ohm)
1.000000e+00,1,73.10,42.55,84.65
...
```

---

## 十一、参考资料

| 资料 | 用途 |
|------|------|
| Gibson, *The Method of Moments in Electromagnetics* (2021) | RWG + CFIE 主要参考 |
| Rao, Wilton, Glisson, IEEE TAP 1982 | RWG 基函数原始论文 |
| Sauter & Schwab, *Boundary Element Methods* (2011), Ch.5 | Sauter-Schwab 奇异积分规则 |
| Duffy, J. Numer. Anal. 1982 | Duffy 变换原始论文 |
| [BEM++ 源码](https://github.com/bempp/bempp-cl) | 奇异积分 Python 参考实现 |
| [FEKO](https://altair.com/feko/) | MoM 对比验证工具（商业） |
| [NEC2](http://www.nec2.org/) | 天线阻抗对比验证（开源） |
