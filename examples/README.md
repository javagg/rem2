# REM Examples

本目录收录 REM 的所有仿真示例配置，分为两类：

- **Palace 对齐示例**：与 [Palace](https://github.com/awslabs/palace) 上游保持字段对应，用于验证 REM 与 Palace 的行为一致性。
- **REM 独立示例**：REM 自定义，采用程序内生成网格，用于单元验证与 Yew demo。

每个 `.json` 配置文件均为独立示例，可单独运行。

---

## 同步状态说明

| 标记 | 含义 |
|------|------|
| `In Sync` | 配置与网格已对齐上游，记录了 Source Commit |
| `Pending` | 本仓库已有配置，尚未与上游逐字段核对 |
| `REM Only` | 无对应 Palace 上游，REM 独立定义 |
| `Diverged` | 与上游存在可见偏差，已记录原因 |

---

## 示例总览

| Config 文件 | 所属目录 | Problem Type | 网格文件 | Palace 对应 | Sync |
|-------------|----------|--------------|----------|-------------|------|
| `adapter/hybrid.json` | adapter | Eigenmode | `adapter.msh` | palace/examples/adapter | Pending |
| `antenna/antenna_halfwave_dipole.json` | antenna | Driven | `antenna.msh` | palace/examples/antenna | Pending |
| `antenna/antenna_short_dipole.json` | antenna | Driven | `antenna.msh` | palace/examples/antenna | Pending |
| `coaxial/coaxial.json` | coaxial | Electrostatic | `coaxial_ascii.msh` | palace/examples/coaxial | Pending |
| `coaxial/coaxial_matched.json` | coaxial | Transient | `coaxial.msh` | palace/examples/coaxial | Pending |
| `coaxial/coaxial_open.json` | coaxial | Transient | `coaxial.msh` | palace/examples/coaxial | Pending |
| `coaxial/coaxial_short.json` | coaxial | Transient | `coaxial.msh` | palace/examples/coaxial | Pending |
| `cpw/cpw_coax_adaptive.json` | cpw | Driven | `cpw_coax_0.msh` | palace/examples/cpw | Pending |
| `cpw/cpw_coax_uniform.json` | cpw | Driven | `cpw_coax_0.msh` | palace/examples/cpw | Pending |
| `cpw/cpw_lumped_adaptive.json` | cpw | Driven | `cpw_lumped_0.msh` | palace/examples/cpw | Pending |
| `cpw/cpw_lumped_uniform.json` | cpw | Driven | `cpw_lumped_0.msh` | palace/examples/cpw | Pending |
| `cpw/cpw_lumped_eigen.json` | cpw | Eigenmode | `cpw_lumped_0.msh` | palace/examples/cpw | Pending |
| `cpw/cpw_wave_adaptive.json` | cpw | Driven | `cpw_wave_0.msh` | palace/examples/cpw | Pending |
| `cpw/cpw_wave_uniform.json` | cpw | Driven | `cpw_wave_0.msh` | palace/examples/cpw | Pending |
| `cpw/cpw_wave_eigen.json` | cpw | Eigenmode | `cpw_wave_0.msh` | palace/examples/cpw | Pending |
| `cylinder/cavity_pec.json` | cylinder | Eigenmode | `cylinder_hex.msh` | palace/examples/cylinder | Pending |
| `cylinder/cavity_impedance.json` | cylinder | Eigenmode | `cylinder_prism.msh` | palace/examples/cylinder | Pending |
| `cylinder/driven_wave.json` | cylinder | Driven | `cylinder_hex.msh` | palace/examples/cylinder | Pending |
| `cylinder/waveguide.json` | cylinder | Eigenmode | `cylinder_tet.msh` | palace/examples/cylinder | Pending |
| `cylinder/floquet.json` | cylinder | Eigenmode | `cylinder_tet.msh` | palace/examples/cylinder | Pending |
| `transmon/transmon_coarse.json` | transmon | Eigenmode | `transmon.msh2` | palace/examples/transmon | Pending |
| `transmon/transmon_amr.json` | transmon | Eigenmode | `transmon.msh2` | palace/examples/transmon | Pending |
| `spheres/spheres.json` | spheres | Electrostatic | `coaxial_2d.msh` | — | REM Only |
| `rings/rings.json` | rings | Magnetostatic | `slab_2d.msh` | — | REM Only |
| `parallel_plate/parallel_plate.json` | parallel_plate | Electrostatic | `plate_2d.msh` | — | REM Only |
| `cylinder/cylinder.json` | cylinder | Magnetostatic | `cylinder_hex.msh` | — | REM Only |
| `cpw/cpw.json` | cpw | Driven | `cpw_coax.msh` | — | REM Only |
| `sbr_sphere/sbr_sphere.json` | sbr_sphere | SBR | `sphere.msh` | — | REM Only |
| `transmon/transmon.json` | transmon | Eigenmode | `transmon.msh2` | — | REM Only |

---

## Palace 对齐示例

### adapter

**网格**：`adapter/mesh/adapter.msh`
**Palace 参考**：`palace/examples/adapter/`

| Config | Problem | Solver | 关键参数 | 网格 |
|--------|---------|--------|----------|------|
| `hybrid.json` | Eigenmode | Order 2 | N=3, Target=6.6 GHz, WavePort×2, PEC | `adapter.msh` |

---

### antenna

**网格**：`antenna/mesh/antenna.msh`（L0=1.0 m，单位为米）
**Palace 参考**：`palace/examples/antenna/`

| Config | Problem | Solver | 关键参数 | 激励 |
|--------|---------|--------|----------|------|
| `antenna_halfwave_dipole.json` | Driven | Order 2 | f=74.9 MHz, Absorbing[4], PEC[1,2] | LumpedPort[3] +Z, R=50Ω |
| `antenna_short_dipole.json` | Driven | Order 2 | f=74.9 MHz, Absorbing[4] (Order 2) | CurrentDipole，无端口 |

- `antenna_halfwave_dipole`：半波偶极子，Feed 在 gap 处，带 FarField 后处理。
- `antenna_short_dipole`：赫兹短偶极子电流源激励，无物理端口，材质属性含 [5,6,7]。

---

### coaxial

**网格**：`coaxial/mesh/coaxial.msh` / `coaxial_ascii.msh`
**Palace 参考**：`palace/examples/coaxial/`

| Config | Problem | Solver | 材质 ε_r / σ | 端口配置 | 网格 |
|--------|---------|--------|--------------|---------|------|
| `coaxial.json` | Electrostatic | Order 1 | ε=2.1，无损 | Terminal[1], Ground[2] | `coaxial_ascii.msh` |
| `coaxial_matched.json` | Transient | Order 3 | ε=2.08, σ=4.629e-2 | LP[3] 激励 +R, LP[4] 匹配 50Ω | `coaxial.msh` |
| `coaxial_open.json` | Transient | Order 3 | ε=2.08, σ=4.629e-2 | LP[3] 激励 +R, PMC[4] 开路 | `coaxial.msh` |
| `coaxial_short.json` | Transient | Order 3 | ε=2.08, σ=4.629e-2 | LP[3] 激励 +R, PEC[2,4] 短路 | `coaxial.msh` |

Transient 三例均使用 `ModulatedGaussian` 激励，f=10 GHz，MaxTime=1 ns，TimeStep=5 ps，验证匹配/开路/短路反射行为。

> **网格注记**：`coaxial.json` 使用 `coaxial_ascii.msh`（Electrostatic），Transient 三例使用 `coaxial.msh`（含端口面 [3][4]）。Yew demo 默认 `coaxial.msh`，与 Transient 一致；Electrostatic 验证使用 ascii 变体。

---

### cpw（共面波导）

**Palace 参考**：`palace/examples/cpw/`  
三套网格对应三种端口拓扑：

| 网格变体 | 端口类型 | 端口数 | Port Z₀ |
|---------|---------|--------|---------|
| `cpw_coax_0.msh` | LumpedPort (+R/-R 方向) | 4 | 56.02 Ω |
| `cpw_lumped_0.msh` | LumpedPort (+Y 方向) | 4 | 56.02 Ω |
| `cpw_wave_0.msh` | WavePort | 4 | — |

| Config | Problem | 网格 | 频率范围 | 求解策略 |
|--------|---------|------|---------|---------|
| `cpw_coax_adaptive.json` | Driven | `cpw_coax_0.msh` | 2–30 GHz | 自适应，AdaptiveTol=1e-3 |
| `cpw_coax_uniform.json` | Driven | `cpw_coax_0.msh` | 2–30 GHz, Step=2 GHz | 均匀扫频 |
| `cpw_lumped_adaptive.json` | Driven | `cpw_lumped_0.msh` | 2–32 GHz | 自适应，AdaptiveTol=1e-3 |
| `cpw_lumped_uniform.json` | Driven | `cpw_lumped_0.msh` | 2–32 GHz, Step=6 GHz | 均匀扫频 |
| `cpw_lumped_eigen.json` | Eigenmode | `cpw_lumped_0.msh` | Target=16 GHz | N=1 |
| `cpw_wave_adaptive.json` | Driven | `cpw_wave_0.msh` | 2–32 GHz | 自适应，AdaptiveTol=1e-3 |
| `cpw_wave_uniform.json` | Driven | `cpw_wave_0.msh` | 2–32 GHz, Step=6 GHz | Sapphire 各向异性，含 Probe/Energy/Dielectric 后处理 |
| `cpw_wave_eigen.json` | Eigenmode | `cpw_wave_0.msh` | Target=10 GHz | N=1 |

`cpw_wave_uniform` 为最完整配置：材质各向异性（蓝宝石），含 `SurfaceFlux`、`Dielectric`（SA/MS/MA 界面损耗）与 `Probe` 后处理。

---

### cylinder

**Palace 参考**：`palace/examples/cylinder/`  
三套网格（hex/prism/tet）对应不同求解场景：

| Config | Problem | 网格 | L0 | 关键参数 |
|--------|---------|------|----|---------|
| `cavity_pec.json` | Eigenmode | `cylinder_hex.msh` | 1 cm | Order 4, N=15, Target=2 GHz, PEC 腔体 |
| `cavity_impedance.json` | Eigenmode | `cylinder_prism.msh` | 1 cm | Order 4, N=15, Target=2 GHz, Impedance Rs=0.0184 Ω |
| `driven_wave.json` | Driven | `cylinder_hex.msh` | 1 cm | Order 4, WavePort[2], 2.5–5 GHz |
| `waveguide.json` | Eigenmode | `cylinder_tet.msh` | 1 cm | Order 4, N=15, ε=2.08, tan δ=4e-4 |
| `floquet.json` | Eigenmode | `cylinder_tet.msh` | 1 cm | Order 4, N=15, Floquet 周期边界（待验证） |

`cavity_pec` vs `cavity_impedance`：同网格族（hex/prism），验证理想 PEC 腔与有损金属腔的 Q 因子差异。

---

### transmon

**Palace 参考**：`palace/examples/transmon/`  
网格：`transmon/mesh/transmon.msh2`（L0=1 μm）

| Config | Problem | Order | Refinement | 材质 | 端口 |
|--------|---------|-------|------------|------|------|
| `transmon_coarse.json` | Eigenmode | 2 | MaxIter=0（无自适应） | Air[2] + Si[1] ε=9.3 | R[6,7], JJ LC[4] |
| `transmon_amr.json` | Eigenmode | 3 | MaxIter=2（2 轮 AMR） | Air[2] + Si[1] ε=9.3 | R[6,7], JJ LC[4] |

两例均求解 N=2 模，Target=4 GHz，Josephson Junction 等效为 C=5.5 fF / L=14.86 nH 集中参数端口。`transmon_amr` 在 `transmon_coarse` 基础上开启自适应网格加密（AMR），验证收敛性。

---

## REM 独立示例

REM 自定义示例，网格由程序生成或专为 rem 定制，无对应 Palace 上游配置。

### spheres — 同轴截面（Electrostatic）

**Config**：`spheres/spheres.json`  
**网格**：`spheres/mesh/coaxial_2d.msh`（程序生成，`annular_msh(1.0, 4.0, 10, 32, 1, 2, 10)`）  
**类比**：Palace spheres 示例（两球电容矩阵）的 2D 轴对称简化版

几何：内导体 r=1 mm，外导体 r=4 mm，真空填充（ε=1）。  
物理标签：[1]=内导体（Terminal，V=1 V），[2]=外导体（Ground），[10]=介质域。  
解析解：C/L = 2πε₀ / ln(r_o/r_i) ≈ 40.12 pF/m。

---

### rings — 双层平板（Magnetostatic）

**Config**：`rings/rings.json`  
**网格**：`rings/mesh/slab_2d.msh`（程序生成，`rect_bimaterial_msh(1.0, 1.0, 20, 20, 1, 2, 10, 20)`）  
**类比**：Palace rings 示例（两环电感矩阵）的 2D 双材料简化版

几何：1 mm × 1 mm 方形，下半（y∈[0, 0.5mm]）铁 μ_r=1000 [10]，上半空气 μ_r=1 [20]。  
边界：[1]=底边（Ground，A_z=0），[2]=顶边（SurfaceCurrent，A_z=1）。  
解析解：A_z(y=0.5mm) = 1000/1001 ≈ 0.999001，验证界面跳变条件。

---

### parallel_plate — 平行板（Electrostatic）

**Config**：`parallel_plate/parallel_plate.json`  
**网格**：`parallel_plate/mesh/plate_2d.msh`

几何：2D 平行板截面，ε=1（真空）。  
物理标签：[1]=底板（Ground），[2]=顶板（Terminal），[10]=介质域。  
用途：验证平行板电容的解析解 C = ε₀·A/d。

---

### cylinder — 磁化柱体（Magnetostatic）

**Config**：`cylinder/cylinder.json`  
**网格**：`cylinder/mesh/cylinder_hex.msh`（L0=1 mm）  
**注**：同目录的其他 json（cavity_pec 等）为 Palace 对齐示例，本配置为 REM 独立验证用

几何：高导磁率（μ_r=1000）柱体。  
边界：[2]=顶面（SurfaceCurrent），[3]=底面（Ground）。  
用途：验证轴对称磁矢量势求解与磁通连续性。

---

### cpw — REM 基线（Driven）

**Config**：`cpw/cpw.json`  
**网格**：`cpw/mesh/cpw_coax.msh`（L0=1 μm，Silicon ε=11.7）

REM demo 用最小可运行版本，简化自 Palace CPW 示例，单 LumpedPort，频率 4–8 GHz。  
对应 Yew demo 中 "CPW (Driven)" 条目。

---

### sbr_sphere — PEC 球体单站 RCS（SBR+）

**Config**：`sbr_sphere/sbr_sphere.json`  
**网格**：`sbr_sphere/mesh/sphere.msh`（UV 球，半径 0.5 m，24×48 面）

频率：3 GHz，ka ≈ 31.4（光学极限区）。  
激励：平面波，θ_inc=0°，θ 极化。  
后处理：单站 RCS，ThetaDeg 0:5:180。  
期望值：σ_mono ≈ πa² = 0.785 m²（≈ −1.05 dBsm），与 Mie 解析解对比。

---

### transmon — REM 简化版（Eigenmode）

**Config**：`transmon/transmon.json`  
**网格**：`transmon/mesh/transmon.msh2`（L0=1 μm）  
**注**：与 `transmon_coarse.json` / `transmon_amr.json` 不同，此版本为 REM 精简配置，仅含基础材质与 PEC 边界，用于 Yew demo

Order=2，N=5，Target=5 GHz，ε_r=11.4（硅），PEC[2]。

---

## 多网格变体默认选择

| 示例 | 可选网格 | Yew demo 默认 | 选择理由 |
|------|---------|--------------|---------|
| `coaxial` | `coaxial.msh`, `coaxial_ascii.msh` | `coaxial.msh` | Transient 三例共用；Electrostatic 配置单独用 ascii 版 |
| `cpw` | `cpw_coax.msh`, `cpw_coax_0.msh`, `cpw_lumped_0.msh`, `cpw_wave_0.msh` | `cpw_coax.msh` | REM 基线配置（`cpw.json`）所指定网格 |
| `cylinder` | `cylinder_hex.msh`, `cylinder_prism.msh`, `cylinder_tet.msh` | `cylinder_hex.msh` | `cavity_pec` 与 `driven_wave` 均使用 hex；prism 仅 `cavity_impedance` 用；tet 用于 waveguide/floquet |

---

## 执行状态追踪

### Palace Source Commit 回填记录

| 示例目录 | Palace Commit | 同步日期 | 核对人 |
|---------|---------------|---------|--------|
| adapter | TBD | — | — |
| antenna | TBD | — | — |
| coaxial | TBD | — | — |
| cpw | TBD | — | — |
| cylinder | TBD | — | — |
| transmon | TBD | — | — |

### 已完成工作

- Yew demo 所有示例 config_json 已与本目录实际 JSON 对齐（2026-04-08）
- 各 config_json 中 `L0` 格式统一为科学记数法（如 `1.0e-3`）
- `adapter` 从 Driven+LumpedPort 更正为 Eigenmode+WavePort（对齐 `hybrid.json`）
- `antenna` 材质属性从 [1] 更正为 [7]，端口配置对齐 halfwave dipole
- `cylinder` Ground 从 [4] 更正为 [3]，移除错误的 SurfaceCurrent Direction
- `coaxial` Terminal 从 [3] 更正为 [1]，材质域从 [1] 更正为 [10]
- `sbr_sphere` 频率从 1 GHz 更新为 3 GHz，RayDensity/MaxBounces/ThetaDeg 对齐

### 待完成

- [ ] 为所有 Palace 对齐示例回填 Source Commit（A2）
- [ ] 逐例输出字段差异摘要（A2）
- [ ] Transient 求解器实现（coaxial_matched/open/short）
- [ ] cylinder/floquet.json Floquet 周期边界实现验证
- [ ] transmon_amr AMR 流程集成
- [ ] MoM sphere 示例文档补充
