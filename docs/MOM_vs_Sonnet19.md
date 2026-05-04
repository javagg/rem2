# REM MoM 求解器 vs Sonnet Suite 19 对比分析

> 版本：2026-05-04（更新至 REM v0.22.0）
> REM 基准版本：v0.22.0（`crates/mom/`）  
> Sonnet 基准版本：Suite 19（商业授权，Sonnet Software Inc.）

---

## 1. 技术定位与适用场景

| 维度 | REM MoM | Sonnet Suite 19 |
|------|---------|-----------------|
| **定位** | 通用全波三维表面积分 MoM，适用于任意封闭/开放 3D 散射体 | 平面化三维（2.5D）MoM，专用于多层基板上的平面导体结构 |
| **几何适用** | 任意三维封闭曲面（球体、飞机、天线等任意曲面网格） | 平面导体层叠（MMIC、PCB、微带、槽线、天线阵列平面结构） |
| **目标用户** | 雷达散射截面（RCS）计算、3D 目标散射与辐射 | 微波/毫米波电路 S 参数提取、无源器件建模、RFIC/MMIC 版图验证 |
| **参考工具** | 对标 FEKO（表面积分方程）、OpenEMS MoM、WIPL-D | 商业"标准"，常作为 ADS、Cadence AWR 的 EM 联合仿真后端 |

---

## 2. 核心算法对比

### 2.1 积分方程公式

| 特性 | REM MoM | Sonnet Suite 19 |
|------|---------|-----------------|
| **电场积分方程（EFIE）** | ✅ RWG 矢量基，标量 Helmholtz 核 | ✅ Rao-Wilton-Glisson（RWG）平面截断版 |
| **磁场积分方程（MFIE）** | ✅ 身份项 + curl-Green 项，完整实现 | 内部使用，不直接暴露 |
| **CFIE（EFIE+MFIE 混合）** | ✅ α 参数可配（`α=0.5` 默认），消除内谐振 | Sonnet 对封闭结构默认采用无内谐振公式 |
| **PMCHWT 介质目标** | ✅ 完整 2N×2N 块矩阵，J+M 未知量，ε_r/μ_r 可配 | ❌ Sonnet 不支持任意三维介质目标散射；仅支持平面多层介质基板 |
| **基函数类型** | RWG 矢量基 + 脉冲（pulse）标量基，可切换 | RWG 变种（平面版），固定正交矩形/三角形网格 |
| **Green 函数** | 均匀自由空间 3D Helmholtz Green 函数 | **分层介质 Green 函数**（Sommerfeld 积分），精确建模多层基板 |

### 2.2 奇异积分处理

| 方法 | REM MoM | Sonnet Suite 19 |
|------|---------|-----------------|
| **自积分（对角块）** | Duffy 变量变换（自积分），4 阶 | 解析预计算（基于矩形单元对称性） |
| **近奇异（共边/共顶点）** | Sauter-Schwab 奇异积分，4 阶 | Sonnet 特有的矩形元素近场修正 |
| **Gauss 求积** | Dunavant 三角形高斯求积，阶次 1/3/5/7/9 可选 | 平面矩形/三角 Gauss，内部自适应选阶 |

### 2.3 快速算法与矩阵压缩

| 加速方法 | REM MoM | Sonnet Suite 19 |
|---------|---------|-----------------|
| **ACA（自适应截面近似）** | ✅ 部分主元 ACA，Z≈U·V^T，O(N·r) 矩阵向量积；`FastSolver: "ACA"` | ❌ 不支持 ACA |
| **FFT 加速 MoM** | ❌ 未实现（待 FMM 路线图） | ✅ **核心优势**：利用平面周期性，FFT 加速矩阵填充和矩阵-向量积，O(N log N) |
| **FMM（快速多极子）** | ❌ 配置项已预留，运行时返回错误 | ❌ Sonnet 不使用 FMM |
| **直接密集 LU** | ✅ nalgebra dense LU，O(N³) | ✅ 内置 dense LU（小问题） |
| **GMRES** | ✅ 重启 GMRES（restart=30，tol=1e-8），O(N²·restart) | ✅ Krylov 迭代求解器 |

### 2.4 线性系统规模

| 问题规模（RWG 基函数数 N）| REM 推荐求解器 | Sonnet 能力 |
|--------------------------|--------------|------------|
| N < 500 | Dense LU（精度最高） | Dense LU |
| 500 ≤ N < 3000 | GMRES | FFT 加速迭代 |
| N ≥ 3000 | **ACA + GMRES**（O(N·r)） | **FFT 加速**（O(N log N)），可处理 N > 100,000 |
| N > 50,000 | 尚未验证（FMM 路线图中） | Sonnet 商业版支持，需多线程或多核机器 |

---

## 3. 几何与网格能力

| 能力 | REM MoM | Sonnet Suite 19 |
|------|---------|-----------------|
| **几何维度** | 任意三维封闭/开放曲面（`.msh` Gmsh 导入） | 平面分层（2.5D）：导体必须在某一层平面上 |
| **网格类型** | 非结构三角网格（Tri3），任意曲率 | 矩形网格为主（自动剖分），也支持 conformal 三角剖分（Suite 高级版） |
| **网格导入** | `.msh` v2/v4（GMSH 完整支持），物理组映射 | 专有 `.son` 格式，支持 DXF/Gerber/ODB++ 导入 |
| **自适应网格** | ✅ v0.20.0 MoM AMR（Dörfler 标记 + 1→4 细分） | ✅ 自适应 Sonnet 单元细化（基于 S 参数收敛） |
| **曲面建模** | 任意曲率（球体、复杂天线、飞机翼面） | 仅平面或近似分段平面 |

---

## 4. 激励与后处理

### 4.1 激励类型

| 激励 | REM MoM | Sonnet Suite 19 |
|------|---------|-----------------|
| **平面波（入射场）** | ✅ θ/φ 极化，任意入射角，频率扫描 | ❌ 不支持平面波激励（Sonnet 专注于端口激励） |
| **集总端口激励** | ✅ v0.17.0 集总端口，内阻可配 | ✅ 集总端口，内阻可配 |
| **波导端口** | ✅ v0.22.0 WavePort 图 Laplacian 模式加权激励（`Type:"WavePort"`, `Mode:N`） | ✅ 矩形/同轴波导端口，去嵌入 |
| **差分端口（混合模）** | ✅ v0.22.0 `PairWith` 字段 + `single_ended_to_mixed_mode` → 2×2 Sdd/Scc 及 2N×2N 全混合模矩阵 | ✅ |
| **自动去嵌入（Deembedding）** | ✅ v0.22.0 参考面相位+衰减去嵌入（`DeembedLength` / `DeembedEpsEff` / `DeembedAlpha`） | ✅ 港口参考面精确去嵌入 |

### 4.2 后处理

| 输出 | REM MoM | Sonnet Suite 19 |
|------|---------|-----------------|
| **RCS 方向图** | ✅ θ/φ 扫描，`rcs.csv`（θ, φ, σ dBsm），全球面 | ❌ 不支持 RCS（非散射工具） |
| **S 参数** | ✅ v0.17.0 集总端口 + Touchstone `.sNp` | ✅ 完整 S/Y/Z 参数矩阵，Touchstone `.sNp` 导出 |
| **Z/Y 矩阵** | ✅ v0.21.0 port-Z.csv / port-Y.csv | ✅ |
| **近场（任意点）** | ✅ v0.22.0 `probe_e_field_portN.csv`（全端口 × 全频率 RWG 辐射积分） | ✅ 2D/3D 近场可视化 |
| **表面电流 VTK** | ✅ v0.22.0 `surface_current_portN_*.vtk`（RWG 矢量电流，J_real/J_imag/J_mag） | ✅ 电流密度 2D 可视化（专有格式） |
| **远场辐射方向图（端口激励）** | ✅ v0.22.0 `far_field_portN.csv`（N_θ/N_φ/方向性 dBi，配置 `Solver.FarField`） | ✅ 平面天线方向图后处理 |
| **传输线参数（R/L/G/C）** | ✅ v0.21.0 `tline_params.csv`（ABCD→RLGC，2-port） | ✅ 传输线参数提取 |
| **等效电路综合** | ✅（Driven FEM：VF 极点-留数 + SPICE `.cir`） | ✅（Circuit Element 模型，导出 SPICE） |

---

## 5. 材料与物理建模

| 特性 | REM MoM | Sonnet Suite 19 |
|------|---------|-----------------|
| **PEC 导体** | ✅ CFIE 公式，完美导体 | ✅ 理想 PEC 或有限电导率导体 |
| **介质目标（均匀）** | ✅ PMCHWT（2N×2N）ε_r/μ_r 可配 | ❌ 不支持任意三维介质散射 |
| **分层介质基板** | ✅ v0.18.0 Sommerfeld/DCIM 分层 Green 函数 | ✅ **核心优势**：Sommerfeld 积分精确建模任意层叠基板 |
| **有损导体（表面阻抗）** | ✅ v0.19.0 Leontovich SIBC（`WallConductivity` 配置） | ✅ 表面电阻（σ 有限），R_s 建模 |
| **各向异性介质** | ❌ MoM 无；FEM 求解器支持 3×3 张量 | ❌ 各向同性基板 |
| **频变材料** | ❌ MoM 中未实现 | ✅ 有限，通过宽带建模（Debye/Lorentz 近似） |

---

## 6. 平台与工程集成

| 特性 | REM MoM | Sonnet Suite 19 |
|------|---------|-----------------|
| **编程语言** | 纯 Rust（零 C++ 依赖） | C++（专有），闭源 |
| **授权** | 开源（MIT/Apache 协议） | 商业授权（年费/节点锁定） |
| **WASM 支持** | ✅ 编译至 `wasm32-unknown-unknown`，浏览器运行 | ❌ 仅本地 Windows/Linux |
| **并行计算** | ✅ Rayon（CPU 多线程）；MPI feature 开关 | ✅ 多线程 + 可选多节点并行（HPC 版） |
| **EDA 工具集成** | ❌（独立工具） | ✅ ADS、AWR、Cadence、Keysight 无缝集成 |
| **参数化扫描 / 优化** | ❌ | ✅ 内置参数扫描 + 梯度优化 |
| **脚本接口** | ✅ JSON/YAML 配置 + Rust API | ✅ Python API（sonnet.exe 命令行） |
| **版本控制友好** | ✅ 纯文本配置 + Git 可跟踪 | JSON-like `.son` 文件，部分版本可比较 |
| **跨平台** | ✅ Windows/Linux/macOS/WASM | ✅ Windows + Linux |

---

## 7. 验证精度对比

| 验证场景 | REM MoM | Sonnet Suite 19 |
|---------|---------|-----------------|
| **PEC 球体 RCS（Mie）** | ka≈10.5，误差 < 0.5 dB（CFIE + RWG） | 不适用（Sonnet 无散射场景） |
| **微带传输线 S11/S21** | 端口已实现（v0.17.0）；精度待与 Sonnet 基准对比 | < 0.1 dB，工业标准 |
| **传输线 RLGC** | v0.21.0 ABCD→RLGC 提取（待实测验证） | 工业标准 |
| **贴片天线 S11** | 端口已实现；AMR 已实现；精度待验证 | < 0.5 dB，业界基准 |
| **二维单站 RCS** | < 0.5 dB vs Mie | 不适用 |
| **PMCHWT 介质球** | 已通过组装有限性验证，精度测试进行中 | 不适用 |

---

## 8. 综合差距与优势分析

### 8.1 REM MoM 相对 Sonnet 19 的**优势**

1. **任意三维曲面目标**：Sonnet 只能处理平面结构，REM 可仿真球体、飞机翼面、导弹头锥等任意三维封闭曲面，是雷达 RCS 计算的核心能力。

2. **PMCHWT 三维介质目标**：支持均匀介质球/任意均匀介质体的散射分析（J+M 双未知量 2N×2N 系统），Sonnet 无此能力。

3. **完整散射场后处理**：提供双/单站 RCS 方向图（dBsm），全球面 θ/φ 扫描，直接输出 `rcs.csv`。Sonnet 没有 RCS 输出。

4. **平面波激励**：支持任意方向入射平面波（θ_inc, φ_inc 可配），适合雷达截面计算；Sonnet 仅支持端口激励。

5. **CFIE 内谐振抑制**：α 参数可调的 CFIE 消除 EFIE 在谐振频率的内谐振问题；Sonnet 利用平面结构物理特性天然避免了此问题。

6. **开源免费 + WASM**：MIT/Apache 协议，可嵌入任意系统；支持浏览器端运行，Sonnet 为昂贵商业软件。

7. **ACA 矩阵压缩**：对大规模三维散射问题，ACA+GMRES 可将内存从 O(N²) 降至 O(N·r)；Sonnet 无 ACA。

8. **Rust 内存安全**：零 C++ 依赖，无内存越界/悬空指针风险；Sonnet 为 C++ 实现。

---

### 8.2 Sonnet 19 相对 REM MoM 的**优势**

1. **分层介质 Green 函数**：Sonnet 使用精确的 Sommerfeld 积分 Green 函数，天然建模多层基板（FR4、Rogers、LTCC 等叠层），无需将基板体积离散化。这是平面电路分析的物理根基——REM 使用自由空间 Green 函数，无法处理嵌入基板的导体。

2. **FFT 加速 MoM（O(N log N)）**：利用平面结构的移位不变性，用 FFT 加速矩阵填充和矩阵-向量积，可处理 N > 100,000 的超大平面电路。REM ACA 虽有效，但对平面问题效率不及 FFT。

3. ~~**波导端口与去嵌入**：Sonnet 波导端口（矩形/同轴）和端口参考面精确去嵌入仍是 REM 的差距~~（✅ v0.22.0 已实现 WavePort 模式加权激励、`PairWith` 差分混合模、`DeembedLength` 参考面去嵌入）

4. **EDA 生态集成**：与 ADS、Cadence、AWR 无缝联动，支持从版图直接驱动 EM 仿真，参数化扫描和优化闭环；REM 当前为独立工具，无 EDA 集成。

5. **工程成熟度与支持**：Sonnet 商业产品有 40 年历史，大量工程师验证案例、技术支持、培训体系；REM MoM 为 v0.22.0 活跃开发版本，主要 EM 功能已覆盖。

6. ~~**自适应网格细化（MoM）**：Sonnet 在 MoM 内部支持自适应单元细化；REM MoM 暂无 AMR~~（✅ v0.20.0 已实现 Dörfler AMR）

7. ~~**有损导体建模**：Sonnet 支持有限电导率导体（表面阻抗边界条件）；REM MoM 仅支持 PEC~~（✅ v0.19.0 已实现 SIBC/WallConductivity）

8. ~~**远场辐射方向图（端口激励）**：Sonnet 支持平面天线远场方向图~~（✅ v0.22.0 已实现端口激励 RWG 远场辐射积分，输出 N_θ/N_φ/方向性 dBi）

---

## 9. 场景选型建议

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

## 10. REM MoM 差距路线图

基于上述对比，弥补与 Sonnet 差距的关键优先项：

| 优先级 | 缺失功能 | 说明 |
|--------|---------|------|
| **高** | **分层介质 Green 函数** | 实现 Sommerfeld 积分（数值/离散复像法），支持 PCB/MMIC 多层基板 |
| **高** | **端口激励 + S 参数** | 集总端口 + 波导端口，输出 Touchstone `.s2p`，打通 EDA 接口 |
| **中** | **FFT 加速矩阵填充** | 针对平面结构利用 FFT，将复杂度从 O(N²) 降至 O(N log N) |
| **中** | **有损导体（表面阻抗）** | 有限电导率 σ 的 SIBC（表面阻抗边界条件） |
| **低** | **FMM（快速多极子）** | 通用三维大规模散射加速（已在配置项预留） |
| **低** | **MoM 内 AMR** | 基于远场/近场误差指示的自适应三角网格细化 |

---

## 参考资料

- REM `crates/mom/src/` 源代码（`assemble.rs`, `aca.rs`, `pmchwt.rs`, `singular.rs`, `green.rs`）
- Rao, Wilton, Glisson, "Electromagnetic Scattering by Surfaces of Arbitrary Shape," IEEE TAP, 1982
- Peterson, Ray, Mittra, *Computational Methods for Electromagnetics*, IEEE Press, 1998
- Sonnet Software, *Sonnet User's Guide*, Suite 19, 2024
- Bebendorf, "Approximation of boundary element matrices," *Numerische Mathematik*, 2000
- [Sonnet Software 官网](https://www.sonnetsoftware.com)
- [Wikipedia: Sonnet Software](https://en.wikipedia.org/wiki/Sonnet_Software)
