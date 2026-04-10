# FEBI (Far-field Evaluation using Boundary Integrals) 设计文档

## 结论：平面 MoM + FFT / MLFMM 加速

### 背景

在 rem2 项目中，`crates/mom` 已实现基于 PMCHWT 的矩量法（MoM）求解器，用于电磁散射问题。
当前 `crates/febi` 需要实现远场评估功能。

### 关键结论

1. **当前 mom crate 的能力**
   - 实现了 PMCHWT（Poggio-Miller-Chang-Harrington-Wu-Tsai）积分方程
   - 支持 ACA（自适应交叉近似）压缩
   - 目前仅支持三维闭合曲面（closed surface）的散射问题
   - **不支持**平面结构（planar MoM）
   - **没有** FFT 或 MLFMM 加速模块

2. **平面 MoM + FFT 加速**
   - 平面结构具有平移不变性，阻抗矩阵为 Toeplitz/Block-Toeplitz 结构
   - 可利用 FFT 将矩阵-向量乘积从 O(N²) 降至 O(N log N)
   - 适用于：印刷电路板、微带天线、FSS（频率选择表面）等平面结构

3. **MLFMM（多层快速多极子方法）**
   - 适用于一般三维结构，非平面限定
   - 将矩阵-向量乘积降至 O(N log N) 或 O(N)
   - 比 ACA 更系统化，适合大规模问题

4. **FEBI crate 的定位**
   - FEBI = Far-field Evaluation using Boundary Integrals（基于边界积分的远场评估）
   - 利用 mom crate 的近场解，计算远场方向图
   - 支持 RCS（雷达散射截面）、辐射方向图等远场量的计算

### 实现方案

#### Phase 1：基础远场评估（当前任务）

```
crates/febi/
├── Cargo.toml
└── src/
    ├── lib.rs          # 公共 API
    ├── far_field.rs    # 远场积分计算
    ├── rcs.rs          # RCS 计算
    └── radiation.rs    # 辐射方向图
```

**核心公式**：
远场电场由等效电流/磁流积分给出：

```
E_far(r̂) = -jk/(4π) * exp(-jkr)/r * [Z₀(r̂ × (r̂ × N)) - L]
```

其中：
- `N` = 电流辐射积分（由表面电流 J 得到）
- `L` = 磁流辐射积分（由表面磁流 M 得到）

#### Phase 2（未来扩展）

- 平面 MoM + FFT 加速（用于平面结构的大规模计算）
- MLFMM 加速（用于一般三维大规模结构）

### 与其他 crate 的依赖关系

```
rem-febi
  ├── rem-mom       (获取表面电流/磁流解)
  ├── rem-mesh      (获取网格信息)
  └── rem-config    (配置参数)
```

### 配置结构

在 `config::schema` 中已预留 `febi: Option<FebiConfig>` 字段（当前为 None）。

`FebiConfig` 需要包含：
- 远场观测点/方向设置
- 输出格式（RCS、方向图等）
- 频率范围

---

*文档创建时间：2026-04-10*
*状态：设计阶段*
