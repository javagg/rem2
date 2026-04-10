# 平面 MoM + FFT / MLFMM 设计文档

> 版本：v1.0，2026-04-11
> 对标：Sonnet Suite 19（平面电路 MoM）、FEKO（MLFMM）

---

## 背景与结论

### 现有 rem-mom 能力（v0.17.0）

| 能力 | 状态 |
|------|------|
| 全波 3D MoM（EFIE/MFIE/CFIE） | ✅ 已实现 |
| PMCHWT 介质目标 | ✅ 已实现 |
| RWG + 脉冲基函数 | ✅ 已实现 |
| ACA 矩阵压缩 O(N log N) | ✅ 已实现 |
| LU / GMRES / ACA+GMRES | ✅ 已实现 |
| S 参数端口扫频 | ✅ 已实现 |
| **开放边界 MoM（Beowulf/MLFMM）** | ❌ 待实现 |
| **平面 MoM + FFT（Sonnet 对标）** | ❌ 待实现 |
| **分层媒质格林函数（Sommerfeld）** | ❌ 待实现 |

---

## 方法一：平面 MoM + FFT（Spectral Domain MoM）

### 适用场景
- PCB 微带线、耦合线、不连续性
- 贴片天线、平面阵列
- 平面滤波器、平衡-不平衡变换器
- 对标：Sonnet Suite、IE3D、Momentum（ADS）

### 物理模型

```
         上方开放空间（辐射边界）
─────────────────────────────────
  金属导体条/patch（xy 平面内）  ← 离散化对象
─────────────────────────────────
  介质层 1：厚度 h₁，ε_r1，tan δ₁
─────────────────────────────────
  介质层 2：厚度 h₂，ε_r2，tan δ₂
─────────────────────────────────
  PEC 地板（或开放下边界）
```

### 核心原理

**1. 谱域格林函数**

分层媒质中，格林函数依赖 ρ = r − r'（平面平移不变性）：

```
G(r, r') = G(ρ, z, z')
```

在谱域（kx, ky）：

```
G̃(kρ, z, z') = 解析闭合公式（传输矩阵法）
```

避免了实空间 Sommerfeld 积分的数值困难。

**2. FFT 加速 O(N log N)**

利用卷积定理：Z[m,n] = G(rm − rn) → 2D FFT

```
电流 J → 2D FFT → G̃ × J̃ → 2D IFFT → 场 E
```

MVP 复杂度：O(N log N)（替代直接矩阵乘法 O(N²)）

### 实现架构（crates/planar/）

```
crates/planar/
├── Cargo.toml
└── src/
    ├── lib.rs              — 入口 run()，Problem.Type = "Planar"
    ├── layer_stack.rs      — 分层媒质描述（LayerStack）
    ├── layered_green.rs    — 谱域格林函数（传输矩阵法）
    ├── spectral_mesh.rs    — 平面网格（矩形 patch，均匀网格）
    ├── spectral_assemble.rs — FFT 加速 MVP（rustfft）
    ├── gmres.rs            — GMRES 迭代求解器（利用 FFT MVP）
    └── postprocess.rs      — S/Y/Z 参数，电流分布，辐射效率
```

### 配置示例

```json
{
  "Problem": { "Type": "Planar" },
  "Solver": {
    "Planar": {
      "FreqMin": 1e9,
      "FreqMax": 10e9,
      "FreqStep": 0.5e9,
      "Layers": [
        { "Thickness": 0.5e-3, "EpsR": 3.5, "MuR": 1.0, "LossTan": 0.002 }
      ],
      "TopOpen": true,
      "BottomPEC": true,
      "Acceleration": "FFT",
      "GmresTol": 1e-6,
      "GmresMaxIter": 500
    }
  }
}
```

### 依赖
- `rustfft`（纯 Rust FFT，支持 WASM）
- 现有 `rem-core`、`rem-config`、`rem-mesh`

---

## 方法二：MLFMM（多层快速多极子）

### 适用场景
- 电大尺寸 3D 目标（N > 10⁵ DOF）
- "Beowulf" 风格大规模 RCS 计算
- 与现有 `crates/mom` 的 ACA 互补

### 核心原理

将远场相互作用展开为多极子（球谐函数）：

```
G(r, r') ≈ Σ_l Σ_m Y_lm(r̂) · α_lm(r') / |r|^(l+1)
```

**八叉树分层（O(N log N) MVP）：**
```
上行（M2M）：子节点多极矩 → 父节点多极矩
平移（M2L）：多极矩 → 局部展开（跨分组）
下行（L2L）：父节点局部展开 → 子节点
近场：直接计算（ACA 处理）
```

### 集成位置
扩展现有 `crates/mom/src/fmm.rs`（新文件），替换 `aca_gmres_solve` 路径。

---

## 开发优先级

| 优先级 | 功能 | 工作量 | 价值 |
|--------|------|--------|------|
| 1 | **平面 MoM + FFT**（`crates/planar/`） | ~2000 行 | ⭐⭐⭐⭐⭐ Sonnet 对标 |
| 2 | **MLFMM**（扩展 `crates/mom/`） | ~3000 行 | ⭐⭐⭐⭐ 电大问题 |
| 3 | 分布式 MoM（MPI） | ~1500 行 | ⭐⭐⭐ HPC 扩展 |

---

## 参考文献

1. Harrington, R.F. "Field Computation by Moment Methods" (1993)
2. Pozar, D.M. "Input impedance and mutual coupling of rectangular microstrip antennas" (1982)
3. Bebendorf, M. "Approximation of boundary element matrices" (2000) — ACA
4. Rokhlin, V. "Rapid solution of integral equations of classical potential theory" (1985) — FMM
5. Chew et al. "Fast and Efficient Algorithms in Computational Electromagnetics" (2001) — MLFMM
