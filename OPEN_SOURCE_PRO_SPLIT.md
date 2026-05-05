# REM 开源版 vs. 专业版分割方案

> 版本：v1.0，2026-05-05
> 适用项目：REM v0.17+

---

## 一、分割逻辑原则

| 纳入开源 | 纳入专业版 |
|---------|---------|
| 学术界/研究者常用功能 | 工业界付费对标功能 |
| 已有成熟开源替代品的功能 | REM 独有差异化功能 |
| 吸引贡献者、建立生态的基础层 | 商业可变现的核心竞争力 |
| 没有护城河的基础设施 | 开发成本高、替代难的算法 |

---

## 二、具体 Crate 分配

### 开源版（Community，MIT/Apache）

| Crate | 纳入原因 |
|-------|---------|
| `core` / `config` / `mesh` / `bc` / `materials` / `result` | 基础设施，必须开放才能有生态 |
| `electrostatic` / `magnetostatic` | 学术竞品（Elmer、FEniCS）均免费，护城河低 |
| `eigenmode` | 谐振腔分析，学术用途为主，竞品多 |
| `driven`（基础单端口，无 ROM） | Palace 本身开源，必须与其保持可比性 |
| `transient`（基础波形） | 同上 |
| `touchstone`（读写，无矩阵转换） | 格式标准，开放有利于生态集成 |
| `convert` | 工具性转换，开放有利于第三方集成 |
| `wasm` / `yew-app` | WASM 浏览器运行是最大差异化宣传点，必须开放展示 |
| `cli`（基础命令） | 工具入口，需要开放让用户能直接使用 |
| `bem`（Laplace P0） | 实现较基础，学术用途为主 |

### 专业版（Pro/Enterprise，私有）

| Crate / 模块 | 商业价值原因 |
|-------------|------------|
| `planar` | Sonnet 最直接竞争对手，平面 MoM 核心商业价值（Sonnet 年费 $10k–30k/节点）|
| `mom`（CFIE/PMCHWT/ACA） | 三维全波 MoM，工业界高度商业化（HFSS 竞品功能）|
| `layered_green` | Sommerfeld 积分分层介质 Green 函数，Sonnet 核心算法 |
| `febi` | FE-BI 混合求解器，高端前沿功能 |
| `ddm` | 区域分解法，大规模并行高端需求 |
| `sbr` | SBR+ 射线追踪，HFSS 收费功能直接对标 |
| `parallel`（MPI） | 企业级 HPC 并行需求 |
| `optim` | Nelder-Mead + Monte Carlo 良率分析，EDA 优化是核心商业功能 |
| `driven/sparams_analysis.rs` | 群延迟、Rollett K 因子、MSG/MAG、无条件稳定性分析 |
| `touchstone/matrix_convert.rs` | Z/Y/ABCD/T 矩阵转换、N 端口级联，工业 RF 设计必备 |
| `driven`（ROM 快速频率扫描） | ROM 是速度差异最大点，高价值功能 |

---

## 三、工程实现方案

### 方案 A：Cargo Feature Flags（推荐早期）

在 `crates/cli/Cargo.toml` 中定义：

```toml
[features]
default = []
pro = [
    "rem-planar",
    "rem-mom",
    "rem-febi",
    "rem-ddm",
    "rem-sbr",
    "rem-parallel",
    "rem-optim",
    "rem-driven/rom",
    "rem-driven/sparams_analysis",
    "rem-touchstone/matrix_convert",
]
```

### 方案 B：双仓库架构（推荐中期）

```
github.com/<org>/rem              ← 开源仓库（MIT/Apache），发布到 crates.io
github.com/<org>/rem-pro          ← 私有仓库，依赖开源 crate 作为 git path 依赖
  crates/
    planar/       ← 平面 MoM（含 Sommerfeld）
    mom-full/     ← 完整 MoM（含 MLFMM 扩展）
    optim-pro/    ← 完整优化框架
    sparams-pro/  ← S 参数高级分析
```

Pro 二进制编译：`cargo build --release --features pro`

### 许可证验证方式（三选一）

| 方式 | 适用阶段 | 实现复杂度 |
|------|---------|-----------|
| **私有 Git 仓库访问控制** | 早期（v1.0）| 极低——Cargo token 控制拉取权限 |
| **离线许可证文件** | 中期（v2.0）| 中——`~/.rem/license.key` + HMAC 签名，无需联网 |
| **联网激活** | 成熟期 | 高——SaaS 场景 |

---

## 四、仓库结构（推荐）

```
rem/                              ← 开源 workspace
  crates/
    core/ config/ mesh/ bc/ materials/ result/
    electrostatic/ magnetostatic/ eigenmode/
    driven/          (无 ROM、无高级 S 参数分析)
    transient/       (基础波形)
    touchstone/      (读写，无矩阵转换)
    bem/
    convert/
    wasm/ yew-app/
    cli/             (--features community)

rem-pro/                          ← 私有 workspace（依赖 rem 开源 crate）
  crates/
    planar/
    mom/
    layered_green/
    febi/
    ddm/
    sbr/
    parallel/
    optim/
    driven-pro/      (ROM + 高级 S 参数分析)
    touchstone-pro/  (矩阵转换)
    cli-pro/         (完整 CLI，依赖上述全部)
```

---

## 五、定价参考

| 版本 | 目标用户 | 核心功能 | 参考定价 |
|------|---------|---------|---------|
| **Community** | 高校/研究者/个人 | 开源全部功能 | 免费 |
| **Pro Individual** | 独立 RF/EM 工程师 | + Planar MoM、S 参数高级分析、矩阵转换 | $500–1,500/年 |
| **Pro Studio** | 小团队（≤5 席） | + MoM/FE-BI/DDM + 优化 + ROM | $3,000–8,000/年/席 |
| **Enterprise** | 公司/研究所 | 全功能 + MPI + 优先技术支持 | 议价 |

---

## 六、优先保护的功能

> **最关键**：`planar`（平面 MoM）是最值得优先保护的专业版功能。
>
> - Sonnet Suite 年费约 $10k–30k/节点
> - REM `planar` crate 已实现 Sommerfeld 积分分层介质 Green 函数 + 2D FFT 卷积加速，直接对标 Sonnet 核心算法
> - 这是最明确、最可量化的商业变现切入点

次优先：`optim`（良率分析）+ `driven` ROM 快速扫描 + S 参数高级分析。

---

## 七、开源策略建议

1. **先全部开源，再逐步商业化**：v0.x 期间保持全开源，积累用户口碑和 GitHub Star
2. **v1.0 时分割**：在有稳定用户基础后再引入 Pro 层，避免过早商业化吓跑社区
3. **学术永久免费**：向学术机构提供免费 Pro 许可证，换取论文引用和背书
4. **开源核心，闭源加速**：开源基础算法（如 MoM RWG），闭源 ACA 压缩和 MLFMM 加速实现
