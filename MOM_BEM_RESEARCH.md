# MoM/BEM 技术方案调研报告

> 项目：rem2 电磁仿真工具  
> 日期：2026-04-05  
> 目标：评估在现有 FEM 框架上实现矩量法（MoM）和边界元法（BEM）的技术方案

---

## 1. 项目现状

| 指标 | 说明 |
|------|------|
| 代码量 | ~6,800 行（14 个 crate） |
| 当前版本 | v0.2 |
| 已实现求解器 | FEM 静电、Eigenmode、Driven 频域 |
| 配置兼容性 | Palace 配置格式兼容 |
| 运行目标 | 原生 + WASM |

---

## 2. 可复用的已有资产

### 2.1 核心基础设施

| 组件 | 文件 | 说明 |
|------|------|------|
| 稀疏矩阵（COO/CSR） | `crates/core/src/sparse.rs` | 装配与 PCG 求解，约 428 行 |
| FEM 装配框架 | `crates/electrostatic/src/assemble.rs` | P1 单元装配，可直接改造为 BEM 表面积分 |
| 边界条件框架 | `crates/electrostatic/src/bc.rs` | 边界条件应用接口 |
| 网格数据结构 | `crates/mesh/src/mesh_data.rs` | `RemMesh` + `boundary_elements` + 物理 tag 映射 |
| GMSH 解析器 | `crates/mesh/src/gmsh.rs` | `.msh` 文件加载，约 635 行 |
| 材料管理 | `crates/materials/src/material.rs` | 介电/磁参数管理 |
| Palace 配置 | `crates/config/src/schema.rs` | Palace 兼容配置格式 |
| MPI 抽象 | `crates/parallel/src/` | 并行通信接口 |

### 2.2 MoM/BEM 直接可用的部分

- `RemMesh` 已有 `boundary_elements` 字段，可直接用于 BEM 表面网格提取
- FEM 装配循环（单元遍历 + 高斯积分调用点）结构可参照改写为 BEM 装配
- 现有 PCG 求解器可用于稀疏化处理后的 BEM 矩阵（配合快速算法）

---

## 3. 必须新增的组件

### P0 — 阻塞项（必须先完成）

| 组件 | 工作量 | 说明 |
|------|--------|------|
| 复数矩阵支持 | 小 | `num-complex` 或自定义 `Complex<f64>`，扩展现有稀疏/密集矩阵 |
| 密集矩阵库 | 中 | `faer` 或 `ndarray` + BLAS 绑定，BEM 阻抗矩阵是稠密的 |
| LU 分解求解器 | 小 | 直接调用 `faer::solvers::Lu` 即可 |

### P1 — 核心功能

| 组件 | 工作量 | 技术难度 | 说明 |
|------|--------|----------|------|
| 3D 标量 Green 函数 | 中 | 低 | $G = e^{-jkr}/(4\pi r)$，含奇异点处理 |
| 2D Green 函数 | 中 | 低 | 第二类 Hankel 函数 $H_0^{(2)}(kr)$ |
| 表面高斯求积 | 中 | 中 | 三角面元上的数值积分，7/12/16 点规则 |
| **奇异积分处理** | 大 | **高** | 自积分/近奇异：Duffy 变换、Sauter-Schwab 规则 |
| RWG 基函数 | 中 | 中 | Rao-Wilton-Glisson 矢量基函数，需边-面拓扑 |
| 边-面拓扑构建 | 中 | 中 | 从 `RemMesh` 构建共享边数据结构 |

### P2 — 后处理与优化

| 组件 | 工作量 | 说明 |
|------|--------|------|
| 远场/RCS 计算 | 中 | 等效电流积分，输出方向图 |
| 快速算法（FMM/ACA） | 大 | 降低 $O(N^2)$ 复杂度，大规模问题必需 |
| MFIE/CFIE 方程 | 中 | 多方程组合提升条件数 |

---

## 4. 技术架构设计

### 4.1 新增 crate 布局

```
crates/
  mom/                     # 矩量法主 crate
    src/
      lib.rs
      green.rs             # Green 函数（标量/矢量）
      quadrature.rs        # 表面高斯求积规则
      singular.rs          # 奇异积分（Duffy/Sauter-Schwab）
      basis/
        rwg.rs             # RWG 基函数
        pulse.rs           # 脉冲基函数（简单情形）
      assemble.rs          # 阻抗矩阵装配
      solver.rs            # 密集矩阵 LU 求解
      postprocess.rs       # 电流分布、远场、RCS
```

### 4.2 核心数据流

```
RemMesh
  └─ 提取表面三角网格
       └─ 构建边-面拓扑
            └─ 生成 RWG 基函数
                 └─ 装配阻抗矩阵 Z[m,n] (密集复数矩阵)
                      ├─ 奇异对角块：Duffy 变换
                      └─ 远场块：标准高斯积分
                           └─ 施加激励向量 V
                                └─ LU 求解 ZI = V
                                     └─ 后处理：RCS / 近场
```

### 4.3 方程选择

| 场景 | 推荐方程 | 原因 |
|------|----------|------|
| PEC 散射（外问题） | EFIE | 实现最简单，适合入门 |
| PEC 散射（大电体） | CFIE = α·EFIE + (1-α)·MFIE | 避免内谐振，条件数好 |
| 介质散射 | PMCHWT | 同时满足切向 E/H 连续条件 |
| 静电 BEM | Laplace EFIE | 频率为零的退化情形 |

---

## 5. 关键技术难点

### 5.1 奇异积分（最高风险）

BEM 阻抗矩阵的奇异性来自 Green 函数在 $r \to 0$ 时的发散：

- **自积分**（$m = n$）：使用 Duffy 变换去奇异化
- **近奇异积分**（邻近单元）：Sauter-Schwab 规则（4 种几何接触情形）
- **参考实现**：可参考 BEM++ / bempp-cl 的 Python 实现移植到 Rust

### 5.2 复数稀疏 vs 密集矩阵

- BEM 产生 $N \times N$ **满**复数矩阵，$N$ 为未知量数
- 1000 个三角面元 → $N \approx 1500$ RWG → 矩阵约 18 MB（可接受）
- 超过 $N > 10000$ 需引入 FMM 或 ACA 压缩

### 5.3 WASM 兼容性

- `faer` 支持 WASM（纯 Rust，无 BLAS 依赖选项）
- 奇异积分计算密集，WASM 性能约为原生的 50-70%，需评估

---

## 6. 依赖选型建议

| 用途 | 推荐 crate | 备选 |
|------|-----------|------|
| 复数类型 | `num-complex` | 已有 |
| 密集矩阵 + LU | `faer` | `nalgebra` |
| 特殊函数（Hankel） | `num-complex` + 手写 | `special` |
| 并行装配 | `rayon` | 已有 |
| 迭代求解（大规模） | `faer` GMRES | 手写 |

---

## 7. 实现路线图

### Phase 1（2-3 周）：基础设施

- [ ] 创建 `crates/mom` crate，加入 workspace
- [ ] 复数矩阵：扩展 `core` 或引入 `faer`
- [ ] 表面网格提取：从 `RemMesh.boundary_elements` 构建三角面元列表
- [ ] 基础高斯求积规则（三角面元 7 点）

### Phase 2（3-4 周）：Green 函数 + 积分

- [ ] 3D 标量 Green 函数实现
- [ ] 非奇异远场积分
- [ ] 脉冲基函数（EFIE 标量版，用于验证）
- [ ] 简单平板散射验证（与解析解对比）

### Phase 3（4-5 周）：RWG + 奇异积分

- [ ] 边-面拓扑构建
- [ ] RWG 基函数生成
- [ ] Duffy 变换（自积分）
- [ ] Sauter-Schwab 规则（4 种近奇异情形）
- [ ] 完整 EFIE 阻抗矩阵装配

### Phase 4（2-3 周）：求解与验证

- [ ] LU 求解 + 电流分布输出
- [ ] 远场积分 + RCS 计算
- [ ] MVP 验证：PEC 球体 @ 1 GHz（vs FEKO/Mie 解析解，误差 < 5%）
- [ ] 集成到 Palace 配置格式

---

## 8. 最小可行产品（MVP）定义

**目标问题**：单个 PEC 球体平面波散射，$f = 1$ GHz

**验收标准**：
- 编译通过（Rust release 模式）
- RCS 与 Mie 级数解析解误差 < 5%
- 支持 ≥ 10×10 表面网格
- 单线程运行时间 < 30 秒

---

## 9. 可行性评估

**总体可行性：4/5**

**优势**：
- 现有 FEM 框架提供了完整的网格、装配、求解基础设施
- Rust 生态（faer、rayon）对 BEM 数值计算友好
- WASM 目标不构成根本性障碍

**主要风险**：
- 奇异积分实现复杂度高，是唯一技术难关
- BEM 全矩阵与现有稀疏矩阵路径不同，需新增密集矩阵路径
- RWG 基函数需要仔细处理边的方向约定

**建议**：先实现脉冲基函数的标量 EFIE（最简单情形）验证整体流程，再逐步升级到 RWG + CFIE。

---

## 10. 参考资料

- Gibson, W. C. *The Method of Moments in Electromagnetics* (2021)
- Rao, Wilton, Glisson. *Electromagnetic scattering by surfaces of arbitrary shape*, IEEE TAP 1982
- Sauter & Schwab. *Boundary Element Methods* (2011) — 奇异积分章节
- [BEM++ 开源库](https://bempp.com) — Python 参考实现
- [OpenEMS](https://openems.de) — 对比工具
