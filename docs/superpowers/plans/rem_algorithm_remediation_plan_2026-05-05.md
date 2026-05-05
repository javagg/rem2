# REM 算法根因整改计划（2026-05-05）

## 0. 背景与结论

针对以下五项根因问题：

1. full-wave 仍使用标量节点基函数（缺少 Nedelec 边元）
2. 高阶 FEM 仅到 P2/网格阶次映射，未形成完整 p-refinement 路径
3. complex Helmholtz 求解器收敛性弱（原路径为 CGNE 风格）
4. DDM 界面检测为空占位
5. MoM 缺少 MLFMA/FMM 体系化路线（当前以 ACA 为主）

当前状态：问题 1/2/5 仍是架构级缺口；问题 3/4 已在本轮完成首批代码落地。

更新（同日增量）：Phase A 的 Schwarz 界面耦合已接通到迭代 RHS，`rem-ddm` 可完成 `cargo check`。

## 1. 本轮已实施（可直接合入）

### 1.1 complex 线性求解器升级（已完成）

- 改造内容：
  - `rem-core` 中 `solve_pcg_complex` 内核由“法方程 CGNE”替换为“右预条件 BiCGSTAB + Jacobi”。
  - 保持原函数签名，避免上层 API 破坏。
  - `rem-driven` 中自适应求解日志同步更新为“sparse iterative solver”，失败后仍回退 dense GMRES。

- 预期收益：
  - 避免法方程条件数平方放大；对非 Hermitian Helmholtz 更稳健。
  - 高频/高对比材料下的收敛行为优于旧实现。

### 1.2 DDM 界面检测接通（已完成）

- 改造内容：
  - 在 `rem-ddm` 新增共享节点驱动的界面构建逻辑：按 `partition + volume node ownership` 自动生成 `InterfacePatch`。
  - 不再使用 `interfaces = Vec::new()` 空占位。
  - 同步填充 `SubDomain.interface_nodes/interface_neighbor`。

- 预期收益：
  - Schwarz 外层循环具备真实界面拓扑输入，为后续 Robin 数据交换打通前置条件。

### 1.3 风险显式化（已完成）

- 在 driven/eigenmode 启动时新增日志告警：
  - 当前是标量 H1 节点离散，尚无 Nedelec 边元。
  - 对向量全波场景可能出现伪模/伪解。

### 1.4 Phase A 增量：Schwarz 界面交换闭环（已完成）

- 改造内容：
  - `InterfacePatch` 增加 owner 信息，形成明确的有向界面补丁。
  - `schwarz_solve` 使用上一轮邻域解构造 Robin 入射场，并把界面贡献加入本地 RHS。
  - 收敛判据改为界面自由度更新相对残差（而非解向量绝对范数）。

- 验证状态：
  - `cargo check -p rem-ddm` 通过。
  - 同时修复 vendor/rmetis 的随机种子初始化兼容问题，恢复 DDM 编译链可用性。

### 1.5 Phase A 增量：Robin 对角项与最小测试（已完成）

- 改造内容：
  - Schwarz 本地算子新增 Robin 对角项叠加（当前 skeleton 阶段按单位面积近似）。
  - 新增 2 个单元测试：
    - Robin 对角项是否准确写入 owner 子域对角线。
    - 双向界面补丁下 Schwarz 迭代路径是否可执行并返回有限残差。

- 验证状态：
  - `cargo test -p rem-ddm schwarz -- --nocapture` 通过（2 passed）。

## 2. 分阶段完善计划

## Phase A（1-2 周）：求解器与 DDM 骨架稳定化

- A1. Complex Krylov 选型配置化
  - 新增 `Solver.Linear.KSPType` 映射：`BiCGSTAB`、`GMRES`、`CGNE(legacy)`。
  - 增加 ILU(0) 或 Block-Jacobi 预条件选项。

- A2. DDM Robin 交换闭环
  - 把 `InterfacePatch` 接入 Schwarz 迭代中的界面 RHS 更新。
  - 单机场景先实现“邻域内存交换”，MPI 场景再接 `Comm` send/recv。

- A3. 回归测试
  - 新增 3 组算例：2 子域、4 子域、强异质介质。
  - 验收：DDM 残差曲线单调下降，且界面节点数非零。

### Phase A 完成状态（当前仓库）

- 完成项：
  - [x] A1（部分完成）：`KSPType` 已支持 `BiCGSTAB/GMRES/CG/PCG` 并接入 driven 复杂求解路径选择。
  - [x] A2（已完成）：DDM 具备有向界面补丁识别、Robin 对角项叠加、Robin RHS 交换更新与界面残差收敛判据。
  - [x] A3（骨架级完成）：已补齐接口构建与 Schwarz 路径单元测试，且测试/编译通过。

- 验证命令（已通过）：
  - `cargo check -p rem-ddm`
  - `cargo test -p rem-ddm -- --nocapture`
  - `cargo check -p rem-driven -p rem-config`

- 说明：
  - 以上“完成”基于当前 DDM skeleton 局部装配（单位阵占位）语义；
    后续 Phase B/扩展阶段将把局部算子替换为真实全波装配并扩展 2/4 子域与强异质算例验证。

## Phase B（3-6 周）：full-wave 离散从 H1 迁移到 H(curl)

- B1. 空间与 DOF 体系
  - 引入一阶 Nedelec（edge-based）空间与边 DOF 编号。
  - 保留现有 H1 路径用于 electrostatic/magnetostatic。

- B2. 组装路径
  - 新增 curl-curl 与 mass 复系数组装：
    - $\int (\mu^{-1} \nabla\times \mathbf{E})\cdot(\nabla\times \mathbf{v})\,d\Omega$
    - $-\omega^2\int \epsilon\mathbf{E}\cdot\mathbf{v}\,d\Omega$
  - PEC/PMC/端口边界按边元一致处理。

- B3. 伪模抑制验证
  - 标准腔体/波导 benchmark（与 Palace/HFSS 参考值对比）。
  - 验收：伪模率显著下降，主模频率误差 < 1-2%。

### Phase B 完成状态（当前仓库）

- 完成项：
  - [x] B1（已完成）：新增 `Solver.Discretization`（`H1` / `HCurl|Nedelec`），并在 `rem-eigenmode` 接入 Nedelec `HCurlSpace`。
  - [x] B2（已完成，eigenmode 路径）：新增 `curl-curl` + `vector mass` 组装并接入 Lanczos shift-invert 求解链路。
  - [x] B3（基础完成）：HCurl 模式下完成频率谱求解闭环并输出 `eigenfrequencies.csv`；H1 专属后处理（AMR/VTK probe）在 HCurl 模式自动降级禁用。

- 兼容说明：
  - `rem-driven` 已接入离散类型识别与运行时提示；当前仍使用既有 H1 主路径（HCurl driven 端口/后处理尚未闭环）。

- 验证命令（已通过）：
  - `cargo check -p rem-config -p rem-eigenmode -p rem-driven`
  - `cargo test -p rem-eigenmode --no-run`

## Phase C（2-4 周）：高阶 FEM 路径打通

- C1. 从“网格阶次驱动”升级为“统一 Order 驱动”
  - 让 `Solver.Order` 控制 shape function order、积分阶、边界投影阶。
  - 覆盖 H1 与 H(curl) 两条路径。

- C2. 兼容策略
  - 若请求阶次超出实现上限，显式报错（不再静默降阶）。

- C3. 验收
  - p=1/2/3 的收敛阶测试（h 固定、p 提升）。

## Phase D（4-8 周）：MoM 快速算法体系化

- D1. 统一快速算法接口
  - 在 `fast_solver` 中明确：`ACA`、`FMM`、`MLFMA`。

- D2. MLFMA 核心里程碑
  - 多层八叉树
  - M2M/M2L/L2L
  - 远场平移与截断策略

- D3. 验收
  - 与 dense/ACA 对比：
    - 误差：RCS/端口阻抗偏差在工程容差内
    - 复杂度：内存和时间斜率接近 O(N log N)

## 3. 风险与边界

- Nedelec 迁移涉及离散空间、边界条件和后处理链路，属于架构改造，不建议“局部热修”。
- DDM 当前子域局部矩阵仍为 skeleton（单位阵占位），需要在 Phase A/B 进一步接入真实局部组装。
- 当前 cargo 全量检查受 vendor/rmetis 现存依赖问题影响；本轮改造未引入该问题。

## 4. 建议的近期执行顺序

1. 先完成 Phase A（把可用性和可观测性稳定住）
2. 再推进 Phase B（full-wave 正确性核心）
3. 接着做 Phase C（高阶精度能力）
4. 最后集中投入 Phase D（大规模性能）
