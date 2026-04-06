# REM Yew Demo & Examples 改造要求

本文档记录 yew demo 与 examples 目录的后续改造要求，作为统一实现与验收依据。

## 范围
- 目录范围：`examples/` 与 `crates/yew-app/`。
- 目标：让示例与 Palace 保持同步、可在浏览器中稳定运行（含 MPI 模式）、并具备可用的输出浏览体验。

## 要求清单

### 1. 与 Palace 示例严格对齐
- Palace 项目中已纳入的对应示例，需在本仓库同步以下内容：
  - `.msh` 网格文件（命名、版本、几何实体标签保持一致）。
  - 配置文件字段（Problem/Model/Domains/Boundaries/Solver/Postprocessing）结构与关键参数一致。
- 同步策略：
  - 增加“来源与版本”注记（来源仓库、提交号、同步日期）。
  - 每次 Palace 示例更新后，需进行差异核对并更新本仓库。
- 验收标准：
  - `examples/*/*.json` 与 Palace 对应示例字段可逐项对照。
  - 缺失网格或字段漂移需在 PR 中明确说明。

### 2. 构建期内嵌示例 msh，运行时加载
- 运行时不得依赖网络下载示例网格文件。
- 所需 `.msh` 在构建阶段打包进 wasm 产物（或通过 Rust `include_bytes!` 固化）。
- 前端运行时从内嵌资源读取并传入求解器。
- 验收标准：
  - 断网状态下仍可运行示例。
  - 发布产物中可追溯所有示例网格来源。

### 3. MPI 模式增加多 rank 输出信息
- demo 在 MPI 模式下应提供更丰富的每 rank 运行信息，至少包含：
  - rank 启动/结束时间点。
  - 当前阶段信息（例如初始化、求解、barrier、收敛/失败）。
  - 关键统计（如迭代次数、局部误差/残差摘要，若可得）。
- 输出要求：
  - 全局日志保留系统事件。
  - per-rank 面板以 rank 视角展示业务日志。
- 验收标准：
  - 2/4/7 ranks 下均能看到每个 rank 的独立输出。

### 4. 使用 OPFS 处理示例文件输出，并支持弹窗浏览
- 前端使用 OPFS（Origin Private File System）承接示例输出文件。
- 运行完成后，允许用户通过弹窗查看输出目录与文件内容（至少支持文本/CSV 预览）。
- 最小能力建议：
  - 输出文件列表。
  - 文件大小与更新时间。
  - 点击预览 + 下载。
- 验收标准：
  - 浏览器刷新后，输出文件在同一 origin 下可继续访问（在浏览器策略允许范围内）。

当前实现进度：
- 已落地最小版本：串行与 MPI 运行结束后，将日志、配置和结果摘要写入 OPFS。
- 已支持弹窗列出文件并预览文本内容。
- Serial 模式已导出 `phi.csv`，并在可用时导出 `e_field.csv` / `b_field.csv`。
- MPI 模式已生成 `ranks_summary.csv`（rank 行数、phase 覆盖、最后状态）。
- 下载能力以文本文件为第一版范围；后续二进制产物可单独扩展。

### 5. 控制面板瘦身 + rank 面板 1 行 4 列
- 示例选择控制区当前占用过大，需要紧凑化：
  - 压缩垂直间距与组件高度。
  - 将次要控件折叠到“高级设置”。
- rank 输出布局调整：
  - 默认采用“每行 4 个 rank panel”的网格（宽屏）。
  - 小屏幕下按断点自适应降列（例如 2 列/1 列）。
- 验收标准：
  - 1920 宽度下可稳定 4 列。
  - 移动端仍可读、可滚动。

## 实施建议（分阶段）
1. 阶段 A：Palace 示例与网格同步机制（含版本注记）。
2. 阶段 B：示例网格构建期内嵌与运行时装载统一化。
3. 阶段 C：MPI 输出增强与前端展示协议稳定化。
4. 阶段 D：OPFS 输出落盘 + 弹窗文件浏览器。
5. 阶段 E：控制面板瘦身与 rank 面板 4 列布局优化。

## 当前执行状态（已启动）

### A1. Palace 同步基线矩阵（当前仓库快照）

> 说明：本表用于后续逐项与 Palace 上游示例核对。
> Source Commit 先留空，待首次对齐时回填。

| Example Key | Config | Mesh Asset(s) | Yew Runtime | Palace Source Commit | Sync Status |
|---|---|---|---|---|---|
| spheres | `examples/spheres/spheres.json` | `examples/spheres/mesh/spheres.msh` | Ready | TBD | Pending Verify |
| rings | `examples/rings/rings.json` | `examples/rings/mesh/rings.msh` | Ready | TBD | Pending Verify |
| adapter | `examples/adapter/adapter.json` | `examples/adapter/mesh/adapter.msh` | Ready | TBD | Pending Verify |
| antenna | `examples/antenna/antenna.json` | `examples/antenna/mesh/antenna.msh` | Ready | TBD | Pending Verify |
| coaxial | `examples/coaxial/coaxial.json` | `examples/coaxial/mesh/coaxial.msh` | Ready | TBD | Pending Verify |
| cpw | `examples/cpw/cpw.json` | `examples/cpw/mesh/cpw_coax.msh` | Ready | TBD | Pending Verify |
| cylinder | `examples/cylinder/cylinder.json` | `examples/cylinder/mesh/cylinder_hex.msh` | Ready | TBD | Pending Verify |
| parallel_plate | `examples/parallel_plate/parallel_plate.json` | `examples/parallel_plate/mesh/plate_2d.msh` | Ready | TBD | Pending Verify |
| sbr_sphere | `examples/sbr_sphere/sbr_sphere.json` | `examples/sbr_sphere/mesh/sphere.msh` | Ready | TBD | Pending Verify |
| transmon | `examples/transmon/transmon.json` | `examples/transmon/mesh/transmon.msh2` | Unimplemented | TBD | Pending Verify |

### A2. 立即执行清单
1. 逐例与 Palace 配置做字段差异核对（Problem/Model/Domains/Boundaries/Solver/Postprocessing）。
2. 对每个示例回填 Palace 来源 commit，并在 PR 描述中附差异摘要。
3. 对同名但多网格变体（如 cpw/coaxial/cylinder）明确 yew 默认使用网格，并记录选择理由。

### A4. 多网格变体默认选择（第一版）

> 目标：避免示例运行时因网格选择不一致导致的“看起来可运行但结果不可比”。

| Example Key | 可选网格文件 | 当前 yew 默认 | 选择理由 | 何时切换 |
|---|---|---|---|---|
| `coaxial` | `coaxial.msh`, `coaxial_ascii.msh` | `coaxial.msh` | 与当前 examples 配置和 yew 内嵌路径一致，优先采用主版本网格 | 若 Palace 主线改为 ascii 版或主版本标注变更时切换 |
| `cpw` | `cpw_coax.msh`, `cpw_coax_0.msh`, `cpw_lumped.msh`, `cpw_lumped_0.msh`, `cpw_wave.msh`, `cpw_wave_0.msh` | `cpw_coax.msh` | 与 `examples/cpw/cpw.json` 当前配置一致，作为 CPW 基线 | 增加“端口类型/边界模式”选项后，可按 UI 模式切换对应网格 |
| `cylinder` | `cylinder_hex.msh`, `cylinder_prism.msh`, `cylinder_tet.msh` | `cylinder_hex.msh` | 当前 yew 示例已走 hex 网格，稳定且与示例配置一致 | 后续增加“网格类型”下拉后，允许在 hex/prism/tet 间切换 |

### A5. 规则补充
- 若示例目录下存在多个 mesh 变体：
  - 必须在本表明确“默认值 + 理由 + 切换条件”。
  - yew 运行态与 `examples/*.json` 默认值保持一致。
  - 若运行态默认与配置默认不一致，必须在 PR 中写明偏差原因。

### A3. 同步状态标记规则
- `Pending Verify`：仅完成本仓库盘点，尚未对齐上游。
- `In Sync`：配置与网格已对齐，且记录了 Source Commit。
- `Diverged`：与上游存在可见偏差，需说明原因与计划。
- `Intentional Fork`：有意识地与上游不同，必须附设计理由。

## 备注
- 本文档为改造基线，后续实现可在对应阶段补充“已完成/待完成”与具体 PR 链接。
