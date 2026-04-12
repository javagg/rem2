# 后端抽象化实现计划

> 目标：统一 rem-mom、rem-febi、rem-ddm 的矩阵-求解器接口
> 优先级：P1（预计 2-3 周）
> 起始日期：2026-04-12

---

## 阶段 1：定义 LinearOperator Trait（第 1 周）

### 1.1 创建 crates/core/src/operator.rs

关键特性：
- 泛型支持 f64 和 Complex64
- 三种矩阵类型：密集、稀疏、矩阵自由
- 适配已有代码最小改动

```rust
// crates/core/src/operator.rs (新文件)

use nalgebra::{DMatrix, DVector};
use num_complex::Complex64;

/// 统一的矩阵-向量操作接口
pub trait LinearOperator<T> {
    /// 矩阵维数
    fn size(&self) -> (usize, usize);
    
    /// 转置矩阵维数
    fn transposed_size(&self) -> (usize, usize) {
        let (m, n) = self.size();
        (n, m)
    }
    
    /// y ← A * x（不改变 y 的零部分）
    fn matvec(&self, x: &DVector<T>, y: &mut DVector<T>) -> Result<(), String>;
    
    /// y ← A^T * x（伴随矩阵用于复数）
    fn matvec_adjoint(&self, x: &DVector<T>, y: &mut DVector<T>) -> Result<(), String> {
        Err("not implemented".to_string())
    }
    
    /// 对角线提取（可选，用于预处理）
    fn diagonal(&self) -> Option<DVector<T>> {
        None
    }
}

/// 求解器 trait
pub trait SparseSolver<T> {
    fn solve(&mut self, op: &dyn LinearOperator<T>, b: &DVector<T>) -> Result<DVector<T>, String>;
}

/// 稠密矩阵的 LinearOperator 适配器
impl LinearOperator<f64> for DMatrix<f64> {
    fn size(&self) -> (usize, usize) {
        (self.nrows(), self.ncols())
    }
    
    fn matvec(&self, x: &DVector<f64>, y: &mut DVector<f64>) -> Result<(), String> {
        if x.len() != self.ncols() || y.len() != self.nrows() {
            return Err(format!("dimension mismatch: matrix {}×{}, x len {}, y len {}",
                self.nrows(), self.ncols(), x.len(), y.len()));
        }
        *y = self * x;
        Ok(())
    }
    
    fn diagonal(&self) -> Option<DVector<f64>> {
        Some(self.diagonal())
    }
}

impl LinearOperator<Complex64> for DMatrix<Complex64> {
    fn size(&self) -> (usize, usize) {
        (self.nrows(), self.ncols())
    }
    
    fn matvec(&self, x: &DVector<Complex64>, y: &mut DVector<Complex64>) -> Result<(), String> {
        if x.len() != self.ncols() || y.len() != self.nrows() {
            return Err(format!("dimension mismatch"));
        }
        *y = self * x;
        Ok(())
    }
    
    fn matvec_adjoint(&self, x: &DVector<Complex64>, y: &mut DVector<Complex64>) -> Result<(), String> {
        if x.len() != self.nrows() || y.len() != self.ncols() {
            return Err(format!("dimension mismatch"));
        }
        *y = self.adjoint() * x;
        Ok(())
    }
    
    fn diagonal(&self) -> Option<DVector<Complex64>> {
        Some(self.diagonal())
    }
}
```

### 1.2 修改 crates/core/src/lib.rs

添加模块导出：
```rust
pub mod operator;
pub use operator::{LinearOperator, SparseSolver};
```

### 1.3 迁移 rem-mom GMRES

**目标**：将 rem-mom 的 GMRES 改成接受 `dyn LinearOperator<Complex64>`

```rust
// crates/mom/src/assemble.rs (修改)

use rem_core::LinearOperator;

/// 改进后的 GMRES 签名
pub fn gmres_solve<'a>(
    op: &dyn LinearOperator<Complex64>,
    rhs: &DVector<Complex64>,
    restart: usize,
    tol: f64,
    max_iter: usize,
) -> RemResult<DVector<Complex64>> {
    let n = op.size().0;
    // ... GMRES 实现，使用 op.matvec() 代替直接矩阵操作
}
```

**变更点**：
- `assemble_efie_pulse()` 返回实现 `LinearOperator` 的 adapter，而非直接 `DMatrix`
- `gmres_solve_aca()` 同样适配

**验证**：
```bash
cargo test -p rem-mom --lib assemble 2>&1 | grep -E "test result|PASSED"
```

---

## 阶段 2：CSR Complex 适配（第 2 周）

### 2.1 扩展 crates/core/src/sparse.rs

为 `CsrMatrix` 添加复数版本：

```rust
// crates/core/src/sparse.rs (添加)

/// Sparse matrix in CSR format with Complex64 entries
#[derive(Debug, Clone)]
pub struct CsrMatrixComplex {
    pub nrows: usize,
    pub ncols: usize,
    pub row_ptr: Vec<usize>,
    pub col_idx: Vec<usize>,
    pub values: Vec<Complex64>,
}

impl LinearOperator<Complex64> for CsrMatrixComplex {
    fn size(&self) -> (usize, usize) {
        (self.nrows, self.ncols)
    }
    
    fn matvec(&self, x: &DVector<Complex64>, y: &mut DVector<Complex64>) -> Result<(), String> {
        // ... CSR SpMV 实现
    }
}
```

### 2.2 rem-febi 集成规划

- 检查 fem-rs 是否已有 `LinearOperator` 接口
- 若无，为 `fem-assembly` 矩阵创建 adapter
- 改进 `FebiSystem` 的求解器调用

---

## 阶段 3：Iterative Solver Framework（第 2-3 周）

### 3.1 实现 GMRES/BiCGSTAB

从 rem-mom 独立出纯粹的求解器供全部 crate 使用：

```rust
// crates/assembly/src/solver.rs (新 crate)

pub struct GMRESSolver {
    restart: usize,
    tol: f64,
    max_iter: usize,
}

impl<T: ComplexField> SparseSolver<T> for GMRESSolver {
    fn solve(&mut self, op: &dyn LinearOperator<T>, b: &DVector<T>) 
        -> Result<DVector<T>, String> {
        // ... 通用 GMRES 实现
    }
}
```

### 3.2 分层求解器选择

根据问题类型自动选择：

```rust
pub enum SolverBackend {
    Direct,        // LU/Cholesky
    GMRESRestart(usize),
    BiCGSTAB,
    CG,            // 仅对称
    AMG,           // 代数多重网格
}

pub fn auto_select_solver(op: &dyn LinearOperator, config: &SolverConfig) -> Box<dyn SparseSolver> {
    match config.solver_type {
        "auto" => {
            if op.size().0 < 5000 {
                Box::new(DirectSolver)
            } else {
                Box::new(GMRESSolver::new(30))
            }
        }
        ...
    }
}
```

---

## 阶段 4：rem-febi / rem-ddm 适配（第 3 周）

### 4.1 rem-febi

改进 `solver::solve_febi` 使用新接口：

```rust
pub fn solve_febi(
    z_bi: &dyn LinearOperator<Complex64>,
    rhs: &DVector<Complex64>,
    config: &FeBiSolverConfig,
) -> RemResult<Vec<Complex64>> {
    let solver = auto_select_solver(z_bi, &config.linear_solver);
    solver.solve(z_bi, &DVector::from_vec(rhs.clone()))
}
```

### 4.2 rem-ddm

在 Schwarz 迭代中每个子域使用通用求解器：

```rust
pub fn schwarz_iteration(
    subdomains: &[SubDomain],
    solver_factory: impl Fn() -> Box<dyn SparseSolver>,
    ...
) -> RemResult<Vec<Complex64>> {
    for iteration in 0..max_iter {
        for subdomain in subdomains {
            let solver = solver_factory();
            let solution = solver.solve(&subdomain.operator(), &subdomain.rhs())?;
            // ...
        }
    }
}
```

---

## 验收标准

### ✅ 代码审查

- [ ] `crates/core/src/operator.rs` 创建，trait 文档完整
- [ ] 为 `DMatrix<Complex64>` 提供 impl
- [ ] rem-mom GMRES 改用 trait（零行为改变）
- [ ] 所有现有单元测试通过

### ✅ 集成测试

```bash
# 运行统一的后端测试
cargo test --test backend_integration -- --nocapture 2>&1 | tail -20
```

### ✅ 性能基准

```bash
cargo bench --features="bench" -- --baseline backend_v0
```

期望：GMRES 性能无退化（<5% overhead）

### ✅ 文档

- [ ] `BACKEND_ABSTRACTION_PLAN.md` (本文件)
- [ ] `crates/core/README_OPERATOR.md` - trait 使用指南
- [ ] examples/backend_abstraction.rs - 示例代码

---

## 里程碑

```
┌─────────────────────────────────────────────────────┐
│ Week 1 (Apr 12-18): LinearOperator trait 定义      │
│ - 定义 operator.rs                                  │
│ - 为 DMatrix<Complex64> 实现                        │
│ - rem-mom GMRES 初步适配                            │
│ - 预期 PR merge                                     │
├─────────────────────────────────────────────────────┤
│ Week 2 (Apr 19-25): Complex CSR / fem-rs 集成      │
│ - CsrMatrixComplex 实现                             │
│ - fem-rs LinearOperator adapter                     │
│ - rem-febi solver 改进                              │
│ - 预期代码覆盖率 > 85%                               │
├─────────────────────────────────────────────────────┤
│ Week 3 (Apr 26-May 2): rem-ddm / 求解框架统一      │
│ - Schwarz 迭代改用通用求解器                        │
│ - DDM 测试通过                                      │
│ - 性能基准建立                                      │
│ - 预期 commit f60509c 后更新所有子模块              │
└─────────────────────────────────────────────────────┘
```

---

## 风险与缓解

| 风险 | 影响 | 缓解 |
|------|------|------|
| trait 太通用导致性能退化 | 5-10% 速度损失 | 使用关联类型 + inline hints |
| fem-rs 接口不兼容 | rem-febi 无法集成 | 提前与 fem-rs 讨论 / 创建 adapter |
| 复数矩阵 trait 造型增加复杂度 | 代码可读性↓ | 类型别名 + 宏减少样板 |
| 现有代码 review 时间长 | 日程延迟 | 优先 rem-mom，后续逐步 |

---

## 外部依赖

- ✅ nalgebra (已有)
- ✅ num-complex (已有)
- ❌ fem-rs LinearOperator 接口（需确认）

## 相关 Issue / PR

- 前置：commit f60509c（fem-rs 最新）
- 后置：Maxwell 验证（P1 优先级 #2）

---

## 参考阅读

- fem-rs backend interface: `vendor/fem-rs/crates/assembly/src/backend.rs`
- rem-mom GMRES 实现: `crates/mom/src/assemble.rs` line 387
- REM 后端抽象化分析：`IMPROVEMENTS_FROM_FEMRS.md`
