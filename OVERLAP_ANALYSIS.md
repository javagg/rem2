SUMMARY: fem-rs to rem2 overlap analysis

EXECUTIVE SUMMARY
=================

fem-rs provides fem2 significant overlap in 5 major areas:

1. Mesh & I/O: SimplexMesh<D>, GMSH reader → ✅ STRONG MATCH
2. Sparse matrices: CooMatrix<T>, CsrMatrix<T> → ✅ FULL MATCH  
3. Error handling: FemError, FemResult → ✅ STRUCTURAL
4. Assembly: Stiffness/mass forms → 🔶 PARTIAL (different paradigms)
5. Linear solvers: CG, PCG, GMRES, AMG, direct → ✅ FULL COVERAGE

DETAILED FINDINGS
=================

PHASE 1 CANDIDATES (Low Risk, Immediate ROI):
==============================================

1. Sparse Matrices
   Replace: rem_core::{TripletMatrix, CsrMatrix}
   With: fem_linalg::{CooMatrix<f64>, CsrMatrix<f64>}
   Benefit: Generic, tested, identical algorithm
   Risk: Minimal
   Effort: 1 week

2. GMSH Reader
   Replace: rmsh-io crate
   With: fem_io::read_msh_file()
   Benefit: Mature v4.1 reader, cleaner API
   Risk: Low (separate I/O layer)
   Effort: 3 days
   NOTE: Needs adapter for BoundaryTag semantics

3. Error Handling  
   Replace: RemError/RemResult
   With: FemError/FemResult + adapter wrapper
   Benefit: Unified error handling
   Risk: Low (error handling is local)
   Effort: 2 days

PHASE 2 CANDIDATES (Medium Risk, High Value):
==============================================

1. Finite Element Spaces
   Add: fem_space::H1Space, fem_element traits
   Benefit: Enable P2/P3 without rewriting element logic
   Risk: Medium (requires testing)
   Effort: 2 weeks

2. Bilinear Form Assembly
   Replace: assemble_stiffness() element loop
   With: fem_assembly::Assembler + standard integrators
   Benefit: Less boilerplate, automatic element handling, mass/Neumann
   Risk: Medium (major refactor)
   Effort: 2-3 weeks

3. Linear Solvers
   Replace: Custom solve_pcg with SSOR
   With: fem_solver options + fem_amg::AmgPrecond
   Benefit: GMRES, direct solvers, AMG without custom code
   Risk: Low-Medium (test convergence)
   Effort: 1 week

PHASE 3 CANDIDATES (Optional, Lower Priority):
==============================================

1. Postprocessing
   Add: fem_assembly::postprocess gradient recovery, error indicators
   Benefit: Professional analysis without custom code
   Risk: Low (supplementary)
   
2. Adaptive Mesh Refinement
   Add: fem_mesh AMR if future projects need it
   Benefit: Uniform/nonconforming refinement, hanging nodes
   Risk: Low (entirely optional)

KEY BLOCKERS
============

1. BoundaryTag SEMANTICS (MEDIUM SEVERITY)
   Problem: fem-rs uses i32 (geometric); rem2 uses enum (physics)
   Solution: Adapter layer in config parsing to map physics tags → i32
   
2. DIRICHLET BC (LOW)
   Problem: Different interfaces
   Solution: Use fem_space::apply_dirichlet; test on benchmarks
   
3. ASSEMBLY PARADIGM (HIGH)
   Problem: rem2 element-centric; fem-rs space-centric
   Solution: Custom integrator wrapper or full refactor (benefit: cleaner)
   
4. SSOR PRECONDITIONING (LOW)
   Problem: fem-rs no SSOR (only Jacobi, ILU, ILDLt)
   Solution: Test ILU(0) first; keep SSOR as fallback
   
5. MPI PARALLELISM (MEDIUM)
   Problem: Different parallel abstractions
   Solution: Check compatibility; may need adapter

fem-rs PACKAGE NAMES (for Cargo.toml)
=====================================

fem-core       → FemError, FemResult, NodeId, ElemId, scalar types
fem-mesh       → SimplexMesh<D>, ElementType, BoundaryTag, AMR functions
fem-element    → ReferenceElement, Lagrange P1/P2/P3, H(curl), H(div), quadrature
fem-space      → H1Space, L2Space, DofManager, constraint helpers
fem-linalg     → CsrMatrix<T>, CooMatrix<T>, Vector<T>, SparsityPattern
fem-assembly   → Assembler, BilinearIntegrator, standard integrators, postprocess
fem-solver     → solve_cg, solve_pcg_jacobi, solve_gmres, solve_sparse_{lu,cholesky,ldlt}
fem-amg        → AmgPrecond, AmgConfig, CoarsenStrategy, solve_amg_cg
fem-io         → read_msh, read_msh_file, VtkWriter, read_matrix_market
fem-parallel   → (MPI support)

RECOMMENDATION
==============

✅ Start with Phase 1 (1-2 weeks, low risk, immediate ROI)
   - Sparse matrices + error handling + GMSH reader
   - Code simplification without risk
   
🔶 Plan Phase 2 after Phase 1 is stable (3-4 weeks)
   - Larger refactor; requires testing
   - High value once complete
   
⏸️  Defer Phase 3 unless explicitly needed

FILE ANALYSIS SUMMARY
====================

Analyzed from fem-rs vendor:
  - Workspace: 13 crates
  - Key exports from core, mesh, element, space, linalg, assembly, solver, amg, io

Analyzed from rem2:
  - Error types: RemError, RemResult
  - Sparse types: TripletMatrix, CsrMatrix
  - Mesh types: RemMesh, Node, Element
  - Assembly: assemble_stiffness (P1 Tri3/Tet4 only)

RESULT: Document saved to c:/Users/lilu/works/rem2/OVERLAP_ANALYSIS.md
