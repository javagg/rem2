//! DDM 后处理：解重组与输出

use num_complex::Complex64;
use rem_core::RemResult;

/// 将各子域解合并为全局解向量
pub fn merge_solutions(
    subdomain_solutions: &[Vec<Complex64>],
    global_to_local: &[Vec<usize>],
    n_global: usize,
) -> RemResult<Vec<Complex64>> {
    let mut global = vec![Complex64::ZERO; n_global];
    for (sub_idx, sol) in subdomain_solutions.iter().enumerate() {
        let mapping = &global_to_local[sub_idx];
        for (local_dof, &val) in sol.iter().enumerate() {
            if local_dof < mapping.len() {
                global[mapping[local_dof]] = val;
            }
        }
    }
    Ok(global)
}
