//! FE-BI 后处理：S 参数提取

use num_complex::Complex64;
use rem_core::RemResult;

/// 从 FEM 求解向量中提取 S 参数。
///
/// 简化实现：将端口节点处的平均电位作为端口电压，
/// 计算 S11 = (Z_in − Z0) / (Z_in + Z0)，Z0 = 50 Ω。
///
/// `port_node_ranges[i]` = 第 i 个端口对应的 DOF 索引范围（起始, 结束）。
/// 当前占位：直接用前 n_ports 个 DOF 的均值作为端口电压。
pub fn extract_sparams(
    solution: &[Complex64],
    n_ports: usize,
    _freq: f64,
) -> RemResult<Vec<Complex64>> {
    let n = solution.len();
    if n == 0 || n_ports == 0 {
        return Ok(vec![Complex64::ZERO; n_ports * n_ports]);
    }

    // 简化：将解向量均分给各端口，取均值作为端口电压
    let chunk = (n / n_ports).max(1);
    let z0 = Complex64::new(50.0, 0.0);
    let mut s = vec![Complex64::ZERO; n_ports * n_ports];

    for i in 0..n_ports {
        let start = (i * chunk).min(n);
        let end   = ((i + 1) * chunk).min(n);
        let v_port: Complex64 = if end > start {
            solution[start..end].iter().sum::<Complex64>()
                / Complex64::new((end - start) as f64, 0.0)
        } else {
            Complex64::ZERO
        };

        // S_ii = (V − Z0) / (V + Z0)  （单端口近似）
        let denom = v_port + z0;
        let s_ii = if denom.norm() > 1e-30 {
            (v_port - z0) / denom
        } else {
            Complex64::ZERO
        };
        s[i * n_ports + i] = s_ii;
    }

    Ok(s)
}
