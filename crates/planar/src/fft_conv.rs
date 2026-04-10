// crates/planar/src/fft_conv.rs
//
// 平面 MoM 的 FFT 加速卷积模块
//
// 原理：对于均匀网格上的平面结构，阻抗矩阵 Z[m,n] = f(m-n)（Toeplitz/循环卷积），
// 可用 FFT 将矩阵向量积从 O(N²) 降至 O(N log N)。
//
// 步骤：
//   1. 将格林函数 g(r) 在网格上采样，得到核向量 h
//   2. 对 h 和电流向量 J 做 FFT
//   3. 频域相乘
//   4. IFFT 得到激励向量片段

use rustfft::{FftPlanner, num_complex::Complex};

/// 执行一维循环卷积：y = IFFT(FFT(h) * FFT(x))
/// h 和 x 长度必须相同（调用方负责补零到合适的 2^k 长度）
pub fn circular_conv_1d(h: &[Complex<f64>], x: &[Complex<f64>]) -> Vec<Complex<f64>> {
    assert_eq!(h.len(), x.len(), "h 和 x 长度必须一致");
    let n = h.len();
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(n);
    let ifft = planner.plan_fft_inverse(n);

    let mut h_buf = h.to_vec();
    let mut x_buf = x.to_vec();
    fft.process(&mut h_buf);
    fft.process(&mut x_buf);

    let mut y: Vec<Complex<f64>> = h_buf.iter().zip(x_buf.iter()).map(|(a, b)| a * b).collect();
    ifft.process(&mut y);

    let scale = 1.0 / n as f64;
    y.iter_mut().for_each(|v| *v *= scale);
    y
}

/// 二维 FFT 卷积（行主序，nx 列 × ny 行）
pub fn circular_conv_2d(
    h: &[Complex<f64>],
    x: &[Complex<f64>],
    nx: usize,
    ny: usize,
) -> Vec<Complex<f64>> {
    assert_eq!(h.len(), nx * ny);
    assert_eq!(x.len(), nx * ny);

    let mut planner = FftPlanner::new();
    let fft_x = planner.plan_fft_forward(nx);
    let fft_y = planner.plan_fft_forward(ny);
    let ifft_x = planner.plan_fft_inverse(nx);
    let ifft_y = planner.plan_fft_inverse(ny);

    // --- 前向 FFT ---
    let fft2 = |buf: &mut Vec<Complex<f64>>| {
        // 按行做 FFT（沿 x 方向）
        for row in buf.chunks_mut(nx) {
            fft_x.process(row);
        }
        // 按列做 FFT（沿 y 方向）
        let mut col = vec![Complex::new(0.0, 0.0); ny];
        for ix in 0..nx {
            for iy in 0..ny {
                col[iy] = buf[iy * nx + ix];
            }
            fft_y.process(&mut col);
            for iy in 0..ny {
                buf[iy * nx + ix] = col[iy];
            }
        }
    };

    let mut h_buf = h.to_vec();
    let mut x_buf = x.to_vec();
    fft2(&mut h_buf);
    fft2(&mut x_buf);

    let mut y: Vec<Complex<f64>> = h_buf.iter().zip(x_buf.iter()).map(|(a, b)| a * b).collect();

    // --- 逆 FFT ---
    // 按行 IFFT
    for row in y.chunks_mut(nx) {
        ifft_x.process(row);
    }
    // 按列 IFFT
    let mut col = vec![Complex::new(0.0, 0.0); ny];
    for ix in 0..nx {
        for iy in 0..ny {
            col[iy] = y[iy * nx + ix];
        }
        ifft_y.process(&mut col);
        for iy in 0..ny {
            y[iy * nx + ix] = col[iy];
        }
    }

    let scale = 1.0 / (nx * ny) as f64;
    y.iter_mut().for_each(|v| *v *= scale);
    y
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustfft::num_complex::Complex;

    #[test]
    fn test_conv_1d_delta() {
        // 与 delta 函数卷积应返回原信号
        let n = 8;
        let mut h = vec![Complex::new(0.0, 0.0); n];
        h[0] = Complex::new(1.0, 0.0);
        let x: Vec<Complex<f64>> = (0..n).map(|i| Complex::new(i as f64, 0.0)).collect();
        let y = circular_conv_1d(&h, &x);
        for (i, v) in y.iter().enumerate() {
            assert!((v.re - x[i].re).abs() < 1e-10, "idx={i}: {v} != {}", x[i]);
        }
    }
}
