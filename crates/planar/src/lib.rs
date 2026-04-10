//! `rem-planar` — 平面分层介质矩量法（Planar MoM + FFT）求解器
//!
//! 对标 Sonnet Suite，支持：
//! - 分层媒质格林函数（谱域传递矩阵法）
//! - 均匀网格 2D FFT 卷积加速 O(N log N)
//! - GMRES 迭代求解

pub mod layered_green;
pub mod grid;
pub mod fft_conv;
pub mod impedance;
pub mod solver;

pub use layered_green::{LayeredMedium, Layer, SpectralGreen};
