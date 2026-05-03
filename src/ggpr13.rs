use ark_bn254::{Bn254, Fr, G1Projective, G2Projective};
use ark_ec::{pairing::Pairing, Group};
use ark_ff::{Field, Zero}; // 【修复】：引入 Field Trait 以支持 inverse() 方法
use ark_poly::{
    univariate::{DenseOrSparsePolynomial, DensePolynomial},
    DenseUVPolynomial, Polynomial, 
};
use ark_std::UniformRand;
use std::time::{Duration, Instant};

// =====================================================================
// GGPR13 的强 QAP (Strong QAP) 转换
// 通过度数膨胀 (3d + 2N) 强制系数一致性
// =====================================================================
fn elevate_to_strong_qap(
    v_polys: &[DensePolynomial<Fr>],
    w_polys: &[DensePolynomial<Fr>],
    y_polys: &[DensePolynomial<Fr>],
    t_poly: &DensePolynomial<Fr>,
    num_gates: usize,
    num_vars: usize,
) -> (Vec<DensePolynomial<Fr>>, Vec<DensePolynomial<Fr>>, Vec<DensePolynomial<Fr>>, DensePolynomial<Fr>) {
    // 强 QAP 的目标阶数 d' ≈ 3d + 2N
    // 在多项式中引入了额外的约束点，确保 Prover 在 V, W, Y 中使用的是同一套系数。
    let k = 3 * num_gates + 2 * num_vars; 

    // 构造位移多项式 S(x) = x^k
    let mut s_coeffs = vec![Fr::zero(); k + 1];
    s_coeffs[k] = Fr::from(1u32);
    let s_poly = DensePolynomial::from_coefficients_vec(s_coeffs);

    // 构造平方位移多项式 S^2(x) = x^{2k}
    let mut s2_coeffs = vec![Fr::zero(); 2 * k + 1];
    s2_coeffs[2 * k] = Fr::from(1u32);
    let s2_poly = DensePolynomial::from_coefficients_vec(s2_coeffs);

    // 确保 (V*S)(W*S) - Y*S^2 = H*(T*S^2) 依然平衡
    let strong_v = v_polys.iter().map(|p| p * &s_poly).collect();
    let strong_w = w_polys.iter().map(|p| p * &s_poly).collect();
    let strong_y = y_polys.iter().map(|p| p * &s2_poly).collect(); 
    let strong_t = t_poly * &s2_poly;

    (strong_v, strong_w, strong_y, strong_t)
}

pub fn run_benchmark(
    v_polys: &[DensePolynomial<Fr>],
    w_polys: &[DensePolynomial<Fr>],
    y_polys: &[DensePolynomial<Fr>],
    t_poly: &DensePolynomial<Fr>,
    witness: &[Fr],
    num_gates: usize,
    num_vars: usize,
) -> (Duration, Duration, Duration, bool) {
    let mut rng = ark_std::test_rng();

    // ---------------------------------------------------------------------
    // 0. 预处理：将 Regular QAP 升级为 Strong QAP
    // ---------------------------------------------------------------------
    let (str_v, str_w, str_y, str_t) = elevate_to_strong_qap(v_polys, w_polys, y_polys, t_poly, num_gates, num_vars);

    // ---------------------------------------------------------------------
    // 1. Setup Phase (系统初始化)
    // ---------------------------------------------------------------------
    let setup_start = Instant::now();
    let secret_s = Fr::rand(&mut rng); // 盲点 s (有毒废料)
    
    // 生成用于系数一致性检查的挑战因子 beta 和 gamma
    let beta_v = Fr::rand(&mut rng);
    let beta_w = Fr::rand(&mut rng);
    let beta_y = Fr::rand(&mut rng);
    let gamma = Fr::rand(&mut rng); 

    let g1 = G1Projective::generator();
    let g2 = G2Projective::generator();

    // 生成验证密钥 VK 中的核心群元素
    let vk_t = g2 * str_t.evaluate(&secret_s);
    let vk_beta_v = g2 * beta_v;
    let vk_beta_w = g1 * beta_w; 
    let vk_beta_y = g2 * beta_y;
    let vk_gamma = g2 * gamma;

    let setup_duration = setup_start.elapsed();

    // ---------------------------------------------------------------------
    // 2. Prove Phase (生成证明)
    // ---------------------------------------------------------------------
    let prove_start = Instant::now();

    // 初始化最终的 V(x), W(x), Y(x)
    let mut v_x = DensePolynomial::<Fr>::zero();
    let mut w_x = DensePolynomial::<Fr>::zero();
    let mut y_x = DensePolynomial::<Fr>::zero();

    // 计算包含 v0, w0, y0 的完整线性组合。
    // 在 witness 数组中，索引 0 的位置存放的是常数 1。
    for i in 0..num_vars {
        let val = witness[i];
        v_x = v_x + DensePolynomial::from_coefficients_vec(str_v[i].coeffs.iter().map(|c| *c * val).collect());
        w_x = w_x + DensePolynomial::from_coefficients_vec(str_w[i].coeffs.iter().map(|c| *c * val).collect());
        y_x = y_x + DensePolynomial::from_coefficients_vec(str_y[i].coeffs.iter().map(|c| *c * val).collect());
    }

    // 计算 p(x) 并求解商多项式 H(x)
    let vw_x = &v_x * &w_x;
    let p_x = &vw_x - &y_x;
    let p_wrap = DenseOrSparsePolynomial::from(&p_x);
    let t_wrap = DenseOrSparsePolynomial::from(&str_t);
    let (h_x, remainder) = p_wrap.divide_with_q_and_r(&t_wrap).expect("多项式除法失败");

    if !remainder.is_zero() {
        return (setup_duration, prove_start.elapsed(), Duration::ZERO, false);
    }

    // 生成基础 Proof 元素 (映射到椭圆曲线)
    let proof_v = g1 * v_x.evaluate(&secret_s);
    let proof_w = g2 * w_x.evaluate(&secret_s);
    let proof_y = g1 * y_x.evaluate(&secret_s);
    let proof_h = g1 * h_x.evaluate(&secret_s);
    
    // 【GGPR13 核心修正】：生成一致性证明项 Z
    // Z = (beta_v * V + beta_w * W + beta_y * Y) / gamma
    // 这一步利用 gamma 的逆元来打包所有的系数一致性检查。
    let consistency_val = (beta_v * v_x.evaluate(&secret_s)) + (beta_w * w_x.evaluate(&secret_s)) + (beta_y * y_x.evaluate(&secret_s));
    let proof_z = g1 * (consistency_val * gamma.inverse().expect("无法计算 gamma 的逆元")); 

    let prove_duration = prove_start.elapsed();

    // ---------------------------------------------------------------------
    // 3. Verify Phase (验证阶段)
    // ---------------------------------------------------------------------
    let verify_start = Instant::now();

    // 校验 1：标准 QAP 整除性检查
    let pairing_left = Bn254::pairing(proof_v, proof_w);
    let pairing_right = Bn254::pairing(proof_h, vk_t) + Bn254::pairing(proof_y, g2);

    // 校验 2：GGPR13 强一致性检查 (Span Test)
    // 验证公式：e(Z, gamma) == e(V, beta_v) * e(beta_w, W) * e(Y, beta_y)
    let check_consistency = Bn254::pairing(proof_z, vk_gamma) == 
        Bn254::pairing(proof_v, vk_beta_v) + 
        Bn254::pairing(vk_beta_w, proof_w) + 
        Bn254::pairing(proof_y, vk_beta_y);

    let is_valid = (pairing_left == pairing_right) && check_consistency;
    let verify_duration = verify_start.elapsed();

    (setup_duration, prove_duration, verify_duration, is_valid)
}