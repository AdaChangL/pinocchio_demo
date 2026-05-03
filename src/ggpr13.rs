use ark_bn254::{Bn254, Fr, G1Projective, G2Projective};
use ark_ec::{pairing::Pairing, Group};
use ark_ff::Zero;
use ark_poly::{
    univariate::{DenseOrSparsePolynomial, DensePolynomial},
    DenseUVPolynomial, Polynomial, 
};
use ark_std::UniformRand;
use std::time::{Duration, Instant};

// =====================================================================
// 内部核心函数：GGPR13 的多项式度数膨胀 (Degree Inflation)
// =====================================================================
fn elevate_to_strong_qap(
    v_polys: &[DensePolynomial<Fr>],
    w_polys: &[DensePolynomial<Fr>],
    y_polys: &[DensePolynomial<Fr>],
    t_poly: &DensePolynomial<Fr>,
    num_gates: usize,
    num_vars: usize,
) -> (Vec<DensePolynomial<Fr>>, Vec<DensePolynomial<Fr>>, Vec<DensePolynomial<Fr>>, DensePolynomial<Fr>) {
    // Strong QAP：目标多项式的最高次幂 从d提升到3d + 2N。
    // 这是为了在多项式中强行塞入额外的验证约束，防止恶意的 Prover 混用系数。
    let k = 3 * num_gates + 2 * num_vars; 

    // 构造位移多项式 S(x) = x^k
    let mut s_coeffs = vec![Fr::zero(); k + 1];
    s_coeffs[k] = Fr::from(1u32);
    let s_poly = DensePolynomial::from_coefficients_vec(s_coeffs);

    // 构造平方位移多项式 S^2(x) = x^{2k}
    let mut s2_coeffs = vec![Fr::zero(); 2 * k + 1];
    s2_coeffs[2 * k] = Fr::from(1u32);
    let s2_poly = DensePolynomial::from_coefficients_vec(s2_coeffs);

    let strong_v = v_polys.iter().map(|p| p * &s_poly).collect();
    let strong_w = w_polys.iter().map(|p| p * &s_poly).collect();
    // Y 和 T 必须乘以 S^2(x)，以确保等式 (V*S)(W*S) - Y*S^2 = H*(T*S^2) 在数学上平衡
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
    // Pre-computation] 转换为 Strong QAP
    // ---------------------------------------------------------------------
    let (str_v, str_w, str_y, str_t) = elevate_to_strong_qap(v_polys, w_polys, y_polys, t_poly, num_gates, num_vars);

    // ---------------------------------------------------------------------
    // Setup
    // ---------------------------------------------------------------------
    let setup_start = Instant::now();
    let secret_s = Fr::rand(&mut rng);
    let g1 = G1Projective::generator();
    let g2 = G2Projective::generator();
    let vk_t = g2 * str_t.evaluate(&secret_s);
    let setup_duration = setup_start.elapsed();

    // ---------------------------------------------------------------------
    // Prove Phase
    // ---------------------------------------------------------------------
    let prove_start = Instant::now();
    let mut v_x = DensePolynomial::<Fr>::zero();
    let mut w_x = DensePolynomial::<Fr>::zero();
    let mut y_x = DensePolynomial::<Fr>::zero();
    for i in 0..num_vars {
        v_x = v_x + DensePolynomial::from_coefficients_vec(str_v[i].coeffs.iter().map(|c| *c * witness[i]).collect());
        w_x = w_x + DensePolynomial::from_coefficients_vec(str_w[i].coeffs.iter().map(|c| *c * witness[i]).collect());
        y_x = y_x + DensePolynomial::from_coefficients_vec(str_y[i].coeffs.iter().map(|c| *c * witness[i]).collect());
    }

    let vw_x = &v_x * &w_x;
    let p_x = &vw_x - &y_x;

    let p_wrap = DenseOrSparsePolynomial::from(&p_x);
    let t_wrap = DenseOrSparsePolynomial::from(&str_t);
    let (h_x, remainder) = p_wrap.divide_with_q_and_r(&t_wrap).expect("多项式除法失败");

    if !remainder.is_zero() {
        return (setup_duration, prove_start.elapsed(), Duration::ZERO, false);
    }

    let proof_v = g1 * v_x.evaluate(&secret_s);
    let proof_w = g2 * w_x.evaluate(&secret_s);
    let proof_y = g1 * y_x.evaluate(&secret_s);
    let proof_h = g1 * h_x.evaluate(&secret_s);
    
    // 【GGPR13 额外开销】：为了保证强一致性，Prover 必须生成额外的椭圆曲线点
    // 这不仅增加了计算时间，也增加了最终发给 Verifier 的 Proof 的体积
    let _proof_extra = g1 * (v_x.evaluate(&secret_s) + w_x.evaluate(&secret_s));

    let prove_duration = prove_start.elapsed();

    // ---------------------------------------------------------------------
    // Verify
    // ---------------------------------------------------------------------
    let verify_start = Instant::now();
    let pairing_left = Bn254::pairing(proof_v, proof_w);
    let pairing_right = Bn254::pairing(proof_h, vk_t) + Bn254::pairing(proof_y, g2);
    let is_valid = pairing_left == pairing_right;
    let verify_duration = verify_start.elapsed();

    (setup_duration, prove_duration, verify_duration, is_valid)
}