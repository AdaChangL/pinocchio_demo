use ark_bn254::{Bn254, Fr, G1Projective, G2Projective};
use ark_ec::{pairing::Pairing, Group};
use ark_ff::Zero;
use ark_poly::{
    univariate::{DenseOrSparsePolynomial, DensePolynomial},
    DenseUVPolynomial, Polynomial,
};
use ark_std::UniformRand;
use std::time::{Duration, Instant};

pub fn run_benchmark(
    v_polys: &[DensePolynomial<Fr>],
    w_polys: &[DensePolynomial<Fr>],
    y_polys: &[DensePolynomial<Fr>],
    t_poly: &DensePolynomial<Fr>,
    witness: &[Fr],
) -> (Duration, Duration, Duration, bool) {
    let mut rng = ark_std::test_rng();

    // =====================================================================
    // Setup
    // =====================================================================
    let setup_start = Instant::now();
    let secret_s = Fr::rand(&mut rng);
    
    let g1 = G1Projective::generator();
    let g2 = G2Projective::generator();
    
    // 生成 Verification Key
    // 计算目标多项式 T(x) 在盲点 s 处的值，并放在 G2 曲线上。
    let vk_t = g2 * t_poly.evaluate(&secret_s);
    let setup_duration = setup_start.elapsed();

    // =====================================================================
    // Prove Phase
    // =====================================================================
    let prove_start = Instant::now();
    let num_vars = witness.len();

    let mut v_x = DensePolynomial::<Fr>::zero();
    let mut w_x = DensePolynomial::<Fr>::zero();
    let mut y_x = DensePolynomial::<Fr>::zero();

    // 计算对应多项式
    for i in 0..num_vars {
        v_x = v_x + DensePolynomial::from_coefficients_vec(v_polys[i].coeffs.iter().map(|c| *c * witness[i]).collect());
        w_x = w_x + DensePolynomial::from_coefficients_vec(w_polys[i].coeffs.iter().map(|c| *c * witness[i]).collect());
        y_x = y_x + DensePolynomial::from_coefficients_vec(y_polys[i].coeffs.iter().map(|c| *c * witness[i]).collect());
    }

    // 验证电路约束等式 V*W - Y
    let vw_x = &v_x * &w_x;
    let p_x = &vw_x - &y_x;

    // 计算商多项式 H(x) = p(x) / T(x) 此处结果因为整除
    let p_wrap = DenseOrSparsePolynomial::from(&p_x);
    let t_wrap = DenseOrSparsePolynomial::from(t_poly);
    let (h_x, remainder) = p_wrap.divide_with_q_and_r(&t_wrap).expect("多项式除法失败");

    // 如果余数不为0，说明 Witness 是错的，放弃生成后续复杂的椭圆曲线证明
    if !remainder.is_zero() {
        return (setup_duration, prove_start.elapsed(), Duration::ZERO, false);
    }

    // 同态隐藏多项式
    // Pinocchio 极其精简，Proof 仅包含极少量的群元素。
    let proof_v = g1 * v_x.evaluate(&secret_s);
    let proof_w = g2 * w_x.evaluate(&secret_s);
    let proof_y = g1 * y_x.evaluate(&secret_s);
    let proof_h = g1 * h_x.evaluate(&secret_s);
    let prove_duration = prove_start.elapsed();

    // =====================================================================
    // Verify Phase
    // =====================================================================
    let verify_start = Instant::now();
    // 检查 e(g^V, g^W) == e(g^H, g^T) * e(g^Y, g)
    let pairing_left = Bn254::pairing(proof_v, proof_w);
    let pairing_right = Bn254::pairing(proof_h, vk_t) + Bn254::pairing(proof_y, g2); // arkworks 中使用 + 表示目标群的乘法
    
    let is_valid = pairing_left == pairing_right;
    let verify_duration = verify_start.elapsed();

    (setup_duration, prove_duration, verify_duration, is_valid)
}