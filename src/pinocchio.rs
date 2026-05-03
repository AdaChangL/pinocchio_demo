use ark_bn254::{Bn254, Fr, G1Projective, G2Projective};
use ark_ec::{pairing::Pairing, CurveGroup, Group};
use ark_ff::{Field, Zero};
use ark_poly::{
    univariate::{DenseOrSparsePolynomial, DensePolynomial},
    DenseUVPolynomial, Polynomial,
};
use ark_std::UniformRand;
use std::time::{Duration, Instant};

// 辅助函数：多项式清理（修剪尾部的0）
fn sanitize(mut p: DensePolynomial<Fr>) -> DensePolynomial<Fr> {
    while let Some(last) = p.coeffs.last() {
        if last.is_zero() { p.coeffs.pop(); } else { break; }
    }
    p
}

pub fn run_benchmark(
    v_polys: &[DensePolynomial<Fr>],
    w_polys: &[DensePolynomial<Fr>],
    y_polys: &[DensePolynomial<Fr>],
    t_poly: &DensePolynomial<Fr>,
    witness: &[Fr],
    num_io: usize, // 公开输入的数量 (Statement)
) -> (Duration, Duration, Duration, bool) {
    let mut rng = ark_std::test_rng();
    let num_vars = witness.len();

    // =====================================================================
    // 1. Setup Phase (可信设置阶段)
    // =====================================================================
    let setup_start = Instant::now();
    let s = Fr::rand(&mut rng);
    let alpha_v = Fr::rand(&mut rng);
    let alpha_w = Fr::rand(&mut rng);
    let alpha_y = Fr::rand(&mut rng);
    let beta = Fr::rand(&mut rng);
    let gamma = Fr::rand(&mut rng);

    // Pinocchio 的核心偏移参数 r_v, r_w 和 r_y
    let r_v = Fr::rand(&mut rng);
    let r_w = Fr::rand(&mut rng);
    let r_y = r_v * r_w; // r_y 必须等于 r_v * r_w

    let g1 = G1Projective::generator();
    let g2 = G2Projective::generator();

    // Verification Key (VK) 的元素
    let vk_g2 = g2;
    let vk_t_s_g2 = g2 * (t_poly.evaluate(&s) * r_y); // 论文中 t(s) 映射到 g_y 上
    let vk_alpha_v_g2 = g2 * alpha_v;
    let vk_alpha_w_g1 = g1 * alpha_w;
    let vk_alpha_y_g2 = g2 * alpha_y;
    let vk_gamma_g2 = g2 * gamma;
    let vk_beta_gamma_g1 = g1 * (beta * gamma);
    let vk_beta_gamma_g2 = g2 * (beta * gamma);

    // [修复4] 实际应用中，VK 还需要包含公开输入 I/O 和常数项的多项式评估值
    // 这里以注释形式标出，在 Benchmark 模拟中暂用明文替代
    // let vk_v_io: Vec<G1Projective> = (0..num_io).map(|k| g1 * (v_polys[k].evaluate(&s) * r_v)).collect();

    let setup_duration = setup_start.elapsed();

    // =====================================================================
    // 2. Prove Phase (证明阶段)
    // =====================================================================
    let prove_start = Instant::now();

    // 零知识盲化因子 (Zero-Knowledge Blinding Factors)
    let delta_v = Fr::rand(&mut rng);
    let delta_w = Fr::rand(&mut rng);
    let delta_y = Fr::rand(&mut rng);

    let mut v_io_x = DensePolynomial::zero();
    let mut v_mid_x = DensePolynomial::zero();
    let mut w_x = DensePolynomial::zero();
    let mut y_x = DensePolynomial::zero();
    let mut z_x = DensePolynomial::zero(); // 用于 Span Check

    // 按 I/O 和 Mid 分离组装多项式
    for i in 0..num_vars {
        let val = witness[i];
        if val.is_zero() { continue; }

        let p_v = &v_polys[i] * val;
        let p_w = &w_polys[i] * val;
        let p_y = &y_polys[i] * val;

        if i < num_io {
            v_io_x = sanitize(&v_io_x + &p_v);
        } else {
            v_mid_x = sanitize(&v_mid_x + &p_v);
        }
        
        w_x = sanitize(&w_x + &p_w);
        y_x = sanitize(&y_x + &p_y);

        // Z 多项式必须包含所有的变量，乘以对应的 r_v, r_w, r_y 缩放参数！
        let z_term_v = &p_v * r_v;
        let z_term_w = &p_w * r_w;
        let z_term_y = &p_y * r_y;
        let z_term_sum = &(&z_term_v + &z_term_w) + &z_term_y;
        z_x = sanitize(&z_x + &z_term_sum);
    }
    let z_x_beta = &z_x * beta;
    // 注入零知识盲化因子，此处直接用没有盲化的代替了
    let v_mid_blinded = sanitize(&v_mid_x + &(t_poly * delta_v));
    let w_blinded = sanitize(&w_x + &(t_poly * delta_w));
    let y_blinded = sanitize(&y_x + &(t_poly * delta_y));
 // 计算 Z 的盲化标量系数: r_v * delta_v + r_w * delta_w + r_y * delta_y
    let delta_z_scalar = (delta_v * r_v) + (delta_w * r_w) + (delta_y * r_y);
    let z_blinding_poly = t_poly * (beta * delta_z_scalar);
    // 现在这里的 &z_x_beta + &z_blinding_poly 就是标准的 & + & 了
    let z_blinded_poly = sanitize(&z_x_beta + &z_blinding_poly);
    // 计算验证方程的左半部分：(V_io + V_mid_blinded) * W_blinded - Y_blinded
    let v_full_blinded = sanitize(&v_io_x + &v_mid_blinded);
    let p_x = sanitize(&(&v_full_blinded * &w_blinded) - &y_blinded);

    // 计算商多项式 H(x)
    let (h_x, remainder) = DenseOrSparsePolynomial::from(p_x.clone())
        .divide_with_q_and_r(&DenseOrSparsePolynomial::from(t_poly.clone()))
        .expect("Division failed");
    
    assert!(remainder.is_zero(), "Witness is invalid!");

    // =====================================================================
    // 生成 Proof 元素 (模拟同态计算：评估后映射到椭圆曲线上)
    // =====================================================================
    let proof_v_mid = g1 * (v_mid_blinded.evaluate(&s) * r_v);
    let proof_w = g2 * (w_blinded.evaluate(&s) * r_w);
    let proof_y = g1 * (y_blinded.evaluate(&s) * r_y);
    let proof_h = g1 * h_x.evaluate(&s); // H(x) 没有偏移，原样映射到 g1

    // =====================================================================
    // KCA (Alpha) 校验项
    // =====================================================================
    let proof_v_mid_alpha = g1 * (v_mid_blinded.evaluate(&s) * r_v * alpha_v);
    let proof_w_alpha = g1 * (w_blinded.evaluate(&s) * r_w * alpha_w); // 注意 W 的 Alpha 仍放在 G1 以优化配对
    let proof_y_alpha = g1 * (y_blinded.evaluate(&s) * r_y * alpha_y);

    // =====================================================================
    // Span/Z 一致性校验项 (模拟使用 CRS 中的 beta 算子)
    // 遵循 Z 的代数定义: \beta * (r_v*V + r_w*W + r_y*Y)
    // =====================================================================
    let z_blinded_val = beta * (
        (v_full_blinded.evaluate(&s) * r_v) + 
        (w_blinded.evaluate(&s) * r_w) + 
        (y_blinded.evaluate(&s) * r_y)
    );
    let proof_z = g1 * z_blinded_poly.evaluate(&s);

    let prove_duration = prove_start.elapsed();

    // =====================================================================
    // 3. Verify Phase (验证阶段)
    // =====================================================================
    let verify_start = Instant::now();

    // Verifier 基于公开输入(Statement)自己计算 V_io 在 s 处的值
    // [核心修复] Verifier 计算的 V_io 也要映射到 g_v 基底 (乘以 r_v)
    let proof_v_io = g1 * (v_io_x.evaluate(&s) * r_v);
    
    // 合并 V (V_io + V_mid_blinded)
    let proof_v_full = proof_v_io + proof_v_mid;

    // A. 核心 QAP 整除性校验： e(g_v^V, g_w^W) == e(g^H, g_y^T) * e(g_y^Y, g)
    // 代码配对: e(V_full(G1), W(G2)) == e(H(G1), T(G2)) + e(Y(G1), G2)
    let check_qap = Bn254::pairing(proof_v_full.into_affine(), proof_w.into_affine()) == 
                    Bn254::pairing(proof_h.into_affine(), vk_t_s_g2.into_affine()) + 
                    Bn254::pairing(proof_y.into_affine(), vk_g2.into_affine());

    // B. KCA Alpha 校验：确保 Prover 用的是合法多项式张成的空间
    let check_alpha_v = Bn254::pairing(proof_v_mid_alpha.into_affine(), vk_g2.into_affine()) == 
                        Bn254::pairing(proof_v_mid.into_affine(), vk_alpha_v_g2.into_affine());
    
    let check_alpha_w = Bn254::pairing(proof_w_alpha.into_affine(), vk_g2.into_affine()) == 
                        Bn254::pairing(vk_alpha_w_g1.into_affine(), proof_w.into_affine()); 
    
    let check_alpha_y = Bn254::pairing(proof_y_alpha.into_affine(), vk_g2.into_affine()) == 
                        Bn254::pairing(proof_y.into_affine(), vk_alpha_y_g2.into_affine());

    // C. Span / Z 一致性校验：确保 V, W, Y 使用相同的变量系数组合
    // 论文要求: e(g^Z, g^\gamma) = e(g_v^V * g_w^W * g_y^Y, g^{\beta\gamma})
    let check_span = Bn254::pairing(proof_z.into_affine(), vk_gamma_g2.into_affine()) == 
                     Bn254::pairing(proof_v_full.into_affine(), vk_beta_gamma_g2.into_affine()) + 
                     Bn254::pairing(vk_beta_gamma_g1.into_affine(), proof_w.into_affine()) + 
                     Bn254::pairing(proof_y.into_affine(), vk_beta_gamma_g2.into_affine());

    let is_valid = check_qap && check_alpha_v && check_alpha_w && check_alpha_y && check_span;

    let verify_duration = verify_start.elapsed();

    (setup_duration, prove_duration, verify_duration, is_valid)
}