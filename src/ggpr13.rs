use ark_bn254::{Bn254, Fr, G1Projective, G2Projective};
use ark_ec::{pairing::Pairing, CurveGroup, Group}; 
use ark_ff::{Field, Zero, One};
use ark_poly::{
    univariate::{DensePolynomial, DenseOrSparsePolynomial}, 
    DenseUVPolynomial, Polynomial,
};
use ark_std::UniformRand;
use std::time::{Duration, Instant};

// =====================================================================
// 辅助函数：多项式清理
// =====================================================================
fn sanitize(mut p: DensePolynomial<Fr>) -> DensePolynomial<Fr> {
    while let Some(last) = p.coeffs.last() {
        if last.is_zero() {
            p.coeffs.pop();
        } else {
            break;
        }
    }
    p
}

// =====================================================================
// 辅助函数：拉格朗日插值 (Lagrange Interpolation)
// =====================================================================
fn lagrange_interpolate(x_vals: &[Fr], y_vals: &[Fr]) -> DensePolynomial<Fr> {
    let mut result = DensePolynomial::zero();
    for i in 0..x_vals.len() {
        let mut basis = DensePolynomial::from_coefficients_vec(vec![Fr::one()]);
        let mut denominator = Fr::one();
        for j in 0..x_vals.len() {
            if i != j {
                let term = DensePolynomial::from_coefficients_vec(vec![-x_vals[j], Fr::one()]);
                basis = &basis * &term;
                denominator *= x_vals[i] - x_vals[j];
            }
        }
        basis = &basis * &DensePolynomial::from_coefficients_vec(vec![y_vals[i] * denominator.inverse().unwrap()]);
        result = &result + &basis;
    }
    sanitize(result)
}

// =====================================================================
// 辅助函数：基于 CRT 的多项式强制插值
// =====================================================================
fn crt_interpolate(
    p_a: &DensePolynomial<Fr>,
    t_poly: &DensePolynomial<Fr>,
    points: &[Fr],
    target_values: &[Fr]
) -> DensePolynomial<Fr> {
    let mut q_targets = Vec::with_capacity(points.len());
    for i in 0..points.len() {
        let p_a_val = p_a.evaluate(&points[i]);
        let t_val = t_poly.evaluate(&points[i]);
        let q_target = (target_values[i] - p_a_val) * t_val.inverse().unwrap();
        q_targets.push(q_target);
    }
    let q_poly = lagrange_interpolate(points, &q_targets);
    sanitize(p_a + &(t_poly * &q_poly))
}

// =====================================================================
// GGPR13 强 QAP (Strong QAP) 严格代数转换
// =====================================================================
fn elevate_to_strong_qap(
    v_polys: &[DensePolynomial<Fr>],
    w_polys: &[DensePolynomial<Fr>],
    y_polys: &[DensePolynomial<Fr>],
    t_poly: &DensePolynomial<Fr>,
    num_vars: usize, // m
) -> (Vec<DensePolynomial<Fr>>, Vec<DensePolynomial<Fr>>, Vec<DensePolynomial<Fr>>, DensePolynomial<Fr>) {
    let mut rng = ark_std::test_rng();
    let m = num_vars-1; 
    
    // 选取 2m 个不是 t(x) 的根的随机点
    let mut rs_points = Vec::with_capacity(2 * m);
    while rs_points.len() < 2 * m {
        let rand_pt = Fr::rand(&mut rng);
        if !t_poly.evaluate(&rand_pt).is_zero() && !rs_points.contains(&rand_pt) {
            rs_points.push(rand_pt);
        }
    }
    
    let mut strong_v = Vec::with_capacity(m + 1);
    let mut strong_w = Vec::with_capacity(m + 1);
    let mut strong_y = Vec::with_capacity(m + 1);
    
    // CRT 映射
    for k in 0..=m {
        let mut v_targets = vec![Fr::zero(); 2 * m];
        let mut w_targets = vec![Fr::zero(); 2 * m];
        let mut y_targets = vec![Fr::zero(); 2 * m];
        
        if k == 0 {
            for i in 0..m {
                v_targets[m + i] = Fr::one();
                w_targets[i] = Fr::one();
            }
        } else {
            let idx = k - 1; 
            v_targets[idx] = Fr::one();
            w_targets[m + idx] = Fr::one();
            y_targets[idx] = Fr::one();
            y_targets[m + idx] = Fr::one();
        }
        
        strong_v.push(crt_interpolate(&v_polys[k], t_poly, &rs_points, &v_targets));
        strong_w.push(crt_interpolate(&w_polys[k], t_poly, &rs_points, &w_targets));
        strong_y.push(crt_interpolate(&y_polys[k], t_poly, &rs_points, &y_targets));
    }
    
    // 构建新的目标多项式 t'(x) = t(x) * \prod (x - r_i)(x - s_i)
    let mut t_prime: DensePolynomial<Fr> = (*t_poly).clone(); 
    
    for pt in &rs_points {
        let term = DensePolynomial::from_coefficients_vec(vec![-*pt, Fr::one()]);
        t_prime = sanitize(&t_prime * &term);
    }
    
    (strong_v, strong_w, strong_y, t_prime)
}

pub fn run_benchmark(
    v_polys: &[DensePolynomial<Fr>],
    w_polys: &[DensePolynomial<Fr>],
    y_polys: &[DensePolynomial<Fr>],
    t_poly: &DensePolynomial<Fr>,
    witness: &[Fr],
    _num_gates: usize,
    num_vars: usize,
    num_io: usize,
) -> (Duration, Duration, Duration, bool, Duration) {
    let mut rng = ark_std::test_rng();

    // 0. Preprocessing: 严格的 GGPR13 Strong QAP 转换
    let setup_start = Instant::now();
    let (str_v, str_w, str_y, str_t) = elevate_to_strong_qap(v_polys, w_polys, y_polys, t_poly, num_vars);
    let preprocessiong_duration = setup_start.elapsed();

    // ---------------------------------------------------------------------
    // 1. Setup Phase
    // ---------------------------------------------------------------------
    let setup_start = Instant::now();
    let s = Fr::rand(&mut rng);
    let alpha = Fr::rand(&mut rng);
    let beta_v = Fr::rand(&mut rng);
    let beta_w = Fr::rand(&mut rng);
    let beta_y = Fr::rand(&mut rng);
    let gamma = Fr::rand(&mut rng);

    let g1 = G1Projective::generator();
    let g2 = G2Projective::generator();

    let vk_g2 = g2;
    let vk_alpha_g1 = g1 * alpha; 
    let vk_alpha_g2 = g2 * alpha; 
    let vk_gamma = g2 * gamma;
    let vk_beta_v_gamma = g2 * (beta_v * gamma);
    let vk_beta_w_gamma = g1 * (beta_w * gamma); 
    let vk_beta_y_gamma = g2 * (beta_y * gamma);
    let vk_t_s = g2 * str_t.evaluate(&s);

    let setup_duration = setup_start.elapsed();
    // ---------------------------------------------------------------------
    // 2. Prove Phase
    // ---------------------------------------------------------------------
    let prove_start = Instant::now();

    let mut v_io_x = DensePolynomial::zero();
    let mut v_mid_x = DensePolynomial::zero();
    let mut w_x = DensePolynomial::zero();
    let mut y_x = DensePolynomial::zero();

    for i in 0..num_vars {
        let val = witness[i];
        if val.is_zero() { continue; }
        
        let p_v = &str_v[i] * val;
        let p_w = &str_w[i] * val;
        let p_y = &str_y[i] * val;

        if i < num_io { 
            v_io_x = sanitize(&v_io_x + &p_v); 
        } else { 
            v_mid_x = sanitize(&v_mid_x + &p_v); 
        }
        
        w_x = sanitize(&w_x + &p_w);
        y_x = sanitize(&y_x + &p_y);
    }

    let v_full_x = sanitize(&v_io_x + &v_mid_x);
    let p_x = sanitize(&(&v_full_x * &w_x) - &y_x);
    
    // 使用 clone 获取所有权，适配 DenseOrSparsePolynomial::from
    let (h_x, _) = DenseOrSparsePolynomial::from(p_x.clone())
        .divide_with_q_and_r(&DenseOrSparsePolynomial::from(str_t.clone()))
        .expect("Division failed: Polynomial does not divide evenly. Check witness!");

    let proof_v_mid = g1 * v_mid_x.evaluate(&s);
    let proof_w = g2 * w_x.evaluate(&s); 
    let proof_y = g1 * y_x.evaluate(&s);
    let proof_h = g1 * h_x.evaluate(&s);

    let proof_v_mid_alpha = proof_v_mid * alpha;
    let proof_w_alpha = g1 * (w_x.evaluate(&s) * alpha); 
    let proof_y_alpha = proof_y * alpha;
    let proof_h_alpha = proof_h * alpha;

    let z_val = (beta_v * v_mid_x.evaluate(&s)) + (beta_w * w_x.evaluate(&s)) + (beta_y * y_x.evaluate(&s));
    let proof_z = g1 * z_val; 

    let prove_duration = prove_start.elapsed();

    // ---------------------------------------------------------------------
    // 3. Verify Phase
    // ---------------------------------------------------------------------
    let verify_start = Instant::now();

    let proof_v_io = g1 * v_io_x.evaluate(&s); 
    let proof_v_full = proof_v_io + proof_v_mid;

    // 所有需要配稳的群元素全部使用 into_affine() 转换为仿射坐标
    let check_qap = Bn254::pairing(proof_v_full.into_affine(), proof_w.into_affine()) == 
                    Bn254::pairing(proof_h.into_affine(), vk_t_s.into_affine()) + Bn254::pairing(proof_y.into_affine(), vk_g2.into_affine());

    let check_alpha_v = Bn254::pairing(proof_v_mid_alpha.into_affine(), vk_g2.into_affine()) == Bn254::pairing(proof_v_mid.into_affine(), vk_alpha_g2.into_affine());
    let check_alpha_w = Bn254::pairing(proof_w_alpha.into_affine(), vk_g2.into_affine()) == Bn254::pairing(vk_alpha_g1.into_affine(), proof_w.into_affine()); 
    let check_alpha_y = Bn254::pairing(proof_y_alpha.into_affine(), vk_g2.into_affine()) == Bn254::pairing(proof_y.into_affine(), vk_alpha_g2.into_affine());
    let check_alpha_h = Bn254::pairing(proof_h_alpha.into_affine(), vk_g2.into_affine()) == Bn254::pairing(proof_h.into_affine(), vk_alpha_g2.into_affine());

    let check_span = Bn254::pairing(proof_z.into_affine(), vk_gamma.into_affine()) == 
                     Bn254::pairing(proof_v_mid.into_affine(), vk_beta_v_gamma.into_affine()) + 
                     Bn254::pairing(vk_beta_w_gamma.into_affine(), proof_w.into_affine()) + 
                     Bn254::pairing(proof_y.into_affine(), vk_beta_y_gamma.into_affine());

    let is_valid = check_qap && 
                   check_alpha_v && 
                   check_alpha_w && 
                   check_alpha_y && 
                   check_alpha_h && 
                   check_span;

    let verify_duration = verify_start.elapsed();

    (setup_duration, prove_duration, verify_duration, is_valid, preprocessiong_duration)
}