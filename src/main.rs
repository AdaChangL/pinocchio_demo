mod pinocchio;
mod ggpr13;

use ark_bn254::Fr;
use ark_ff::Zero;
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, Evaluations, GeneralEvaluationDomain, Polynomial,
};
use serde::Deserialize;
use std::fs;

#[derive(Deserialize)]
struct CircuitConfig {
    num_gates: usize,
    num_vars: usize,
    r1cs: R1CSMatrices,
    witness: Vec<u32>,
}

#[derive(Deserialize)]
#[allow(non_snake_case)]
struct R1CSMatrices {
    A: Vec<Vec<u32>>,
    B: Vec<Vec<u32>>,
    C: Vec<Vec<u32>>,
}

fn main() {
    println!("=== ZK-SNARK 对比测试 ===");

    // 读取并解析 JSON
    let config_data = fs::read_to_string("circuit.json").expect("无法读取 circuit.json");
    let config: CircuitConfig = serde_json::from_str(&config_data).expect("JSON 解析失败");
    println!("成功加载电路配置：{} 变量, {} 门", config.num_vars, config.num_gates);

    // 构造列矩阵并进行拉格朗日插值
    let domain = GeneralEvaluationDomain::<Fr>::new(config.num_gates).unwrap();
    let mut a_cols = vec![vec![Fr::zero(); config.num_gates]; config.num_vars];
    let mut b_cols = vec![vec![Fr::zero(); config.num_gates]; config.num_vars];
    let mut c_cols = vec![vec![Fr::zero(); config.num_gates]; config.num_vars];

    for (g, row) in config.r1cs.A.iter().enumerate() { for (v, &val) in row.iter().enumerate() { a_cols[v][g] = Fr::from(val); } }
    for (g, row) in config.r1cs.B.iter().enumerate() { for (v, &val) in row.iter().enumerate() { b_cols[v][g] = Fr::from(val); } }
    for (g, row) in config.r1cs.C.iter().enumerate() { for (v, &val) in row.iter().enumerate() { c_cols[v][g] = Fr::from(val); } }

    let mut reg_v = Vec::new(); let mut reg_w = Vec::new(); let mut reg_y = Vec::new();
    for k in 0..config.num_vars {
        reg_v.push(Evaluations::from_vec_and_domain(a_cols[k].clone(), domain).interpolate());
        reg_w.push(Evaluations::from_vec_and_domain(b_cols[k].clone(), domain).interpolate());
        reg_y.push(Evaluations::from_vec_and_domain(c_cols[k].clone(), domain).interpolate());
    }
    let reg_t: DensePolynomial<Fr> = domain.vanishing_polynomial().into();
    let witness: Vec<Fr> = config.witness.iter().map(|&w| Fr::from(w)).collect();
    
    println!("[Compiler] 基础 R1CS 已成功编译为 Regular QAP");

    // 运行基准测试
    println!("\n正在运行 Pinocchio (Regular QAP) ...");
    let (p_setup, p_prove, p_verify, p_valid) = pinocchio::run_benchmark(&reg_v, &reg_w, &reg_y, &reg_t, &witness);

    println!("正在运行 GGPR13 (Strong QAP) ...");
    let (g_setup, g_prove, g_verify, g_valid) = ggpr13::run_benchmark(&reg_v, &reg_w, &reg_y, &reg_t, &witness, config.num_gates, config.num_vars);
    println!("\n=======================================================");
    println!("Pinocchio vs GGPR13");
    println!("=======================================================");
    println!("指标                | Pinocchio (Regular) | GGPR13 (Strong)");
    println!("--------------------|---------------------|----------------");
    println!("密钥生成 (Setup)    | {:<19?} | {:<15?}", p_setup, g_setup);
    println!("证明生成 (Prove)    | {:<19?} | {:<15?}", p_prove, g_prove);
    println!("证明验证 (Verify)   | {:<19?} | {:<15?}", p_verify, g_verify);
    println!("验证结果            | {:<19} | {:<15}", if p_valid {"通过"} else {"失败"}, if g_valid {"通过"} else {"失败"});
    println!("=======================================================\n");
}