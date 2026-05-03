use serde_json::json;
use std::fs::File;
use std::io::Write;

fn main() {
    println!("=== R1CS 电路生成器 ===");

    // 调节电路规模。要求为2^n
    let num_gates = 1024; 
    let num_vars = num_gates + 2; // [常数1, 输入x, v1, v2, ..., v_n]

    println!("正在生成包含 {} 个门, {} 个变量的电路...", num_gates, num_vars);

    let mut a_matrix = vec![vec![0u32; num_vars]; num_gates];
    let mut b_matrix = vec![vec![0u32; num_vars]; num_gates];
    let mut c_matrix = vec![vec![0u32; num_vars]; num_gates];
    
    // 生成一条巨大的线性链：v_{i} * 1 = v_{i+1}
    // 索引映射：常数1(idx 0), x(idx 1), v1(idx 2), v2(idx 3)...
    for i in 0..num_gates {
        let left_input_idx = i + 1;  // 当前门左输入 (x, v1, v2...)
        let right_input_idx = 0;     // 当前门右输入 (永远是常数 1)
        let output_idx = i + 2;      // 当前门输出 (v1, v2, v3...)

        a_matrix[i][left_input_idx] = 1;
        b_matrix[i][right_input_idx] = 1;
        c_matrix[i][output_idx] = 1;
    }

    // 生成对应的 Witness：只要输入 x=3，后面全是 3
    let mut witness = vec![1u32, 3u32];
    for _ in 0..num_gates {
        witness.push(3u32);
    }

    // 构造 JSON 对象
    let circuit_json = json!({
        "num_gates": num_gates,
        "num_vars": num_vars,
        "r1cs": {
            "A": a_matrix,
            "B": b_matrix,
            "C": c_matrix
        },
        "witness": witness
    });

    // 写入文件
    let mut file = File::create("circuit.json").expect("无法创建 circuit.json");
    let json_string = serde_json::to_string(&circuit_json).expect("JSON 序列化失败");
    file.write_all(json_string.as_bytes()).expect("写入文件失败");

    println!("✅ 成功生成巨型电路并保存至 circuit.json");
    println!("文件大小约为: {:.2} MB", json_string.len() as f64 / 1_048_576.0);
}