use serde_json::json;
use std::fs::File;
use std::io::Write;

fn main() {
    println!("=== 巨型 R1CS 电路生成器 ===");

    // 【可调节电路规模】：建议使用 1024 门来获取顺畅的测试体验和明显的数据对比
    let num_gates = 64; 
    
    // 变量总数：[常数1, 输入x, v1, v2, ..., v_n]
    let num_vars = num_gates + 2; 
    
    // 【新增】：公开输入的数量 (Statement)
    // 索引 0: 常数 1
    // 索引 1: 输入变量 x
    let num_io = 2; 

    println!("正在生成包含 {} 个门, {} 个变量, {} 个公开输入的电路...", num_gates, num_vars, num_io);

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

    // 生成对应的 Witness：只要输入 x=3，后面全推导为 3
    let mut witness = vec![1u32, 3u32];
    for _ in 0..num_gates {
        witness.push(3u32);
    }

    // 构造 JSON 对象，加入新参数 num_io
    let circuit_json = json!({
        "num_gates": num_gates,
        "num_vars": num_vars,
        "num_io": num_io, // 传递给协议后端的 I/O 拆分参数
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

    println!("✅ 成功生成带 num_io 的巨型电路并保存至 circuit.json");
    println!("文件大小约为: {:.2} MB", json_string.len() as f64 / 1_048_576.0);
}