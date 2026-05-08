# Pinocchio vs GGPR13 — ZK-SNARK 对比实现

一个基于 Rust 的零知识证明（ZK-SNARK）演示项目，实现了 **Pinocchio** 和 **GGPR13** 两种经典协议的完整流程（Setup → Prove → Verify），并在 **BN254** 椭圆曲线上进行性能对比。

---

## 📚 背景

本项目实现了两篇经典论文的核心协议：

| 协议 | 论文 | QAP 类型 | 特点 |
|------|------|----------|------|

| **GGPR13** | [Quadratic Span Programs and Succinct NIZKs without PCPs] | Strong QAP | 通过 CRT 插值构造Strong QAP |
| **Pinocchio** | [Pinocchio: Nearly Practical Verifiable Computation] | Regular QAP | 使用 `r_v, r_w, r_y` 缩放参数实现一致性校验（Span Check） |
---

## 📂 项目结构

```
pinocchio_demo/
├── Cargo.toml              # Rust 项目配置 & 依赖
├── circuit.json            # 预生成的 R1CS 电路配置（64 门, 66 变量, 2 公开输入）
├── src/
│   ├── main.rs             # 程序入口：读取电路、编译 QAP、运行双协议对比
│   ├── pinocchio.rs        # Pinocchio 协议实现
│   ├── ggpr13.rs           # GGPR13 协议实现（含 Strong QAP 转换）
│   └── bin/
│       └── generate_circuit.rs  # 独立电路生成器（可重新生成 circuit.json）
└── README.md
```

---

## 快速开始

### 环境要求

- **Rust** 2024 edition 或以上
- 已安装 Cargo

### 运行

```bash
# 克隆项目
git clone https://github.com/AdaChangL/pinocchio_demo.git
cd pinocchio_demo

# 运行对比测试（使用预生成的 circuit.json）
cargo run --release --bin pinocchio_demo
```

### 重新生成电路

```bash
cargo run --release --bin generate_circuit
```
可修改电路参数