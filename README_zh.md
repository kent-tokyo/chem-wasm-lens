# chem-wasm-lens

[English](README.md) | [日本語](README_ja.md)

用 Pure Rust 编写的超轻量分子分析内核，编译为 WebAssembly。设计用于在浏览器 Web Worker 中运行，在不阻塞 UI 线程的情况下，对大型分子结构（10k+ 原子）执行高性能的拓扑分析、距离查询和几何计算。

**[在线演示](https://kent-tokyo.github.io/chem-wasm-lens/examples/)** — [SMILES → SVG](https://kent-tokyo.github.io/chem-wasm-lens/examples/smiles_svg.html) · [3D 查看器](https://kent-tokyo.github.io/chem-wasm-lens/examples/viewer.html) · [真实蛋白质（RCSB）](https://kent-tokyo.github.io/chem-wasm-lens/examples/fetch_demo.html) · [结构式编辑器](https://kent-tokyo.github.io/chem-wasm-lens/examples/editor.html)

## 为什么需要 chem-wasm-lens？

### 问题：浏览器中的分子计算为何困难

基于 Web 的分子工具面临共同的瓶颈。

**1. JavaScript 不适合大规模计算**
现实中的分子——蛋白质、核酸——通常超过 10,000 个原子。原子对距离计算为 O(N²)。在 JS 主线程上运行会导致数秒的 UI 冻结。

**2. 服务端工具引入延迟且无法离线使用**
Python + RDKit 等服务端方案需要搭建和维护后端，存在网络延迟，无法在离线环境或高频交互场景中使用。

**3. 现有 Wasm 分子库体积过大**
RDKit.js（C++ RDKit 的 Wasm 移植版）功能丰富，但 **Bundle 大小超过 ~10MB**。对于只需要距离查询和近邻搜索的应用来说，这远超所需。

### 解决方案

| | chem-wasm-lens | RDKit.js | Python + RDKit | Pure JS |
|---|:---:|:---:|:---:|:---:|
| 浏览器中运行 | 是 | 是 | 否 | 是 |
| 离线可用 | 是 | 是 | 否 | 是 |
| 支持 10k+ 原子 | 是 | 是 | 是 | 否 |
| Bundle 大小 | 小 | 大(~10MB+) | — | 小 |
| 不阻塞 UI | 是 (Web Worker) | 部分支持 | — | 否 |
| 无需安装 | 是 | 是 | 否 | 是 |
| C/C++ 依赖 | 无 | 有 | 有 | 无 |

**chem-wasm-lens 的方法：**

- **Rust → Wasm 实现接近原生的速度** — 无垃圾回收器，在浏览器中实现可预测的低延迟
- **Web Worker 优先设计** — 重计算在独立线程运行，UI 始终保持响应
- **零 C/C++ 依赖** — Pure Rust，`wasm-pack build` 一步完成构建，无需复杂的交叉编译配置
- **专注于分析内核** — 不含 3D 渲染或 UI，可与现有可视化工具（3Dmol.js、NGL Viewer 等）组合使用

### 适用场景

- 在浏览器中直接加载 PDB/XYZ 文件，实时计算原子间距离或邻近残基
- 构建无需后端的离线分子查看器或教育工具
- 介于「RDKit.js 太重」与「Pure JS 太慢」之间的需求

---

## 特性

- **Wasm 优先** — 零原生 C/C++ 依赖；通过 `wasm-pack` 干净构建
- **高性能** — 使用扁平 `Vec<f32>` 坐标布局，实现缓存友好的距离计算
- **零/最小拷贝** — 从 JS 传递原始文件内容；解析和状态保持完全在 Wasm 内完成
- **安全解析** — 解析逻辑中无 `unwrap()`；全面使用 `Result` 类型进行显式错误处理
- **键检测** — Cordero 2008 共价半径表（18 种元素）；通过 `compute_bonds()` 按需计算
- **体素网格空间索引** — 均匀网格将近邻查询加速至平均 O(1)；未构建索引时自动回退到 O(N) 线性扫描
- **serde JSON 输出** — `get_atom_info()` 和 `get_neighbors_info()` 通过 `serde-wasm-bindgen` 返回结构化 JS 对象
- **结构式编辑器内核** — 提供原子/键的增删改、环模板、缩环、矩形选择、价键检查等 API，可在浏览器中构建完整的结构式编辑器

## 状态

| 阶段 | 内容 | 状态 |
|------|------|------|
| 1 | XYZ 解析器、`MolecularSystem` 结构体、Wasm 暴露 | 完成 |
| 1 | PDB 解析器（`ATOM`/`HETATM`）、键检测 | 完成 |
| 2 | 距离查询、半径近邻搜索 | 完成 |
| 2 | 20k+ 原子的体素网格空间索引 | 完成 |
| 3 | CI（GitHub Actions: `cargo test`、`clippy`、`wasm-pack build`） | 完成 |
| 3 | JS/TS 使用示例、浏览器测试、npm 发布、基准测试 | 完成 |
| 4 | SMILES、SVG 2D 渲染、SMARTS、指纹、文件格式 I/O | 完成 |
| 5 | 结构式编辑器内核（P42–P49：23 个方法） | 完成 |
| 6 | RDKit.js 功能对齐（P50–P53：描述符·3D 构象·反应·SMARTS 强化） | 完成 |

## 快速开始

### 前置条件

```sh
# Rust 工具链
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# wasm-pack
cargo install wasm-pack
```

### 构建

```sh
# 原生构建（用于 cargo test）
cargo build

# WebAssembly — 打包器目标（Webpack, Vite 等）
wasm-pack build --target bundler

# WebAssembly — 原生 Web 目标（无需打包器）
wasm-pack build --target web
```

### 测试

```sh
# 原生单元测试
cargo test

# 代码检查
cargo clippy

# 浏览器测试（需要 Chrome）
wasm-pack test --headless --chrome
```

## 使用示例（JavaScript / TypeScript）

### 基础 — 解析并查看原子

```js
import init, { MolecularSystem } from './pkg/chem_wasm_lens.js';

await init();

const xyzData = `3
water
O   0.000  0.000  0.119
H   0.000  0.757 -0.477
H   0.000 -0.757 -0.477
`;

const mol = MolecularSystem.from_xyz_string(xyzData);

console.log(mol.atom_count());   // 3
console.log(mol.get_symbol(0));  // "O"
console.log(mol.get_x(1));      // 0.0
```

### 空间查询 — "lens" 核心用例

```js
import init, { MolecularSystem } from './pkg/chem_wasm_lens.js';

await init();

// 将完整 PDB 文件内容作为字符串传入 — 解析完全在 Wasm 内完成
const mol = MolecularSystem.from_pdb_string(pdbFileContent);

// 构建体素网格空间索引（推荐网格大小：3–5 Å）
// 此步骤可选，但可将 get_atoms_within_radius 加速至平均 O(1)
mol.build_spatial_index(5.0);

// 查询原子索引 42 半径 5 Å 内的所有原子
// 返回 JS Array（普通对象数组）: [{index, symbol, x, y, z, atom_name, ...}, ...]
const centerAtom = 42;
const neighbors = mol.get_neighbors_info(centerAtom, 5.0);

neighbors.forEach(atom => {
  console.log(`${atom.chain_id}:${atom.residue_name}${atom.residue_id} ${atom.atom_name} — dist ≈ 5 Å`);
});

// 也可以只获取半径内的残基标签
const residues = mol.get_residues_within_radius(centerAtom, 5.0);
console.log(residues); // 例如: ["A:ALA:10", "A:GLY:11", ...]
```

## 编辑器内核 API

用于在浏览器中构建结构式编辑器 UI 的编辑原语。

| 方法 | 返回值 | 说明 |
|------|--------|------|
| `add_atom(symbol, x, y)` | `number` | 添加原子，返回新索引 |
| `remove_atom(idx)` | `void` | 删除原子（自动重映射键索引） |
| `set_atom_symbol(idx, symbol)` | `void` | 修改元素符号 |
| `set_atom_position(idx, x, y)` | `void` | 修改 2D 坐标 |
| `set_atom_charge(idx, charge)` | `void` | 设置形式电荷 |
| `add_bond(a, b, order)` | `void` | 添加键（若已存在则更新键级） |
| `remove_bond(a, b)` | `void` | 删除键 |
| `set_bond_order(a, b, order)` | `void` | 修改键级 |
| `closest_atom(x, y, tol)` | `number\|undefined` | 命中测试：最近邻原子 |
| `bond_at(x, y, tol)` | `Uint32Array` | 命中测试：最近邻键 `[a, b]` |
| `normalize_bond_length(target)` | `void` | 将平均键长缩放至 target |
| `translate_atoms(dx, dy)` | `void` | 平移所有原子 |
| `implicit_h_count(idx)` | `number` | 隐式氢数（标准价 − 键级之和） |
| `add_ring_template(n, cx, cy, bond_len)` | `Uint32Array` | 放置正 n 元环，返回新原子索引 |
| `attach_ring_to_bond(a, b, n)` | `Uint32Array` | 在键 a–b 上缩合 n 元环 |
| `get_bounds()` | `Float32Array` | 包围盒 `[min_x, min_y, max_x, max_y]` |
| `rotate_atoms(angle, cx, cy)` | `void` | 绕点 (cx,cy) 旋转所有原子 |
| `flip_horizontal(cx)` | `void` | 以 x=cx 为轴水平镜像 |
| `flip_vertical(cy)` | `void` | 以 y=cy 为轴垂直镜像 |
| `select_atoms_in_rect(x1, y1, x2, y2)` | `Uint32Array` | 矩形框选（自动规范化） |
| `move_atoms(indices, dx, dy)` | `void` | 仅平移指定原子 |
| `check_valence()` | `Uint32Array` | 返回键级超出标准价的原子索引 |
| `copy_atoms(indices)` | `MolecularSystem` | 将原子子集提取为新实例 |

### RDKit.js 功能对齐 — 描述符与分析

| 方法 | 返回值 | 说明 |
|------|--------|------|
| `largest_fragment()` | `MolecularSystem` | 最大连通分量（脱盐处理） |
| `murcko_scaffold()` | `MolecularSystem` | 返回 Murcko 骨架的新实例 |
| `num_heavy_atoms()` | `number` | 不含 H 的原子数 |
| `fraction_csp3()` | `number` | sp3 碳占比（Fsp3） |
| `molar_refractivity()` | `number` | Wildman-Crippen MR 估算值 |
| `embed_molecule(seed)` | `boolean` | 使用距离几何算法生成 3D 坐标 |
| `Reaction.run_reaction(reactant)` | `MolecularSystem[]` | 对每处子结构匹配应用反应模板，返回各匹配位点的产物 |
| `ring_sizes_for_atom(idx)` | `number[]` | 包含原子 idx 的 SSSR 环大小列表 |
| `ring_info()` | `{num_rings, ring_sizes}` | 整个分子的环数和环大小 |
| `aliphatic_ring_count()` | `number` | 非芳香环数量 |
| `fingerprint_atom_pair()` | `Uint8Array` | 2048-bit Atom Pair 指纹（256 字节） |
| `generate_aligned_coords(template)` | `void` | 通过子结构匹配 + 2D Kabsch 旋转将坐标对齐至模板 |

---

## 架构

```
src/lib.rs
├── ParseError          — 错误枚举（实现 Display，不依赖 JsValue）
├── SpatialGrid         — 私有；基于 HashMap 的均匀网格，实现平均 O(1) 近邻查找
├── AtomInfo            — serde::Serialize；get_atom_info / get_neighbors_info 的输出形状
├── MolecularSystem     — 核心结构体（#[wasm_bindgen]）
│   ├── symbols / x / y / z: Vec<f32>  — 分离扁平向量，适合 SIMD
│   ├── atom_names / residue_names / residue_ids / chain_ids / hetatm_flags  — PDB 元数据
│   ├── bonds: Vec<Vec<usize>>          — 通过 compute_bonds() 按需计算
│   └── spatial_grid: Option<SpatialGrid>  — 通过 build_spatial_index() 按需构建
├── parse_xyz()         — 纯 Rust XYZ 解析器，可通过 cargo test 测试
├── parse_pdb()         — 固定宽度列 PDB 解析器（ATOM/HETATM 记录）
├── covalent_radius()   — Cordero 2008 表（18 种元素）
└── impl MolecularSystem (#[wasm_bindgen])
    ├── 解析器:     from_xyz_string(), from_pdb_string()
    ├── 访问器:     atom_count(), get_symbol/x/y/z(), get_atom_name(),
    │               get_residue_name/id(), get_chain_id(), is_hetatm()
    ├── 批量导出:   get_positions_flat() → Float32Array, get_symbols_json()
    ├── 键检测:     compute_bonds(), get_bonds(), bond_count(), has_bonds_computed()
    ├── 空间查询:   distance(), get_atoms_within_radius(), get_residues_within_radius()
    ├── 空间索引:   build_spatial_index(), has_spatial_index()
    ├── JSON 输出:  get_atom_info() → JsValue, get_neighbors_info() → JsValue
    └── 编辑器:     add/remove_atom, add/remove_bond, set_atom_*, set_bond_order,
                    implicit_h_count, add_ring_template, attach_ring_to_bond,
                    get_bounds, rotate_atoms, flip_*, select_atoms_in_rect,
                    move_atoms, check_valence, copy_atoms
```

**主要设计决策：**

- **解析器与 Wasm 边界分离** — 纯 Rust 函数（`parse_xyz`、`parse_pdb`）不含 `JsValue`，可通过 `cargo test` 在无浏览器环境下完整测试
- **坐标向量分离** — `x`、`y`、`z` 作为独立的 `Vec<f32>`，最大化向量化距离计算的缓存局部性
- **传字符串而非对象** — 整个文件内容作为单个 `&str` 跨越 JS/Wasm 边界，解析完全在 Rust 侧进行
- **按需计算** — 键和空间网格均为懒加载；调用方控制计算时机和开销
- **体素网格自动回退** — `get_atoms_within_radius` 在索引已构建时使用网格（平均 O(1)），未调用 `build_spatial_index` 时自动回退到 O(N) 线性扫描

## 许可证

MIT
