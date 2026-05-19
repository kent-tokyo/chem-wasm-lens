# chem-wasm-lens

[English](README.md) | [中文](README_zh.md)

Pure Rust で書かれた超軽量の分子解析カーネルで、WebAssembly にコンパイルされます。ブラウザの Web Worker 内で動作し、UI スレッドをブロックすることなく、大規模な分子構造（10k+ 原子）に対するトポロジー解析・距離クエリ・幾何計算を高速に処理します。

**[ライブデモ](https://kent-tokyo.github.io/chem-wasm-lens/examples/)** — [SMILES → SVG](https://kent-tokyo.github.io/chem-wasm-lens/examples/smiles_svg.html) · [3D ビューア](https://kent-tokyo.github.io/chem-wasm-lens/examples/viewer.html) · [実タンパク質（RCSB）](https://kent-tokyo.github.io/chem-wasm-lens/examples/fetch_demo.html) · [構造式エディタ](https://kent-tokyo.github.io/chem-wasm-lens/examples/editor.html) · [フィンガープリント](https://kent-tokyo.github.io/chem-wasm-lens/examples/fingerprint.html) · [2D アライメント](https://kent-tokyo.github.io/chem-wasm-lens/examples/alignment.html)

## なぜ chem-wasm-lens が必要か

### 課題：ブラウザでの分子計算はなぜ難しいのか

分子構造を扱うウェブアプリには、共通の壁があります。

**1. JavaScript は大規模計算に向かない**
タンパク質や核酸など実用的な分子は 10,000 原子を超えることが多く、すべての原子ペア間の距離計算は O(N²) になります。これを JavaScript のメインスレッドで実行すると、数秒〜数十秒の UI フリーズが発生します。

**2. サーバーサイド処理はオフライン不可・レイテンシが生じる**
Python + RDKit などのサーバー依存ソリューションは、バックエンドの構築・維持が必要で、ネットワーク遅延も避けられません。オフライン環境や高速インタラクティブ操作には向きません。

**3. 既存の Wasm 分子ライブラリはバンドルが重い**
RDKit.js は C++ の RDKit を Wasm ポートしたもので、機能は豊富ですが **バンドルサイズが ~10MB 超**。軽量なウェブアプリや「距離計算・近傍探索だけしたい」ユースケースには過剰です。

### 解決策：chem-wasm-lens のアプローチ

| | chem-wasm-lens | RDKit.js | Python + RDKit | Pure JS |
|---|:---:|:---:|:---:|:---:|
| ブラウザで動作 | Yes | Yes | No | Yes |
| オフライン対応 | Yes | Yes | No | Yes |
| 10k+ 原子対応 | Yes | Yes | Yes | No |
| バンドルサイズ | 小 | 大(~10MB+) | — | 小 |
| UI をブロックしない | Yes (Web Worker) | 条件付き | — | No |
| インストール不要 | Yes | Yes | No | Yes |
| C/C++ 依存 | なし | あり | あり | なし |

**chem-wasm-lens が選ぶアプローチ：**

- **Rust → Wasm でネイティブに近い速度** — ガベージコレクションのない Rust はブラウザ内でも予測可能な低レイテンシを実現
- **Web Worker 前提の設計** — 重い計算はメインスレッドとは別スレッドで実行。スクロールやアニメーションが途切れない
- **C/C++ 依存ゼロ** — Pure Rust なので `wasm-pack build` 一発でビルド完結。複雑なクロスコンパイル設定が不要
- **「解析カーネル」に特化** — 3D レンダリングや UI は持たない。既存のビジュアライザ（3Dmol.js, NGL Viewer など）と組み合わせて使う

### こんなケースに向いています

- PDB / XYZ ファイルをブラウザで直接読み込み、原子間距離や近傍残基をリアルタイムに計算したい
- バックエンドなしで動く分子ビューアや教育ツールを作りたい
- RDKit.js は重すぎる、でも JavaScript だけでは遅すぎる、というニーズ

---

## 特徴

- **Wasm ファースト** — ネイティブの C/C++ 依存ゼロ。`wasm-pack` でクリーンにビルド
- **高パフォーマンス** — キャッシュ効率の高い距離計算ループのために、座標をフラットな `Vec<f32>` で保持
- **ゼロ/最小コピー** — JS からは生のファイル内容を渡すだけ。パースと状態保持はすべて Wasm 側の Rust で完結
- **安全なパース** — パースロジックに `unwrap()` なし。全面的に `Result` 型で明示的なエラーハンドリング
- **ボンド検出** — Cordero 2008 の共有結合半径テーブル（18 元素）。`compute_bonds()` でオンデマンドに計算
- **ボクセルグリッド空間インデックス** — 均一グリッドにより近傍クエリを平均 O(1) に高速化。インデックス未構築時は O(N) 線形探索にフォールバック
- **serde JSON 出力** — `get_atom_info()` と `get_neighbors_info()` が `serde-wasm-bindgen` を介して構造化 JS オブジェクトを返す
- **構造式エディタカーネル** — 原子・結合の追加・削除・変更、環テンプレート、縮環、矩形選択、バレンスチェックなど、ブラウザ上で動く構造式エディタを構築するための API セットを提供

## ステータス

| フェーズ | 内容 | 状態 |
|---------|------|------|
| 1 | XYZ パーサー・`MolecularSystem` 構造体・Wasm 公開 | 完了 |
| 1 | PDB パーサー（`ATOM`/`HETATM`）・ボンド検出 | 完了 |
| 2 | 距離クエリ・半径ベースの近傍探索 | 完了 |
| 2 | 20k+ 原子向けのボクセルグリッド空間インデックス | 完了 |
| 3 | CI（GitHub Actions: `cargo test`・`clippy`・`wasm-pack build`） | 完了 |
| 3 | JS/TS 使用例・ブラウザテスト・npm 公開・ベンチマーク | 完了 |
| 4 | SMILES・SVG 2D 描画・SMARTS・フィンガープリント・ファイル形式 I/O | 完了 |
| 5 | 構造式エディタカーネル（P42–P49: 23 メソッド） | 完了 |
| 6 | RDKit.js パリティ（P50–P53: 記述子・3D コンフォーマー・Reaction・SMARTS 強化） | 完了 |

## クイックスタート

### 前提条件

```sh
# Rust ツールチェーン
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# wasm-pack
cargo install wasm-pack
```

### ビルド

```sh
# ネイティブビルド（cargo test 用）
cargo build

# WebAssembly — バンドラーターゲット（Webpack, Vite など）
wasm-pack build --target bundler

# WebAssembly — バニラウェブターゲット（バンドラー不要）
wasm-pack build --target web
```

### テスト

```sh
# ネイティブユニットテスト
cargo test

# リント
cargo clippy

# ブラウザテスト（Chrome が必要）
wasm-pack test --headless --chrome
```

## 使用例（JavaScript / TypeScript）

### 基本 — 原子のパースと参照

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

### 空間クエリ — "lens" の中核ユースケース

```js
import init, { MolecularSystem } from './pkg/chem_wasm_lens.js';

await init();

// PDB ファイルの内容を文字列で渡す — パースはすべて Wasm 内で完結
const mol = MolecularSystem.from_pdb_string(pdbFileContent);

// ボクセルグリッド空間インデックスを構築（推奨セルサイズ: 3〜5 Å）
// この手順はオプションだが、get_atoms_within_radius の速度を平均 O(1) に引き上げる
mol.build_spatial_index(5.0);

// 原子インデックス 42 から半径 5 Å 以内の全原子を取得
// JS の Array（プレーンオブジェクトの配列）を返す: [{index, symbol, x, y, z, atom_name, ...}, ...]
const centerAtom = 42;
const neighbors = mol.get_neighbors_info(centerAtom, 5.0);

neighbors.forEach(atom => {
  console.log(`${atom.chain_id}:${atom.residue_name}${atom.residue_id} ${atom.atom_name} — dist ≈ 5 Å`);
});

// 半径内の残基ラベルだけ取得することも可能
const residues = mol.get_residues_within_radius(centerAtom, 5.0);
console.log(residues); // 例: ["A:ALA:10", "A:GLY:11", ...]
```

## エディタカーネル API

構造式エディタ UI をブラウザ内で構築するための編集プリミティブです。

| メソッド | 返り値 | 説明 |
|---------|--------|------|
| `add_atom(symbol, x, y)` | `number` | 原子追加（インデックス返却） |
| `remove_atom(idx)` | `void` | 原子削除（結合インデックス自動再マッピング） |
| `set_atom_symbol(idx, symbol)` | `void` | 元素変更 |
| `set_atom_position(idx, x, y)` | `void` | 座標変更 |
| `set_atom_charge(idx, charge)` | `void` | 形式電荷変更 |
| `add_bond(a, b, order)` | `void` | 結合追加（既存なら次数更新） |
| `remove_bond(a, b)` | `void` | 結合削除 |
| `set_bond_order(a, b, order)` | `void` | 結合次数変更 |
| `closest_atom(x, y, tol)` | `number\|undefined` | ヒットテスト：最近傍原子 |
| `bond_at(x, y, tol)` | `Uint32Array` | ヒットテスト：最近傍結合 `[a, b]` |
| `normalize_bond_length(target)` | `void` | 平均結合長を target に正規化 |
| `translate_atoms(dx, dy)` | `void` | 全体平行移動 |
| `implicit_h_count(idx)` | `number` | 暗黙的 H 数（標準原子価 − 結合次数合計） |
| `add_ring_template(n, cx, cy, bond_len)` | `Uint32Array` | n 員環を配置、新原子インデックスを返す |
| `attach_ring_to_bond(a, b, n)` | `Uint32Array` | 既存結合に n 員環を縮環 |
| `get_bounds()` | `Float32Array` | バウンディングボックス `[min_x, min_y, max_x, max_y]` |
| `rotate_atoms(angle, cx, cy)` | `void` | 点 (cx,cy) を中心に全原子を回転 |
| `flip_horizontal(cx)` | `void` | 垂直軸 x=cx に鏡映 |
| `flip_vertical(cy)` | `void` | 水平軸 y=cy に鏡映 |
| `select_atoms_in_rect(x1, y1, x2, y2)` | `Uint32Array` | 矩形内原子インデックス（自動正規化） |
| `move_atoms(indices, dx, dy)` | `void` | 指定原子のみ平行移動 |
| `check_valence()` | `Uint32Array` | バレンス違反原子のインデックス |
| `copy_atoms(indices)` | `MolecularSystem` | 原子サブセットを新インスタンスとして抽出 |

### RDKit.js パリティ — 記述子 & 解析

| メソッド | 返り値 | 説明 |
|---------|--------|------|
| `largest_fragment()` | `MolecularSystem` | 最大連結成分（脱塩処理） |
| `murcko_scaffold()` | `MolecularSystem` | Murcko スキャフォールドを新インスタンスで返す |
| `num_heavy_atoms()` | `number` | H を除いた原子数 |
| `fraction_csp3()` | `number` | sp3 炭素割合（Fsp3） |
| `molar_refractivity()` | `number` | Wildman-Crippen MR 推定値 |
| `embed_molecule(seed)` | `boolean` | 距離幾何学による 3D 座標生成 |
| `Reaction.run_reaction(reactant)` | `MolecularSystem[]` | 反応テンプレートを適用、マッチごとに生成物を返す |
| `ring_sizes_for_atom(idx)` | `number[]` | 原子 idx を含む SSSR 環サイズ一覧 |
| `ring_info()` | `{num_rings, ring_sizes}` | 分子全体の環数と環サイズ |
| `aliphatic_ring_count()` | `number` | 非芳香族環の数 |
| `fingerprint_atom_pair()` | `Uint8Array` | 2048-bit Atom Pair フィンガープリント（256 バイト） |
| `generate_aligned_coords(template)` | `void` | サブストラクチャーマッチ + 2D Kabsch 回転でテンプレート座標系に揃える |

---

## アーキテクチャ

```
src/lib.rs
├── ParseError          — エラー列挙型（Display 実装。JsValue 非依存）
├── SpatialGrid         — 非公開。HashMap ベースの均一グリッド。平均 O(1) の近傍探索を実現
├── AtomInfo            — serde::Serialize。get_atom_info / get_neighbors_info の出力形状
├── MolecularSystem     — コア構造体（#[wasm_bindgen]）
│   ├── symbols / x / y / z: Vec<f32>  — 分離フラットベクタ。SIMD フレンドリー
│   ├── atom_names / residue_names / residue_ids / chain_ids / hetatm_flags  — PDB メタデータ
│   ├── bonds: Vec<Vec<usize>>          — compute_bonds() でオンデマンドに計算
│   └── spatial_grid: Option<SpatialGrid>  — build_spatial_index() でオンデマンドに構築
├── parse_xyz()         — 純 Rust XYZ パーサー。cargo test でテスト可能
├── parse_pdb()         — 固定幅カラム形式 PDB パーサー（ATOM/HETATM レコード）
├── covalent_radius()   — Cordero 2008 テーブル（18 元素）
└── impl MolecularSystem (#[wasm_bindgen])
    ├── パーサー:         from_xyz_string(), from_pdb_string()
    ├── アクセサ:         atom_count(), get_symbol/x/y/z(), get_atom_name(),
    │                    get_residue_name/id(), get_chain_id(), is_hetatm()
    ├── 一括エクスポート: get_positions_flat() → Float32Array, get_symbols_json()
    ├── ボンド検出:       compute_bonds(), get_bonds(), bond_count(), has_bonds_computed()
    ├── 空間クエリ:       distance(), get_atoms_within_radius(), get_residues_within_radius()
    ├── 空間インデックス: build_spatial_index(), has_spatial_index()
    ├── JSON 出力:        get_atom_info() → JsValue, get_neighbors_info() → JsValue
    └── エディタ:         add/remove_atom, add/remove_bond, set_atom_*, set_bond_order,
                          implicit_h_count, add_ring_template, attach_ring_to_bond,
                          get_bounds, rotate_atoms, flip_*, select_atoms_in_rect,
                          move_atoms, check_valence, copy_atoms
```

**設計上の主な判断：**

- **パーサーと Wasm 境界を分離** — 純 Rust 関数（`parse_xyz`, `parse_pdb`）には `JsValue` を含めず、ブラウザなしで `cargo test` による完全なテストを実現
- **座標ベクタの分離** — `x`, `y`, `z` を独立した `Vec<f32>` として保持することで、ベクタライズされた距離計算においてキャッシュ局所性を最大化
- **オブジェクトではなく文字列を渡す** — ファイル内容全体を単一の `&str` として JS/Wasm 境界を越えて渡し、パースはすべて Rust 側で行う
- **オンデマンド計算** — ボンドと空間グリッドは遅延評価。呼び出し側がコストとタイミングを制御できる
- **ボクセルグリッドの自動フォールバック** — `get_atoms_within_radius` はグリッドが構築済みなら平均 O(1)、`build_spatial_index` 未呼び出し時は O(N) 線形探索にフォールバック

## ライセンス

MIT
