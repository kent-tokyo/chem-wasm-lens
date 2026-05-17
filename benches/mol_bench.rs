use criterion::{criterion_group, criterion_main, BatchSize, Criterion};
use chem_wasm_lens::{parse_pdb, parse_xyz};

// ── Synthetic fixtures ────────────────────────────────────────────────────────
//
// Atoms are distributed in 3D (50 Å cube, ~0.08 atoms/Å³ ≈ real protein
// density) so the voxel grid covers all cells meaningfully.
// Distribution: LCG-like scatter via distinct prime-ratio strides.

fn make_pdb(n: usize) -> String {
    let mut s = String::with_capacity(n * 82);
    for i in 0..n {
        let x = (i as f32 * 0.93) % 50.0;
        let y = (i as f32 * 1.71) % 50.0;
        let z = (i as f32 * 2.37) % 50.0;
        let serial = (i + 1) % 100_000;
        let res    = (i / 10 + 1) % 10_000;
        // Fixed-width PDB ATOM line; element column at cols 76-77 = " C"
        s.push_str(&format!(
            "ATOM  {serial:5}  CA  ALA A{res:4}    {x:8.3}{y:8.3}{z:8.3}  1.00  0.00           C  \n"
        ));
    }
    s
}

fn make_xyz(n: usize) -> String {
    let mut s = format!("{n}\nbench\n");
    s.reserve(n * 25);
    for i in 0..n {
        let x = (i as f32 * 0.93) % 50.0;
        let y = (i as f32 * 1.71) % 50.0;
        let z = (i as f32 * 2.37) % 50.0;
        s.push_str(&format!("C {x:.3} {y:.3} {z:.3}\n"));
    }
    s
}

// ── Parse benchmarks ──────────────────────────────────────────────────────────

fn bench_parse_pdb_10k(c: &mut Criterion) {
    let data = make_pdb(10_000);
    c.bench_function("parse_pdb_10k", |b| {
        b.iter(|| parse_pdb(&data).unwrap())
    });
}

fn bench_parse_xyz_10k(c: &mut Criterion) {
    let data = make_xyz(10_000);
    c.bench_function("parse_xyz_10k", |b| {
        b.iter(|| parse_xyz(&data).unwrap())
    });
}

// ── Spatial index ─────────────────────────────────────────────────────────────

fn bench_build_spatial_index_10k(c: &mut Criterion) {
    let data = make_pdb(10_000);
    c.bench_function("build_spatial_index_10k", |b| {
        b.iter_batched(
            || parse_pdb(&data).unwrap(),
            |mut mol| mol.build_spatial_index(5.0),
            BatchSize::SmallInput,
        )
    });
}

// ── Radius query: linear O(N) vs. grid O(1 avg) ───────────────────────────────

fn bench_radius_query_linear_10k(c: &mut Criterion) {
    let data = make_pdb(10_000);
    let mol = parse_pdb(&data).unwrap(); // no spatial index
    c.bench_function("radius_query_linear_10k", |b| {
        b.iter(|| mol.get_atoms_within_radius(5_000, 5.0))
    });
}

fn bench_radius_query_grid_10k(c: &mut Criterion) {
    let data = make_pdb(10_000);
    let mut mol = parse_pdb(&data).unwrap();
    mol.build_spatial_index(5.0);
    c.bench_function("radius_query_grid_10k", |b| {
        b.iter(|| mol.get_atoms_within_radius(5_000, 5.0))
    });
}

// ── Bond detection (O(N²) — capped at 500 atoms) ─────────────────────────────

fn bench_compute_bonds_500(c: &mut Criterion) {
    let data = make_xyz(500);
    c.bench_function("compute_bonds_500", |b| {
        b.iter_batched(
            || parse_xyz(&data).unwrap(),
            |mut mol| mol.compute_bonds(),
            BatchSize::SmallInput,
        )
    });
}

// ── Bond detection: sequential vs parallel comparison ────────────────────────

#[cfg(feature = "parallel")]
fn bench_bonds_comparison(c: &mut Criterion) {
    use criterion::BenchmarkId;

    let mut group = c.benchmark_group("bonds");
    for n in [500usize, 1_000, 2_000] {
        let data = make_xyz(n);
        group.bench_with_input(BenchmarkId::new("sequential", n), &data, |b, d| {
            b.iter_batched(
                || parse_xyz(d).unwrap(),
                |mut mol| mol.compute_bonds(),
                BatchSize::SmallInput,
            )
        });
        group.bench_with_input(BenchmarkId::new("parallel", n), &data, |b, d| {
            b.iter_batched(
                || parse_xyz(d).unwrap(),
                |mut mol| mol.compute_bonds_parallel(),
                BatchSize::SmallInput,
            )
        });
    }
    group.finish();
}

// ── Registration ──────────────────────────────────────────────────────────────

#[cfg(not(feature = "parallel"))]
criterion_group!(
    benches,
    bench_parse_pdb_10k,
    bench_parse_xyz_10k,
    bench_build_spatial_index_10k,
    bench_radius_query_linear_10k,
    bench_radius_query_grid_10k,
    bench_compute_bonds_500,
);

#[cfg(feature = "parallel")]
criterion_group!(
    benches,
    bench_parse_pdb_10k,
    bench_parse_xyz_10k,
    bench_build_spatial_index_10k,
    bench_radius_query_linear_10k,
    bench_radius_query_grid_10k,
    bench_compute_bonds_500,
    bench_bonds_comparison,
);

criterion_main!(benches);
