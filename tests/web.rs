use wasm_bindgen::JsValue;
use wasm_bindgen_test::*;
use chem_wasm_lens::MolecularSystem;

wasm_bindgen_test_configure!(run_in_browser);

const WATER_XYZ: &str = "3\nwater\nO  0.000  0.000  0.119\nH  0.000  0.757 -0.477\nH  0.000 -0.757 -0.477\n";
const ATOM_N_GLY: &str =
    "ATOM      1  N   GLY A   1       1.000   2.000   3.000  1.00  0.00           N  ";
const THREE_CARBONS: &str =
    "3\ntest\nC  0.000  0.000  0.000\nC  1.000  0.000  0.000\nC  5.000  0.000  0.000\n";

// ── from_xyz_string: Result<_, JsValue> boundary ─────────────────────────────

#[wasm_bindgen_test]
fn from_xyz_string_ok_parses_atom_count() {
    let mol = MolecularSystem::from_xyz_string(WATER_XYZ).unwrap();
    assert_eq!(mol.atom_count(), 3);
    assert_eq!(mol.get_symbol(0), Some("O".to_string()));
}

#[wasm_bindgen_test]
fn from_xyz_string_err_is_js_string() {
    let err = MolecularSystem::from_xyz_string("").unwrap_err();
    assert!(err.is_string());
    assert!(err.as_string().unwrap().contains("empty input"));
}

// ── from_pdb_string: Result<_, JsValue> boundary ─────────────────────────────

#[wasm_bindgen_test]
fn from_pdb_string_ok_parses_fields() {
    let mol = MolecularSystem::from_pdb_string(ATOM_N_GLY).unwrap();
    assert_eq!(mol.atom_count(), 1);
    assert_eq!(mol.get_symbol(0), Some("N".to_string()));
}

#[wasm_bindgen_test]
fn from_pdb_string_err_is_js_string() {
    let err = MolecularSystem::from_pdb_string("").unwrap_err();
    assert!(err.is_string());
    assert!(err.as_string().unwrap().contains("empty input"));
}

// ── get_positions_flat → Vec<f32> (becomes Float32Array on JS side) ───────────

#[wasm_bindgen_test]
fn positions_flat_correct_length_and_values() {
    let mol = MolecularSystem::from_xyz_string(WATER_XYZ).unwrap();
    let flat = mol.get_positions_flat();
    assert_eq!(flat.len(), 9); // 3 atoms × 3 coords
    assert!((flat[2] - 0.119_f32).abs() < 1e-4); // z of O
    assert!((flat[4] - 0.757_f32).abs() < 1e-4); // y of first H
}

// ── get_atom_info → JsValue (serde_wasm_bindgen) ─────────────────────────────

#[wasm_bindgen_test]
fn get_atom_info_pdb_fields_accessible_via_reflect() {
    let mol = MolecularSystem::from_pdb_string(ATOM_N_GLY).unwrap();
    let info: JsValue = mol.get_atom_info(0);
    assert!(!info.is_null());

    let symbol = js_sys::Reflect::get(&info, &"symbol".into()).unwrap();
    assert_eq!(symbol.as_string().unwrap(), "N");

    let res_name = js_sys::Reflect::get(&info, &"residue_name".into()).unwrap();
    assert_eq!(res_name.as_string().unwrap(), "GLY");

    let is_hetatm = js_sys::Reflect::get(&info, &"is_hetatm".into()).unwrap();
    assert!(!is_hetatm.as_bool().unwrap());
}

#[wasm_bindgen_test]
fn get_atom_info_returns_null_for_out_of_bounds() {
    let mol = MolecularSystem::from_xyz_string(WATER_XYZ).unwrap();
    assert!(mol.get_atom_info(99).is_null());
}

// ── get_neighbors_info → JsValue (JS Array of objects) ───────────────────────

#[wasm_bindgen_test]
fn get_neighbors_info_returns_array_with_correct_length() {
    let mol = MolecularSystem::from_xyz_string(THREE_CARBONS).unwrap();
    let neighbors = mol.get_neighbors_info(0, 2.0);
    let arr = js_sys::Array::from(&neighbors);
    assert_eq!(arr.length(), 1);
    let item = arr.get(0);
    let symbol = js_sys::Reflect::get(&item, &"symbol".into()).unwrap();
    assert_eq!(symbol.as_string().unwrap(), "C");
}

// ── get_atoms_within_radius → Vec<u32> ───────────────────────────────────────

#[wasm_bindgen_test]
fn atoms_within_radius_serializes_indices_correctly() {
    let mol = MolecularSystem::from_xyz_string(THREE_CARBONS).unwrap();
    assert_eq!(mol.get_atoms_within_radius(0, 2.0), vec![1u32]);
    assert_eq!(mol.get_atoms_within_radius(0, 6.0), vec![1u32, 2u32]);
}

// ── build_spatial_index → state + query end-to-end ───────────────────────────

#[wasm_bindgen_test]
fn spatial_index_state_and_query_end_to_end() {
    let mut mol = MolecularSystem::from_xyz_string(THREE_CARBONS).unwrap();
    assert!(!mol.has_spatial_index());
    mol.build_spatial_index(3.0);
    assert!(mol.has_spatial_index());
    assert_eq!(mol.get_atoms_within_radius(0, 2.0), vec![1u32]);
}

// ── compute_bonds + get_bonds → Vec<u32> ─────────────────────────────────────

#[wasm_bindgen_test]
fn bonds_round_trip_through_wasm_boundary() {
    let mut mol = MolecularSystem::from_xyz_string(
        "2\ntest\nH  0.000  0.000  0.000\nH  0.740  0.000  0.000\n",
    )
    .unwrap();
    assert!(!mol.has_bonds_computed());
    mol.compute_bonds();
    assert!(mol.has_bonds_computed());
    assert_eq!(mol.bond_count(), 1);
    assert_eq!(mol.get_bonds(0), vec![1u32]);
    assert_eq!(mol.get_bonds(1), vec![0u32]);
}
