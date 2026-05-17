import { loadChem } from '@kent-tokyo/chem-wasm-lens';
import type { MolecularSystem } from '@kent-tokyo/chem-wasm-lens';
import { useEffect, useState } from 'react';

type ChemExports = { MolecularSystem: typeof MolecularSystem };

let _promise: Promise<ChemExports> | null = null;

export function useChem() {
  const [chem, setChem] = useState<ChemExports | null>(null);

  useEffect(() => {
    if (!_promise) _promise = loadChem();
    const p = _promise;
    p.then(setChem);
  }, []);

  return chem;
}
