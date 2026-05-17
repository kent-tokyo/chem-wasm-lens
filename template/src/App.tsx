import { useState } from 'react';
import { MolViewer } from './components/MolViewer';

const PRESETS = [
  { label: 'Ethanol', smiles: 'CCO' },
  { label: 'Benzene', smiles: 'c1ccccc1' },
  { label: 'Aspirin', smiles: 'CC(=O)Oc1ccccc1C(=O)O' },
  { label: 'Caffeine', smiles: 'Cn1cnc2c1c(=O)n(c(=O)n2C)C' },
];

export default function App() {
  const [smiles, setSmiles] = useState('c1ccccc1');

  return (
    <div style={{ fontFamily: 'monospace', maxWidth: 600, margin: '2rem auto', padding: '0 1rem' }}>
      <h1 style={{ fontSize: '1.1rem' }}>chem-wasm-lens — React starter</h1>
      <div style={{ marginBottom: '0.5rem' }}>
        {PRESETS.map(p => (
          <button key={p.smiles} onClick={() => setSmiles(p.smiles)} style={{ marginRight: 6 }}>
            {p.label}
          </button>
        ))}
      </div>
      <input
        value={smiles}
        onChange={e => setSmiles(e.target.value)}
        style={{ width: '100%', fontFamily: 'monospace', padding: '0.3rem', boxSizing: 'border-box' }}
        placeholder="Enter SMILES..."
      />
      <MolViewer smiles={smiles} width={560} height={400} />
    </div>
  );
}
