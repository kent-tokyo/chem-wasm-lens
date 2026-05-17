import { useEffect, useRef } from 'react';
import { useChem } from '../hooks/useChem';

interface Props {
  smiles: string;
  width?: number;
  height?: number;
}

export function MolViewer({ smiles, width = 400, height = 300 }: Props) {
  const chem = useChem();
  const containerRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    if (!chem || !containerRef.current) return;
    const el = containerRef.current;
    while (el.firstChild) el.removeChild(el.firstChild);
    try {
      const mol = chem.MolecularSystem.from_smiles(smiles);
      mol.compute_2d_coords();
      const svg = mol.to_svg_string(width, height);
      const svgEl = new DOMParser()
        .parseFromString(svg, 'image/svg+xml').documentElement;
      el.appendChild(svgEl);
    } catch {
      el.textContent = 'Invalid SMILES';
    }
  }, [chem, smiles, width, height]);

  if (!chem) {
    return (
      <div style={{ width, height, display: 'flex', alignItems: 'center', justifyContent: 'center' }}>
        Loading...
      </div>
    );
  }
  return <div ref={containerRef} />;
}
