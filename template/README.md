# chem-wasm-lens React/Vite starter

Get a working molecular viewer in 3 commands:

    npm install
    npm run dev
    # open http://localhost:5173

## Usage

Edit `src/App.tsx` to add your own SMILES strings or integrate
the `MolViewer` component into your existing app.

## Key patterns

`useChem()` hook — singleton `loadChem`, avoids re-initializing Wasm on re-renders.

`MolViewer` — accepts a `smiles` prop and renders an SVG via DOMParser (no innerHTML).

## Why `optimizeDeps.exclude`?

`vite.config.ts` excludes `@kent-tokyo/chem-wasm-lens` from Vite's
pre-bundler because the package contains a `.wasm` file that must
be handled at runtime, not at build time.

## Scripts

| Command | Description |
|---------|-------------|
| `npm run dev` | Start dev server at http://localhost:5173 |
| `npm run build` | TypeScript check + production build |
| `npm run preview` | Preview production build locally |
