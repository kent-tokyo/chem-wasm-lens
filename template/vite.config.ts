import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';

export default defineConfig({
  plugins: [react()],
  optimizeDeps: {
    exclude: ['@kent-tokyo/chem-wasm-lens'],
  },
  server: {
    fs: {
      // allow serving files from one level up (for local file: link in dev)
      allow: ['..'],
    },
  },
});
