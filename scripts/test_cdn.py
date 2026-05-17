"""
CDN bundle smoke test via Playwright.

Serves pkg/dist/chem-wasm-lens.esm.js on a local HTTP server,
loads it in a real browser (Chromium), and verifies the public API.

Usage:
    python3 scripts/test_cdn.py
"""

import asyncio
import threading
import http.server
import os
from pathlib import Path
from playwright.async_api import async_playwright

ROOT = Path(__file__).parent.parent
PORT = 18765

# --------------------------------------------------------------------------- #
# Minimal HTTP server for the test
# --------------------------------------------------------------------------- #

class QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *_):
        pass

def start_server():
    os.chdir(ROOT)
    server = http.server.HTTPServer(("127.0.0.1", PORT), QuietHandler)
    t = threading.Thread(target=server.serve_forever, daemon=True)
    t.start()
    return server

# --------------------------------------------------------------------------- #
# Inline test page (served as data URL via page.goto)
# --------------------------------------------------------------------------- #

TEST_HTML = f"""<!DOCTYPE html>
<html>
<head><meta charset="UTF-8"></head>
<body>
<script type="module">
import {{ loadChem }} from 'http://127.0.0.1:{PORT}/pkg/dist/chem-wasm-lens.esm.js';

async function runTests() {{
  const results = [];

  function ok(label, cond) {{
    results.push({{ label, pass: Boolean(cond) }});
  }}

  try {{
    const {{ MolecularSystem }} = await loadChem();

    // 1. from_smiles — atom count (ethanol CCO → 2C+1O+6H = 9)
    const eth = MolecularSystem.from_smiles('CCO');
    ok('from_smiles atom_count CCO', eth.atom_count() === 9);

    // 2. from_smiles — bond count
    ok('from_smiles bond_count CCO', eth.bond_count() === 8);

    // 3. compute_2d_coords — positions change from 0
    eth.compute_2d_coords();
    const pos = eth.get_positions_flat();
    const nonzero = Array.from(pos).some(v => Math.abs(v) > 0.01);
    ok('compute_2d_coords sets nonzero positions', nonzero);

    // 4. to_svg_string — returns valid SVG
    const svg = eth.to_svg_string(400, 300);
    ok('to_svg_string returns SVG string', typeof svg === 'string' && svg.includes('<svg'));
    ok('to_svg_string contains oxygen label', svg.includes('>O<'));
    ok('to_svg_string no carbon label', !svg.includes('>C<'));

    // 5. Benzene aromatic
    const benz = MolecularSystem.from_smiles('c1ccccc1');
    benz.compute_2d_coords();
    ok('benzene is_aromatic(0)', benz.is_aromatic(0) === true);
    const svgBenz = benz.to_svg_string(300, 300);
    ok('benzene SVG has aromatic dashes', svgBenz.includes('stroke-dasharray'));

    // 6. from_pdb_string
    const pdbSnip = [
      'ATOM      1  N   GLY A   1       1.885  22.498   3.903  1.00  0.00           N  ',
      'ATOM      2  CA  GLY A   1       2.849  22.100   2.916  1.00  0.00           C  ',
    ].join('\\n');
    const prot = MolecularSystem.from_pdb_string(pdbSnip);
    ok('from_pdb_string atom_count', prot.atom_count() === 2);

    // 7. loadChem idempotent (second call returns same object)
    const {{ MolecularSystem: MS2 }} = await loadChem();
    ok('loadChem idempotent', MS2 === MolecularSystem);

  }} catch (e) {{
    results.push({{ label: 'EXCEPTION', pass: false, error: String(e) }});
  }}

  window.__TEST_RESULTS__ = results;
}}

runTests();
</script>
</body>
</html>
"""

# --------------------------------------------------------------------------- #
# Main test runner
# --------------------------------------------------------------------------- #

PASS = "\033[32mPASS\033[0m"
FAIL = "\033[31mFAIL\033[0m"

async def main():
    server = start_server()
    print(f"[server] http://127.0.0.1:{PORT}")

    async with async_playwright() as pw:
        browser = await pw.chromium.launch()
        page = await browser.new_page()

        # Write test page and navigate
        test_page_path = ROOT / "pkg" / "dist" / "_test.html"
        test_page_path.write_text(TEST_HTML, encoding="utf-8")
        url = f"http://127.0.0.1:{PORT}/pkg/dist/_test.html"

        await page.goto(url)

        # Wait for results (up to 15 s)
        await page.wait_for_function("window.__TEST_RESULTS__ !== undefined", timeout=15000)
        results = await page.evaluate("window.__TEST_RESULTS__")

        await browser.close()
        test_page_path.unlink(missing_ok=True)

    server.shutdown()

    # Print results
    print()
    passed = 0
    failed = 0
    for r in results:
        status = PASS if r["pass"] else FAIL
        print(f"  {status}  {r['label']}")
        if r["pass"]:
            passed += 1
        else:
            failed += 1
            if "error" in r:
                print(f"         {r['error']}")

    print()
    print(f"  {passed} passed, {failed} failed")
    print()
    if failed:
        raise SystemExit(1)

asyncio.run(main())
