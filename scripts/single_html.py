#!/usr/bin/env python3
"""Pack the web build (dist/web) into one self-contained HTML file: docs/index.html.

The JavaScript glue is inlined and the WebAssembly module is embedded as base64,
so the file runs from a web server or straight from disk with no other files.
"""
import base64
import re
import sys
from pathlib import Path

root = Path(__file__).resolve().parent.parent
web = root / "dist" / "web"
js_files = sorted(web.glob("ticalc-*.js"))
wasm_files = sorted(web.glob("ticalc-*_bg.wasm"))
if len(js_files) != 1 or len(wasm_files) != 1:
    sys.exit("run scripts/build-dist.sh first: expected one .js and one .wasm in dist/web")

js = js_files[0].read_text()
js, n = re.subn(r"^export \{[^}]*\};?\s*$", "", js, flags=re.M)
if n != 1:
    sys.exit("unexpected export layout in the generated JavaScript")
js = js.replace("</script", "<\\/script")
wasm = base64.b64encode(wasm_files[0].read_bytes()).decode()

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>TI-Calc CAS</title>
<style>
html, body {{ margin: 0; height: 100%; overflow: hidden; background: #1b1b1b; color: #ccc; font-family: sans-serif; }}
canvas {{ width: 100%; height: 100%; display: block; }}
#loading {{ position: absolute; top: 45%; width: 100%; text-align: center; }}
</style>
</head>
<body>
<div id="loading">Loading calculator...</div>
<canvas id="ticalc_canvas"></canvas>
<script type="module">
{js}
const b64 = "{wasm}";
const bin = atob(b64);
const bytes = new Uint8Array(bin.length);
for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
try {{
    await __wbg_init({{ module_or_path: bytes }});
    document.getElementById("loading").remove();
}} catch (err) {{
    document.getElementById("loading").textContent = "Could not start: " + err;
}}
</script>
</body>
</html>
"""
out = root / "docs" / "index.html"
out.write_text(html)
print(f"{out.relative_to(root)}: {out.stat().st_size / 1e6:.1f} MB")
