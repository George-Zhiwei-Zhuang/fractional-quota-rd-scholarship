"""Regenerate all deliverables in a clean copy and compare file hashes.

The same installed Python environment is used. This checks the package's
self-contained workflow, not dependency installation on a new machine.
"""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[2]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    paths=[]
    for directory,pattern in [("output/tables","*"),("output/figures","*"),("data/raw","*"),("docs","*")]:
        paths.extend(p for p in (ROOT/directory).glob(pattern) if p.is_file() and not p.name.startswith("."))
    paths.extend([ROOT/"output/Scholarship_Programme_Figure_Guide.pdf",ROOT/"README.md",ROOT/"environment/versions.json",ROOT/"output/logs/validation.json"])
    expected={str(p.relative_to(ROOT)):digest(p) for p in sorted(paths)}
    with tempfile.TemporaryDirectory(prefix="scholarship_programme_clean_") as temporary:
        clean=Path(temporary)/"Scholarship_Programme_Replication"
        shutil.copytree(ROOT,clean,ignore=shutil.ignore_patterns("output","data","__pycache__","SHA256SUMS.txt"))
        # Run from outside the package to detect accidental dependence on cwd.
        subprocess.run([sys.executable,str(clean/"run.py")],cwd=temporary,check=True)
        rows=[]
        for relative,wanted in expected.items():
            got=digest(clean/relative)
            rows.append(dict(file=relative,match=got==wanted,sha256=wanted))
    result=dict(check="Full workflow in a clean source copy, from an external working directory",
                environment="Same installed environment; see environment/versions.json",
                passed=all(row["match"] for row in rows),files_checked=len(rows),results=rows)
    (ROOT/"output/logs/clean_rerun_check.json").write_text(json.dumps(result,indent=2)+"\n")
    if not result["passed"]:
        raise AssertionError([row["file"] for row in rows if not row["match"]])
    print(f"Clean rerun: {len(rows)} generated files match byte for byte.",flush=True)


if __name__=="__main__":main()
