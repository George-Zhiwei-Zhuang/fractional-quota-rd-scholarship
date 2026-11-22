"""Run from any working directory: python /path/to/package/run.py."""
from pathlib import Path
import argparse
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser()
p.add_argument("--figures-only",action="store_true",help="Use the included simulation tables.")
p.add_argument("--simulation-only",action="store_true",help="Regenerate numerical outputs only.")
args=p.parse_args()
if args.figures_only and args.simulation_only:
    p.error("Choose at most one partial-run option.")
if not args.figures_only:
    subprocess.run([sys.executable,str(ROOT/"code/03_analyze/simulate.py")],check=True,cwd=ROOT)
if not args.simulation_only:
    subprocess.run([sys.executable,str(ROOT/"code/05_figures/make_figures.py")],check=True,cwd=ROOT)
    subprocess.run([sys.executable,str(ROOT/"code/04_tables/make_report.py")],check=True,cwd=ROOT)
