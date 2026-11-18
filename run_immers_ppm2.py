#!/usr/bin/env python3
"""Run immers in an isolated short-name work directory."""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_IMMERS = SCRIPT_DIR / "immers"
TABLE_OUTPUTS = {"datapar1", "datasub1"}
DEFAULT_HETERO = 0


def make_short_stem(index: int) -> str:
    alphabet = "0123456789abcdefghijklmnopqrstuvwxyz"
    value = index
    digits = []
    for _ in range(3):
        digits.append(alphabet[value % len(alphabet)])
        value //= len(alphabet)
    if value:
        raise ValueError("too many PDB files for short temporary names")
    return "p" + "".join(reversed(digits))


def render_input(hetero: int, topology: str, pdb_name: str) -> str:
    #ppm2 format: " 0 out filename.pdb"  (no membrane code, only in ppm3 this exists)
    return f"{hetero:2d} {topology:3s} {pdb_name}\n"


def rewrite_table_names(path: Path, short_name: str, original_name: str) -> None:
    text = path.read_text(errors="replace")
    text = text.replace(short_name, original_name)
    text = text.replace(Path(short_name).stem, Path(original_name).stem)
    path.write_text(text)


def run_one(pdb_path: Path, topology: str, immers: Path, res_lib: Path, debug: bool) -> list[Path]:
    pdb_path = pdb_path.resolve()
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    output_dir = pdb_path.parent
    short_name = f"{make_short_stem(1)}.pdb"

    #write .inp next to the original PDB for reference
    generated_input = output_dir / f"{pdb_path.stem}.inp"
    generated_input.write_text(render_input(DEFAULT_HETERO, topology, pdb_path.name))

    temp_dir = Path(tempfile.mkdtemp(prefix="opm2-"))
    try:
        if debug:
            print(f"[debug] temp dir: {temp_dir}", file=sys.stderr)

        #fill temp dir
        shutil.copy2(pdb_path, temp_dir / short_name)
        try:
            (temp_dir / "immers").symlink_to(immers.resolve())
        except OSError:
            shutil.copy2(immers.resolve(), temp_dir / "immers")
        try:
            (temp_dir / "res.lib").symlink_to(res_lib.resolve())
        except OSError:
            shutil.copy2(res_lib.resolve(), temp_dir / "res.lib")

        rewritten_input = temp_dir / "input.inp"
        rewritten_input.write_text(render_input(DEFAULT_HETERO, topology, short_name))

        initial_names = {p.name for p in temp_dir.iterdir()}

        stdout_path = temp_dir / "immers.stdout"
        if debug:
            print("[debug] running: ./immers < input.inp", file=sys.stderr)

        with rewritten_input.open("rb") as stdin, stdout_path.open("wb") as stdout:
            completed = subprocess.run(
                ["./immers"],
                cwd=temp_dir,
                stdin=stdin,
                stdout=stdout,
                stderr=subprocess.STDOUT,
                check=False,
            )

        if debug:
            print(f"[debug] exit code: {completed.returncode}", file=sys.stderr)

        #collect and rename outputs
        output_dir.mkdir(parents=True, exist_ok=True)
        copied: list[Path] = [generated_input]

        stdout_dest = output_dir / f"{pdb_path.stem}.opm2.stdout"
        shutil.copy2(stdout_path, stdout_dest)
        copied.append(stdout_dest)

        short_stem = Path(short_name).stem  # e.g. "p001"
        orig_stem  = pdb_path.stem          # e.g. "c0"

        for path in sorted(temp_dir.iterdir()):
            if path.name in initial_names or path.name == stdout_path.name or not path.is_file():
                continue

            #rename to not rewrtie files p001out.pdb -> c0out.pdb,  datapar1 -> c0.datapar1
            if path.name in TABLE_OUTPUTS:
                dest = output_dir / f"{orig_stem}.{path.name}"
            else:
                dest = output_dir / path.name.replace(short_stem, orig_stem)

            shutil.copy2(path, dest)
            if path.name in TABLE_OUTPUTS:
                rewrite_table_names(dest, short_name, pdb_path.name)
            copied.append(dest)

        if completed.returncode != 0:
            raise RuntimeError(
                f"immers (ppm2) failed for {pdb_path} with exit code {completed.returncode}"
            )
        return copied

    finally:
        if not debug:
            shutil.rmtree(temp_dir)


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description="Run ppm2 immers on a single PDB file.")
    parser.add_argument("pdb_file", type=Path, help="Input PDB file.")
    parser.add_argument("topology", choices=["in", "out"], help="N-terminus topology: in or out.")
    parser.add_argument("--immers", type=Path, default=DEFAULT_IMMERS, help="Path to immers executable.")
    parser.add_argument("--res-lib", type=Path, default=None, help="Path to res.lib.")
    parser.add_argument("--debug", action="store_true", help="Keep temp dir and print debug info.")
    args = parser.parse_args()

    immers  = args.immers.resolve()
    res_lib = (args.res_lib or immers.parent / "res.lib").resolve()

    if not immers.exists():
        print(f"immers executable not found: {immers}", file=sys.stderr); sys.exit(1)
    if not os.access(immers, os.X_OK):
        print(f"immers is not executable: {immers}", file=sys.stderr); sys.exit(1)
    if not res_lib.exists():
        print(f"res.lib not found: {res_lib}", file=sys.stderr); sys.exit(1)

    try:
        outputs = run_one(args.pdb_file, args.topology.lower(), immers, res_lib, args.debug)
        for o in outputs:
            print(o)
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
