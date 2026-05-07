#!/usr/bin/env python3
"""Strip notebook outputs on commit and block obvious credential tokens."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

SENSITIVE_PATTERNS = [
    re.compile(r"akid[0-9a-z]{16,}", re.IGNORECASE),
    re.compile(r"secretid", re.IGNORECASE),
    re.compile(r"tencent", re.IGNORECASE),
]


def run_git(args: list[str]) -> str:
    return subprocess.check_output(["git", *args], text=True).strip()


def staged_notebooks() -> list[Path]:
    out = run_git(["diff", "--cached", "--name-only", "--diff-filter=ACMR"])
    if not out:
        return []
    files = [Path(line) for line in out.splitlines() if line.endswith(".ipynb")]
    return [f for f in files if f.exists()]


def stringify_outputs(cell: dict) -> str:
    chunks: list[str] = []
    for output in cell.get("outputs", []):
        text = output.get("text")
        if isinstance(text, list):
            chunks.append("".join(str(x) for x in text))
        elif isinstance(text, str):
            chunks.append(text)

        data = output.get("data", {})
        if isinstance(data, dict):
            for value in data.values():
                if isinstance(value, list):
                    chunks.append("".join(str(x) for x in value))
                else:
                    chunks.append(str(value))
    return "\n".join(chunks)


def has_sensitive(text: str) -> bool:
    return any(p.search(text) for p in SENSITIVE_PATTERNS)


def clean_notebook(path: Path) -> tuple[bool, bool]:
    raw = path.read_text(encoding="utf-8")
    nb = json.loads(raw)

    changed = False
    blocked = False

    for cell in nb.get("cells", []):
        source = "".join(cell.get("source", []))
        if has_sensitive(source):
            blocked = True

        if cell.get("cell_type") == "code":
            out_text = stringify_outputs(cell)
            if has_sensitive(out_text):
                changed = True
            if cell.get("outputs"):
                cell["outputs"] = []
                changed = True
            if cell.get("execution_count") is not None:
                cell["execution_count"] = None
                changed = True

    if changed:
        new_raw = json.dumps(nb, ensure_ascii=False, indent=1) + "\n"
        if new_raw != raw:
            path.write_text(new_raw, encoding="utf-8")

    return changed, blocked


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="*", help="Notebook files to clean")
    parser.add_argument("--staged", action="store_true", help="Use staged notebook files")
    args = parser.parse_args()

    if args.staged:
        targets = staged_notebooks()
    else:
        targets = [Path(f) for f in args.files]

    if not targets:
        return 0

    modified: list[Path] = []
    blocked: list[Path] = []

    for path in targets:
        changed, has_source_secret = clean_notebook(path)
        if changed:
            modified.append(path)
        if has_source_secret:
            blocked.append(path)

    if modified:
        subprocess.check_call(["git", "add", *[str(p) for p in modified]])
        print("Sanitized notebook outputs:")
        for p in modified:
            print(f"  - {p}")

    if blocked:
        print("\nCommit blocked: sensitive token pattern found in notebook source cell(s).")
        for p in blocked:
            print(f"  - {p}")
        print("Remove credentials from source cells and retry commit.")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
