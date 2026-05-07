#!/usr/bin/env python3
"""Notebook secret sanitization utilities.

This tool intentionally preserves notebook outputs and execution counts.
It only redacts/blocks sensitive token patterns.

Modes:
  --sanitize-secrets-staged     Redact secrets in staged notebooks and re-stage them.
  --check-push-secrets REMOTE   Check pushed commit range(s) for secrets; block push if found.
  [files]                       Redact secrets in the given notebook files in-place.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

NULL_SHA = "0" * 40

SENSITIVE_PATTERNS = [
    re.compile(r"akid[0-9a-z]{16,}", re.IGNORECASE),
    re.compile(r"secretid", re.IGNORECASE),
    re.compile(r"tencent", re.IGNORECASE),
]

REDACTED = "[REDACTED_SECRET]"


def run_git(args: list[str]) -> str:
    return subprocess.check_output(["git", *args], text=True).strip()


def run_git_bytes(args: list[str]) -> bytes:
    return subprocess.check_output(["git", *args])


def has_sensitive(text: str) -> bool:
    return any(p.search(text) for p in SENSITIVE_PATTERNS)


def redact_sensitive(text: str) -> str:
    out = text
    for pat in SENSITIVE_PATTERNS:
        out = pat.sub(REDACTED, out)
    return out


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


def redact_value(value):
    if isinstance(value, str):
        return redact_sensitive(value)
    if isinstance(value, list):
        return [redact_value(v) for v in value]
    if isinstance(value, dict):
        return {k: redact_value(v) for k, v in value.items()}
    return value


def sanitize_notebook_object(nb: dict) -> tuple[dict, bool, bool]:
    """Return (notebook, changed, still_sensitive_after_sanitize)."""
    changed = False

    for cell in nb.get("cells", []):
        src = cell.get("source")
        new_src = redact_value(src)
        if new_src != src:
            cell["source"] = new_src
            changed = True

        if cell.get("cell_type") == "code":
            outs = cell.get("outputs", [])
            new_outs = redact_value(outs)
            if new_outs != outs:
                cell["outputs"] = new_outs
                changed = True

    still_sensitive = False
    for cell in nb.get("cells", []):
        source = "".join(cell.get("source", []))
        if has_sensitive(source):
            still_sensitive = True
            break
        if cell.get("cell_type") == "code" and has_sensitive(stringify_outputs(cell)):
            still_sensitive = True
            break

    return nb, changed, still_sensitive


def sanitize_notebook_file(path: Path) -> tuple[bool, bool]:
    """Sanitize a notebook file in place. Returns (changed, blocked)."""
    raw = path.read_text(encoding="utf-8")
    nb = json.loads(raw)
    nb, changed, blocked = sanitize_notebook_object(nb)

    if changed:
        path.write_text(json.dumps(nb, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")

    return changed, blocked


def staged_notebooks() -> list[Path]:
    out = run_git(["diff", "--cached", "--name-only", "--diff-filter=ACMR"])
    if not out:
        return []
    files = [Path(line) for line in out.splitlines() if line.endswith(".ipynb")]
    return [f for f in files if f.exists()]


def sanitize_secrets_staged() -> int:
    targets = staged_notebooks()
    if not targets:
        return 0

    modified: list[Path] = []
    blocked: list[Path] = []

    for path in targets:
        try:
            changed, still_sensitive = sanitize_notebook_file(path)
        except Exception:
            # If parsing fails, do not silently pass a potentially unsafe notebook.
            blocked.append(path)
            continue
        if changed:
            modified.append(path)
        if still_sensitive:
            blocked.append(path)

    if modified:
        subprocess.check_call(["git", "add", *[str(p) for p in modified]])
        print("Sanitized notebook secrets:")
        for p in modified:
            print(f"  - {p}")

    if blocked:
        print("\nCommit blocked: sensitive token pattern still present in notebook(s).")
        for p in blocked:
            print(f"  - {p}")
        print("Remove credentials from source/output cells and retry commit.")
        return 1

    return 0


def check_blob_for_secret(blob_sha: str) -> bool:
    content = run_git_bytes(["cat-file", "blob", blob_sha])
    try:
        nb = json.loads(content.decode("utf-8"))
    except Exception:
        return False

    for cell in nb.get("cells", []):
        source = "".join(cell.get("source", []))
        if has_sensitive(source):
            return True
        if cell.get("cell_type") == "code" and has_sensitive(stringify_outputs(cell)):
            return True

    return False


def check_push_secrets(stdin_text: str) -> int:
    """Inspect notebook blobs in pushed commit ranges; block push on any secret pattern."""
    refs: list[tuple[str, str, str, str]] = []
    for line in stdin_text.strip().splitlines():
        parts = line.strip().split()
        if len(parts) == 4:
            refs.append((parts[0], parts[1], parts[2], parts[3]))

    if not refs:
        return 0

    offenders: set[str] = set()

    for _local_ref, local_sha, _remote_ref, remote_sha in refs:
        if local_sha == NULL_SHA:
            continue

        rev_range = local_sha if remote_sha == NULL_SHA else f"{remote_sha}..{local_sha}"
        commits_out = run_git(["rev-list", rev_range])
        if not commits_out:
            continue

        for commit_sha in commits_out.splitlines():
            names_out = run_git(["diff-tree", "--no-commit-id", "--name-only", "-r", commit_sha])
            nb_paths = [p for p in names_out.splitlines() if p.endswith(".ipynb")]
            for nb_path in nb_paths:
                try:
                    blob_sha = run_git(["rev-parse", f"{commit_sha}:{nb_path}"])
                except subprocess.CalledProcessError:
                    continue
                if check_blob_for_secret(blob_sha):
                    offenders.add(nb_path)

    if offenders:
        print("\nPush blocked: sensitive token pattern found in notebook content.")
        for p in sorted(offenders):
            print(f"  - {p}")
        print("Run commit again to sanitize, or manually redact these cells.")
        return 1

    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="*", help="Notebook files to sanitize secrets in-place")
    parser.add_argument(
        "--sanitize-secrets-staged", action="store_true",
        help="(pre-commit) Redact sensitive patterns from staged notebooks and re-stage.",
    )
    parser.add_argument(
        "--check-push-secrets", metavar="REMOTE",
        help="(pre-push) Check pushed commit ranges for notebook secrets; refs read from stdin.",
    )
    args = parser.parse_args()

    if args.sanitize_secrets_staged:
        return sanitize_secrets_staged()

    if args.check_push_secrets is not None:
        stdin_text = sys.stdin.read()
        return check_push_secrets(stdin_text)

    targets = [Path(f) for f in args.files]
    if not targets:
        return 0

    blocked: list[Path] = []
    for path in targets:
        changed, still_sensitive = sanitize_notebook_file(path)
        if changed:
            print(f"Sanitized notebook secrets: {path}")
        if still_sensitive:
            blocked.append(path)

    if blocked:
        print("\nSensitive token pattern still present in notebook(s):")
        for p in blocked:
            print(f"  - {p}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
