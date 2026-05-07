#!/usr/bin/env python3
"""Notebook sanitization utilities.

Modes:
  --check-credentials        (pre-commit) Block staged notebooks containing credential patterns.
  --pre-push REMOTE          (pre-push)   Rewrite pushed commits to strip outputs, push clean
                                          copies, update local refs, restore working-tree outputs.
  [files]                    (manual)     Strip outputs from the given notebook files in-place.
"""

from __future__ import annotations

import argparse
import json
import os
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


# ---------------------------------------------------------------------------
# Git helpers
# ---------------------------------------------------------------------------

def run_git(args: list[str]) -> str:
    return subprocess.check_output(["git", *args], text=True).strip()


def run_git_bytes(args: list[str]) -> bytes:
    return subprocess.check_output(["git", *args])


# ---------------------------------------------------------------------------
# Credential detection
# ---------------------------------------------------------------------------

def has_sensitive(text: str) -> bool:
    return any(p.search(text) for p in SENSITIVE_PATTERNS)


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


def staged_notebooks() -> list[Path]:
    out = run_git(["diff", "--cached", "--name-only", "--diff-filter=ACMR"])
    if not out:
        return []
    files = [Path(line) for line in out.splitlines() if line.endswith(".ipynb")]
    return [f for f in files if f.exists()]


def check_credentials_staged() -> int:
    """Check staged notebooks for credential patterns. Returns 1 if any found."""
    targets = staged_notebooks()
    if not targets:
        return 0
    blocked: list[Path] = []
    for path in targets:
        try:
            raw = path.read_text(encoding="utf-8")
            nb = json.loads(raw)
        except Exception:
            continue
        for cell in nb.get("cells", []):
            source = "".join(cell.get("source", []))
            if has_sensitive(source):
                blocked.append(path)
                break
            if cell.get("cell_type") == "code":
                if has_sensitive(stringify_outputs(cell)):
                    blocked.append(path)
                    break
    if blocked:
        print("\nCommit blocked: sensitive token pattern found in notebook(s).")
        for p in blocked:
            print(f"  - {p}")
        print("Remove credentials from source cells and retry commit.")
        return 1
    return 0


# ---------------------------------------------------------------------------
# In-place notebook stripping (manual / legacy use)
# ---------------------------------------------------------------------------

def clean_notebook(path: Path) -> tuple[bool, bool]:
    """Strip outputs from a notebook file in-place. Returns (changed, blocked)."""
    raw = path.read_text(encoding="utf-8")
    nb = json.loads(raw)
    changed = False
    blocked = False
    for cell in nb.get("cells", []):
        source = "".join(cell.get("source", []))
        if has_sensitive(source):
            blocked = True
        if cell.get("cell_type") == "code":
            if has_sensitive(stringify_outputs(cell)):
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


# ---------------------------------------------------------------------------
# Git-object-level stripping for pre-push commit rewriting
# ---------------------------------------------------------------------------

def strip_blob(sha: str) -> tuple[str, bool]:
    """Strip outputs from a notebook blob. Returns (new_sha, changed)."""
    content = run_git_bytes(["cat-file", "blob", sha])
    try:
        nb = json.loads(content.decode("utf-8"))
    except Exception:
        return sha, False
    changed = False
    for cell in nb.get("cells", []):
        if cell.get("cell_type") == "code":
            if cell.get("outputs"):
                cell["outputs"] = []
                changed = True
            if cell.get("execution_count") is not None:
                cell["execution_count"] = None
                changed = True
    if not changed:
        return sha, False
    new_content = (json.dumps(nb, ensure_ascii=False, indent=1) + "\n").encode("utf-8")
    new_sha = subprocess.check_output(
        ["git", "hash-object", "-w", "--stdin"], input=new_content
    ).decode("utf-8").strip()
    return new_sha, True


def strip_tree(tree_sha: str, cache: dict[str, tuple[str, bool]]) -> tuple[str, bool]:
    """Recursively strip notebooks from a git tree. Returns (new_tree_sha, changed)."""
    if tree_sha in cache:
        return cache[tree_sha]
    lines = run_git(["ls-tree", tree_sha]).splitlines()
    new_entries: list[str] = []
    changed = False
    for line in lines:
        tab = line.index("\t")
        mode, type_, sha = line[:tab].split()
        name = line[tab + 1:]
        if type_ == "blob" and name.endswith(".ipynb"):
            new_sha, blob_changed = strip_blob(sha)
            if blob_changed:
                sha = new_sha
                changed = True
        elif type_ == "tree":
            new_sha, sub_changed = strip_tree(sha, cache)
            if sub_changed:
                sha = new_sha
                changed = True
        new_entries.append(f"{mode} {type_} {sha}\t{name}")
    if not changed:
        cache[tree_sha] = (tree_sha, False)
        return tree_sha, False
    new_tree = subprocess.check_output(
        ["git", "mktree"], input="\n".join(new_entries).encode("utf-8")
    ).decode("utf-8").strip()
    cache[tree_sha] = (new_tree, True)
    return new_tree, True


def rewrite_commits(base_sha: str, tip_sha: str) -> str:
    """Rewrite base_sha..tip_sha stripping notebook outputs. Returns new tip sha."""
    rev_range = tip_sha if base_sha == NULL_SHA else f"{base_sha}..{tip_sha}"
    commits_out = run_git(["rev-list", "--reverse", rev_range])
    if not commits_out:
        return tip_sha
    commits = commits_out.splitlines()
    mapping: dict[str, str] = {}
    tree_cache: dict[str, tuple[str, bool]] = {}
    for commit_sha in commits:
        tree_sha = run_git(["rev-parse", f"{commit_sha}^{{tree}}"])
        new_tree, _ = strip_tree(tree_sha, tree_cache)
        parents_raw = run_git(["log", "--format=%P", "-1", commit_sha])
        parents = parents_raw.split() if parents_raw else []
        new_parents = [mapping.get(p, p) for p in parents]
        parent_args: list[str] = []
        for p in new_parents:
            parent_args += ["-p", p]
        commit_msg = run_git(["log", "--format=%B", "-1", commit_sha])
        env = {
            **os.environ,
            "GIT_AUTHOR_NAME":     run_git(["log", "--format=%an", "-1", commit_sha]),
            "GIT_AUTHOR_EMAIL":    run_git(["log", "--format=%ae", "-1", commit_sha]),
            "GIT_AUTHOR_DATE":     run_git(["log", "--format=%aI", "-1", commit_sha]),
            "GIT_COMMITTER_NAME":  run_git(["log", "--format=%cn", "-1", commit_sha]),
            "GIT_COMMITTER_EMAIL": run_git(["log", "--format=%ce", "-1", commit_sha]),
            "GIT_COMMITTER_DATE":  run_git(["log", "--format=%cI", "-1", commit_sha]),
        }
        new_commit = subprocess.check_output(
            ["git", "commit-tree", new_tree, *parent_args, "-m", commit_msg], env=env
        ).decode("utf-8").strip()
        mapping[commit_sha] = new_commit
    return mapping.get(tip_sha, tip_sha)


# ---------------------------------------------------------------------------
# Pre-push mode
# ---------------------------------------------------------------------------

def pre_push_mode(remote: str, stdin_text: str) -> int:
    """
    Rewrite pushed commits to strip notebook outputs.
    Pushes clean commits to remote, updates local refs to match, and restores
    working-tree notebook content so outputs remain visible locally.
    Returns 1 after the inner push (to abort the original unclean push), or 0
    if nothing needed stripping (original push proceeds normally).
    """
    git_dir = Path(run_git(["rev-parse", "--git-dir"]))
    lock_path = git_dir / "SANITIZE_PUSH_LOCK"

    if lock_path.exists():
        # Inner push triggered by this very hook; let it through.
        return 0

    refs: list[tuple[str, str, str, str]] = []
    for line in stdin_text.strip().splitlines():
        parts = line.strip().split()
        if len(parts) == 4:
            refs.append((parts[0], parts[1], parts[2], parts[3]))
    if not refs:
        return 0

    # Save working-tree content of all tracked notebooks before any modification.
    tracked_out = run_git(["ls-files", "*.ipynb"])
    tracked_nbs: list[Path] = [Path(p) for p in tracked_out.splitlines() if p and Path(p).exists()]
    saved: dict[Path, bytes] = {p: p.read_bytes() for p in tracked_nbs}

    push_specs: list[str] = []
    local_ref_updates: list[tuple[str, str]] = []
    any_changed = False

    for local_ref, local_sha, remote_ref, remote_sha in refs:
        if local_sha == NULL_SHA:
            push_specs.append(f":{remote_ref}")
            continue
        base = remote_sha if remote_sha != NULL_SHA else NULL_SHA
        new_tip = rewrite_commits(base, local_sha)
        if new_tip != local_sha:
            any_changed = True
        push_specs.append(f"{new_tip}:{remote_ref}")
        local_ref_updates.append((local_ref, new_tip))

    if not any_changed:
        # No notebooks had outputs; let the original push proceed unmodified.
        return 0

    lock_path.touch()
    push_rc = 1
    try:
        push_result = subprocess.run(["git", "push", remote, *push_specs])
        push_rc = push_result.returncode
    finally:
        try:
            lock_path.unlink()
        except FileNotFoundError:
            pass

    if push_rc == 0:
        # Align local refs with the cleaned tips so local and remote agree.
        for local_ref, new_tip in local_ref_updates:
            if local_ref.startswith("refs/"):
                subprocess.run(["git", "update-ref", local_ref, new_tip], check=False)
        # Restore working-tree notebooks so outputs remain visible locally.
        for nb_path, content in saved.items():
            nb_path.write_bytes(content)
        print("pre-push: notebook outputs stripped for remote; working tree restored.")

    # Return non-zero so git aborts the original (unstripped) push.
    return 1 if push_rc == 0 else push_rc


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="*", help="Notebook files to strip in-place")
    parser.add_argument(
        "--check-credentials", action="store_true",
        help="(pre-commit) Check staged notebooks for credential patterns only; do not strip outputs",
    )
    parser.add_argument(
        "--pre-push", metavar="REMOTE",
        help="(pre-push) Rewrite pushed commits to strip outputs. Pass remote name as argument; refs read from stdin.",
    )
    args = parser.parse_args()

    if args.check_credentials:
        return check_credentials_staged()

    if args.pre_push:
        stdin_text = sys.stdin.read()
        return pre_push_mode(args.pre_push, stdin_text)

    # Manual / legacy mode: strip files given on command line.
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
        print("Stripped notebook outputs:")
        for p in modified:
            print(f"  - {p}")
    if blocked:
        print("\nSensitive token pattern found in notebook source cell(s).")
        for p in blocked:
            print(f"  - {p}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
