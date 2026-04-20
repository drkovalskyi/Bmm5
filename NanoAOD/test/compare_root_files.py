#!/usr/bin/env python3
"""Compare two ROOT files branch-by-branch.

Auto-discovers all TTrees. Default mode: entry-by-entry exact
comparison using uproot + numpy with NaN-aware matching.
"""

import sys
import os
import re
import json
import argparse
import numpy as np
import uproot


def discover_trees(f):
    """Return sorted list of TTree names in a ROOT file."""
    return sorted(k.split(";")[0] for k in f.keys(filter_classname="TTree"))


def compute_stats(arr):
    """Compute min/max/mean/std/n for a (possibly jagged) array."""
    if arr.dtype.kind == "O":
        parts = [np.asarray(x, dtype=float).ravel() for x in arr if len(x) > 0]
        flat = np.concatenate(parts) if parts else np.array([], dtype=float)
    else:
        flat = arr.ravel().astype(float)
    if len(flat) == 0:
        return {"n": 0, "min": None, "max": None, "mean": None, "std": None}
    valid = flat[~np.isnan(flat)]
    if len(valid) == 0:
        return {"n": len(flat), "min": float("nan"), "max": float("nan"),
                "mean": float("nan"), "std": float("nan")}
    return {
        "n": len(flat),
        "min": float(np.min(valid)),
        "max": float(np.max(valid)),
        "mean": float(np.mean(valid)),
        "std": float(np.std(valid)),
    }


def arrays_match(a, b, tolerance):
    """Compare two arrays (scalar or array) with optional tolerance."""
    a = np.asarray(a)
    b = np.asarray(b)
    if a.shape != b.shape:
        return False
    if a.dtype.kind in ("U", "S", "O") or b.dtype.kind in ("U", "S", "O"):
        return np.array_equal(a, b)
    if tolerance is not None:
        return np.allclose(a, b, rtol=0, atol=tolerance, equal_nan=True)
    # Exact match with NaN-awareness for floats
    if a.dtype.kind == "f" or b.dtype.kind == "f":
        return np.allclose(a, b, rtol=0, atol=0, equal_nan=True)
    return np.array_equal(a, b)


def compare_branch_deep(arr1, arr2, tolerance, labels=("file1", "file2")):
    """Entry-by-entry comparison. Returns (match, first_mismatch_info)."""
    l1, l2 = labels
    if len(arr1) != len(arr2):
        return False, f"different entry counts ({l1}={len(arr1)}, {l2}={len(arr2)})"

    if arr1.dtype.kind == "O" or arr2.dtype.kind == "O":
        # Jagged arrays: compare entry by entry
        for i in range(len(arr1)):
            if not arrays_match(arr1[i], arr2[i], tolerance):
                return False, f"entry {i}:\n      {l1}: {arr1[i]}\n      {l2}: {arr2[i]}"
        return True, None

    # Flat arrays: vectorized comparison
    if arrays_match(arr1, arr2, tolerance):
        return True, None

    # Find first mismatch for reporting
    if tolerance is not None:
        diff = np.where(~np.isclose(arr1, arr2, rtol=0, atol=tolerance, equal_nan=True))[0]
    elif arr1.dtype.kind == "f":
        diff = np.where(~np.isclose(arr1, arr2, rtol=0, atol=0, equal_nan=True))[0]
    else:
        diff = np.where(arr1 != arr2)[0]
    if len(diff) > 0:
        i = diff[0]
        return False, f"entry {i}:\n      {l1}: {arr1[i]}\n      {l2}: {arr2[i]}"
    return True, None


def compare_fast(tree1, tree2, branch_map, pattern, ignore):
    """Fast comparison: entry counts and compressed sizes only."""
    branches1 = set(tree1.keys())
    branches2 = set(tree2.keys())
    results = {"ok": 0, "different": 0, "added": [], "removed": []}

    for b in sorted(branches1):
        mapped = branch_map.get(b, b)
        if pattern and not re.search(pattern, b):
            continue
        if mapped not in branches2:
            if not ignore:
                results["removed"].append(b)
            continue
        n1 = tree1[b].num_entries
        n2 = tree2[mapped].num_entries
        s1 = tree1[b].compressed_bytes
        s2 = tree2[mapped].compressed_bytes
        if n1 == n2 and s1 == s2:
            results["ok"] += 1
        else:
            print(f"  [DIFFERENT] {b}")
            if n1 != n2:
                print(f"    entries: {n1} vs {n2}")
            if s1 != s2:
                print(f"    compressed bytes: {s1} vs {s2}")
            results["different"] += 1

    for b in sorted(branches2):
        reverse = {v: k for k, v in branch_map.items()}
        orig = reverse.get(b, b)
        if pattern and not re.search(pattern, b):
            continue
        if orig not in branches1 and b not in branches1:
            if not ignore:
                results["added"].append(b)

    return results


def compare_deep(tree1, tree2, branch_map, pattern, ignore, tolerance, verbose, labels=("file1", "file2")):
    """Deep entry-by-entry comparison."""
    branches1 = set(tree1.keys())
    branches2 = set(tree2.keys())
    results = {"ok": 0, "different": 0, "added": [], "removed": []}

    for b in sorted(branches1):
        mapped = branch_map.get(b, b)
        if pattern and not re.search(pattern, b):
            continue
        if mapped not in branches2:
            if not ignore:
                results["removed"].append(b)
            continue

        try:
            arr1 = tree1[b].array(library="np")
        except Exception:
            arr1 = None
        try:
            arr2 = tree2[mapped].array(library="np")
        except Exception:
            arr2 = None
        if arr1 is None and arr2 is None:
            results["ok"] += 1  # both unreadable - same
            continue
        if arr1 is None or arr2 is None:
            print(f"  [DIFFERENT] {b} (readable in only one file)")
            results["different"] += 1
            continue

        match, info = compare_branch_deep(arr1, arr2, tolerance, labels)
        if match:
            results["ok"] += 1
        else:
            print(f"  [DIFFERENT] {b}")
            if info:
                print(f"    {info}")
            if verbose:
                s1 = compute_stats(arr1)
                s2 = compute_stats(arr2)
                l1, l2 = labels
                print(f"    {l1}: n={s1['n']} min={s1['min']} max={s1['max']} "
                      f"mean={s1['mean']:.6g} std={s1['std']:.6g}"
                      if s1['mean'] is not None else f"    {l1}: n={s1['n']} (empty)")
                print(f"    {l2}: n={s2['n']} min={s2['min']} max={s2['max']} "
                      f"mean={s2['mean']:.6g} std={s2['std']:.6g}"
                      if s2['mean'] is not None else f"    {l2}: n={s2['n']} (empty)")
            results["different"] += 1

    reverse = {v: k for k, v in branch_map.items()}
    for b in sorted(branches2):
        orig = reverse.get(b, b)
        if pattern and not re.search(pattern, b):
            continue
        if orig not in branches1 and b not in branches1:
            if not ignore:
                results["added"].append(b)

    return results


def main():
    parser = argparse.ArgumentParser(description="Compare ROOT files branch-by-branch")
    parser.add_argument("files", nargs=2, help="ROOT files to compare")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="show statistics for differing branches")
    parser.add_argument("-f", "--fast", action="store_true",
                        help="compare only entry counts and compressed sizes")
    parser.add_argument("-p", "--pattern", type=str,
                        help="only compare branches matching regex pattern")
    parser.add_argument("-t", "--tree", type=str,
                        help="compare specific tree only")
    parser.add_argument("-m", "--map", type=str, metavar="FILE",
                        help="JSON file mapping renamed branches (old->new)")
    parser.add_argument("-i", "--ignore", action="store_true",
                        help="suppress missing/added branch warnings")
    parser.add_argument("--tolerance", type=float, default=None,
                        help="absolute tolerance for float comparison")
    parser.add_argument("--labels", nargs=2, default=["file1", "file2"],
                        metavar=("LABEL1", "LABEL2"),
                        help="labels for the two files in diff output")

    args = parser.parse_args()

    branch_map = {}
    if args.map:
        with open(args.map) as f:
            branch_map = json.load(f)

    try:
        f1 = uproot.open(args.files[0])
        f2 = uproot.open(args.files[1])
    except Exception as e:
        print(f"Error opening files: {e}", file=sys.stderr)
        sys.exit(2)

    trees1 = discover_trees(f1)
    trees2 = discover_trees(f2)
    all_trees = sorted(set(trees1) | set(trees2))

    any_diff = False

    for tree_name in all_trees:
        if args.tree and tree_name != args.tree:
            continue

        if tree_name not in trees1:
            print(f"Tree {tree_name}: only in file2")
            any_diff = True
            continue
        if tree_name not in trees2:
            print(f"Tree {tree_name}: only in file1")
            any_diff = True
            continue

        print(f"Comparing tree: {tree_name}")
        tree1 = f1[tree_name]
        tree2 = f2[tree_name]

        if args.fast:
            results = compare_fast(tree1, tree2, branch_map, args.pattern, args.ignore)
        else:
            results = compare_deep(tree1, tree2, branch_map, args.pattern,
                                   args.ignore, args.tolerance, args.verbose,
                                   tuple(args.labels))

        for b in results["removed"]:
            print(f"  [REMOVED]   {b}")
        for b in results["added"]:
            print(f"  [ADDED]     {b}")

        total_diff = results["different"] + len(results["added"]) + len(results["removed"])
        print(f"  Matching branches: {results['ok']}")
        print(f"  Different branches: {results['different']}")
        if not args.ignore:
            print(f"  Added branches: {len(results['added'])}")
            print(f"  Removed branches: {len(results['removed'])}")

        if total_diff > 0:
            any_diff = True

    sys.exit(1 if any_diff else 0)


if __name__ == "__main__":
    main()
