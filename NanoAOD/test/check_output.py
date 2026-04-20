#!/usr/bin/env python3
"""Validate a NanoAOD output file for basic correctness."""

import sys
import os
import argparse
import numpy as np
import uproot


REQUIRED_TREES = ["Events", "Runs", "LuminosityBlocks"]

REQUIRED_BRANCHES = [
    "mm_kin_mass",
    "mm_kin_vtx_prob",
    "mm_kin_lxy",
    "mm_kin_slxy",
    "mm_kin_alpha",
    "mm_kin_l3d",
    "mm_kin_sl3d",
]


def flatten_jagged(arr):
    """Flatten a potentially jagged numpy object array."""
    if arr.dtype.kind == "O":
        parts = [np.asarray(x).ravel() for x in arr if len(x) > 0]
        return np.concatenate(parts) if parts else np.array([], dtype=np.float64)
    return arr.ravel()


def main():
    parser = argparse.ArgumentParser(description="Validate NanoAOD output file")
    parser.add_argument("rootfile", help="NanoAOD ROOT file to validate")
    args = parser.parse_args()

    passed = 0
    failed = 0

    print(f"Checking: {args.rootfile}")

    # Check 1: file exists and is non-empty
    if not os.path.exists(args.rootfile):
        print("  [FAIL] File does not exist")
        sys.exit(1)
    if os.path.getsize(args.rootfile) == 0:
        print("  [FAIL] File is empty")
        sys.exit(1)

    # Check 2: valid ROOT file
    try:
        f = uproot.open(args.rootfile)
    except Exception as e:
        print(f"  [FAIL] Not a valid ROOT file: {e}")
        sys.exit(1)
    print("  [PASS] File structure OK")
    passed += 1

    # Check 3: expected trees
    trees = [k.split(";")[0] for k in f.keys(filter_classname="TTree")]
    missing_trees = [t for t in REQUIRED_TREES if t not in trees]
    if missing_trees:
        print(f"  [FAIL] Missing trees: {', '.join(missing_trees)}")
        failed += 1
    else:
        print(f"  [PASS] Trees present: {', '.join(REQUIRED_TREES)}")
        passed += 1

    if "Events" not in trees:
        print("  [FAIL] Cannot continue without Events tree")
        sys.exit(1)

    events = f["Events"]
    nevents = events.num_entries
    branches = set(events.keys())

    # Check 4: key Bmm branches
    missing_branches = [b for b in REQUIRED_BRANCHES if b not in branches]
    if missing_branches:
        print(f"  [FAIL] Missing branches: {', '.join(missing_branches)}")
        failed += 1
    else:
        print(f"  [PASS] Key branches present ({len(REQUIRED_BRANCHES)}/{len(REQUIRED_BRANCHES)})")
        passed += 1

    # Check 5: sanity checks on physics values (only for valid kinematic fits)
    issues = []

    # Build mask for valid kinematic fits
    kin_valid_mask = None
    if "mm_kin_valid" in branches:
        kin_valid_arr = events["mm_kin_valid"].array(library="np")
        kin_valid_mask = np.array([np.asarray(x) == 1 for x in kin_valid_arr], dtype=object)

    def flatten_valid(arr, mask):
        """Flatten jagged array, keeping only entries where mask is True."""
        if mask is None:
            return flatten_jagged(arr)
        parts = []
        for vals, m in zip(arr, mask):
            v = np.asarray(vals)
            mm = np.asarray(m)
            if len(v) > 0 and len(mm) > 0:
                parts.append(v[mm])
        return np.concatenate(parts) if parts else np.array([], dtype=np.float64)

    if "mm_kin_mass" in branches:
        mass = flatten_valid(events["mm_kin_mass"].array(library="np"), kin_valid_mask)
        if len(mass) > 0:
            n_nan = int(np.sum(np.isnan(mass)))
            if n_nan > 0:
                issues.append(f"mm_kin_mass has {n_nan} NaN values (valid fits only)")

    if "mm_kin_vtx_prob" in branches:
        vtx_prob = flatten_valid(events["mm_kin_vtx_prob"].array(library="np"), kin_valid_mask)
        if len(vtx_prob) > 0:
            valid = vtx_prob[~np.isnan(vtx_prob)]
            if len(valid) > 0 and (np.any(valid < 0) or np.any(valid > 1)):
                issues.append("mm_kin_vtx_prob has values outside [0, 1] (valid fits only)")

    if "mm_kin_lxy" in branches:
        lxy = flatten_valid(events["mm_kin_lxy"].array(library="np"), kin_valid_mask)
        if len(lxy) > 0:
            valid = lxy[~np.isnan(lxy)]
            if len(valid) > 0 and np.any(valid < 0):
                issues.append("mm_kin_lxy has negative values (valid fits only)")

    if issues:
        for issue in issues:
            print(f"  [FAIL] {issue}")
        failed += 1
    else:
        print("  [PASS] Sanity checks passed")
        passed += 1

    print(f"  Events: {nevents}")

    if failed > 0:
        print(f"\n  Result: {failed} check(s) FAILED")
        sys.exit(1)
    else:
        print(f"\n  Result: all {passed} checks passed")
        sys.exit(0)


if __name__ == "__main__":
    main()
