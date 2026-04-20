#!/usr/bin/env python3
"""Check event yields for selections defined in test_config.py.

Uses PyROOT's TTree::GetEntries(selection) to count events passing
each cut. Optionally compares against a reference file.
"""

import argparse
import json
import sys

import ROOT
ROOT.gROOT.SetBatch(True)


def get_yields(root_file, selections):
    """Return {name: count} for each selection."""
    f = ROOT.TFile.Open(root_file)
    if not f or f.IsZombie():
        print(f"Error: cannot open {root_file}", file=sys.stderr)
        sys.exit(2)
    tree = f.Get("Events")
    if not tree:
        print(f"Error: no Events tree in {root_file}", file=sys.stderr)
        sys.exit(2)
    yields = {}
    for name, cut in selections.items():
        yields[name] = tree.GetEntries(cut)
    f.Close()
    return yields


def main():
    parser = argparse.ArgumentParser(description="Check selection yields")
    parser.add_argument("rootfile", help="NanoAOD ROOT file")
    parser.add_argument("--selections", required=True,
                        help="JSON string of {name: cut} selections")
    parser.add_argument("--ref", help="reference ROOT file for comparison")
    args = parser.parse_args()

    selections = json.loads(args.selections)
    if not selections:
        return

    current = get_yields(args.rootfile, selections)
    reference = get_yields(args.ref, selections) if args.ref else None

    max_name = max(len(n) for n in selections)
    for name in selections:
        c = current[name]
        if reference:
            r = reference[name]
            diff = c - r
            sign = "+" if diff > 0 else ""
            print(f"  {name:{max_name}s}  current: {c:6d}   reference: {r:6d}   (diff: {sign}{diff})")
        else:
            print(f"  {name:{max_name}s}  {c}")


if __name__ == "__main__":
    main()
