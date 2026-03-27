#!/usr/bin/env python3

import json
import argparse

def count_total_lumis(json_path):
    with open(json_path, "r") as f:
        data = json.load(f)

    total_lumis = 0

    # data is expected to be { "run": [ [start, end], [start, end], ... ], ... }
    for run, lumi_ranges in data.items():
        for lumi_range in lumi_ranges:
            if len(lumi_range) != 2:
                raise ValueError(f"Invalid lumi range in run {run}: {lumi_range}")
            start, end = lumi_range
            # assuming ranges are inclusive
            total_lumis += end - start + 1

    return total_lumis


def main():
    parser = argparse.ArgumentParser(
        description="Compute total number of lumi sections from a CMS style JSON file"
    )
    parser.add_argument("json_file", help="Path to the JSON file")
    args = parser.parse_args()

    total = count_total_lumis(args.json_file)
    print(f"Total number of lumi sections: {total}")


if __name__ == "__main__":
    main()
