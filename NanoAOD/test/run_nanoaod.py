#!/usr/bin/env python3
"""Run a single NanoAOD production job for testing.

Reads process definitions from test_config.py.
Output: /tmp/$USER/bmm_test/<sample>/
Writes test_info.json alongside output for reference verification.
"""

import argparse
import json
import os
import shutil
import socket
import subprocess
import sys
from datetime import datetime
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from test_config import (
    processes, get_nevents, get_git_version, get_git_commit, LOCAL_OUTPUT,
)


def list_samples():
    for name in sorted(processes):
        p = processes[name]
        status = "TBD" if not p.get("input") else "ready"
        print(f"  {name:25s} {p.get('type','?'):5s} {p.get('era','?'):30s} [{status}]")


def build_cmsdriver_cmd(proc, nevents, output_root, config_py, level="nominal"):
    cmd = ["cmsDriver.py", "RECO"]
    cmd += ["--conditions", proc["conditions"]]
    cmd += ["--datatier", "NANOAOD"]
    cmd += ["--era", proc["era"]]
    cmd += ["--eventcontent", "NANOAODSIM" if proc["type"] == "mc" else "NANOAOD"]
    cmd += ["--filein", proc["input"][0]]
    cmd += ["--fileout", f"file:{output_root}"]
    # Single-threaded for smoke to preserve event order for
    # entry-by-entry content comparison
    cmd += ["--nThreads", "1" if level == "smoke" else "16"]
    cmd += ["-n", str(nevents)]
    cmd += ["--no_exec"]
    cmd += ["--python_filename", str(config_py)]
    cmd += ["--scenario", "pp"]
    cmd += ["--step", "NANO"]
    if proc["type"] == "mc":
        cmd += ["--mc"]
    for c in proc.get("customise", []):
        cmd += [f"--customise={c}"]
    for cc in proc.get("customise_commands", []):
        cmd += [f"--customise_commands={cc}"]
    return cmd


def main():
    parser = argparse.ArgumentParser(description="Run NanoAOD production for testing")
    parser.add_argument("-p", "--sample", help="sample name from test_config.py")
    parser.add_argument("-l", "--level", default="nominal", choices=["smoke", "nominal", "pre-release", "release"])
    parser.add_argument("-o", "--output", help="output directory (overrides default)")
    parser.add_argument("-n", "--nevents",
                        type=lambda s: -1 if s.lower() == "all" else int(s),
                        help="number of events (overrides level); 'all' or -1 = all events in input")
    parser.add_argument("--list", action="store_true", help="list available samples and exit")
    args = parser.parse_args()

    if args.list or not args.sample:
        list_samples()
        return

    if args.sample not in processes:
        print(f"Error: unknown sample '{args.sample}'", file=sys.stderr)
        list_samples()
        sys.exit(1)

    proc = processes[args.sample]
    if not proc.get("input"):
        print(f"Error: sample '{args.sample}' has no input files configured", file=sys.stderr)
        sys.exit(1)

    if not os.environ.get("CMSSW_BASE"):
        print("Error: CMSSW environment not set. Run: eval `scramv1 runtime -sh`", file=sys.stderr)
        sys.exit(1)

    nevents = args.nevents or get_nevents(proc["nevents"], args.level)
    output_dir = Path(args.output) if args.output else Path(LOCAL_OUTPUT) / args.sample
    output_dir.mkdir(parents=True, exist_ok=True)

    output_root = output_dir / "nanoaod_test.root"
    output_log = output_dir / "nanoaod_test.log"
    config_py = output_dir / "nanoaod_test_cfg.py"
    test_info = output_dir / "test_info.json"

    # Remove stale outputs up front so an early cmsRun failure can't be
    # masked by a leftover ROOT file from a previous successful run.
    for stale in (output_root, output_log, config_py, test_info):
        if stale.exists():
            stale.unlink()

    git_version = get_git_version()
    git_commit = get_git_commit()

    print(f"=== Sample: {args.sample} ===")
    print(f"Version: {git_version}")
    print(f"Level:   {args.level}")
    print(f"Events:  {nevents}")
    print(f"Output:  {output_dir}")

    # Generate cmsDriver config
    print("Generating cmsDriver config...")
    cmd = build_cmsdriver_cmd(proc, nevents, output_root, config_py, args.level)
    print(" ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)
    print(result.stdout)

    # Append extra input files if any
    if len(proc["input"]) > 1:
        with open(config_py, "a") as f:
            f.write("\n# Additional input files\n")
            for inp in proc["input"][1:]:
                f.write(f"process.source.fileNames.append('{inp}')\n")

    # Run cmsRun
    print("Running cmsRun...")
    with open(output_log, "w") as log:
        ret = subprocess.run(["cmsRun", str(config_py)], stdout=log, stderr=subprocess.STDOUT,
                             cwd=str(output_dir))
    if ret.returncode != 0:
        # Check if output is valid despite non-zero exit (shutdown crashes)
        valid = output_root.exists() and output_root.stat().st_size > 0
        if valid:
            try:
                import uproot
                f = uproot.open(str(output_root))
                n = f["Events"].num_entries
                if n >= nevents:
                    print(f"WARNING: cmsRun exit code {ret.returncode} but output valid ({n} events)")
                    valid = True
                else:
                    valid = False
            except Exception:
                valid = False

        if not valid:
            print(f"ERROR: cmsRun failed with exit code {ret.returncode}")
            print(f"See log: {output_log}")
            with open(output_log) as f:
                lines = f.readlines()
                for line in lines[-20:]:
                    print(line, end="")
            sys.exit(ret.returncode)

    print("cmsRun completed successfully")

    # Write test_info.json
    info = {
        "sample": args.sample,
        "type": proc["type"],
        "era": proc["era"],
        "conditions": proc["conditions"],
        "input": proc["input"],
        "nevents": nevents,
        "level": args.level,
        "git_version": git_version,
        "git_commit": git_commit,
        "timestamp": datetime.now().isoformat(),
        "hostname": socket.gethostname(),
    }
    with open(output_dir / "test_info.json", "w") as f:
        json.dump(info, f, indent=2)

    print(f"Output: {output_root}")
    print(f"Log:    {output_log}")
    print(f"Info:   {output_dir / 'test_info.json'}")


if __name__ == "__main__":
    main()
