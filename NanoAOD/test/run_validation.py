#!/usr/bin/env python3
"""Validation orchestrator for Bmm5/NanoAOD.

Produces NanoAOD output, validates it, compares against reference,
and saves a report.

Directory structure: <base_path>/<version>/<sample>/
  version = commit hash, tag, or "current"

Reference depends on test level:
  smoke/nominal - last commit (auto-produced via stash/build cycle)
  release       - tagged reference on EOS

The library mtime is recorded in test_info.json to detect whether
output needs regeneration.
"""

import argparse
import glob
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from functools import partial
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from test_config import (
    processes, REFERENCE_BASE, LOCAL_OUTPUT, DEFAULT_SAMPLE,
    PRE_RELEASE_SAMPLE, BMM5_DIR,
    get_git_version, get_git_commit, get_latest_tag,
)


class Report:
    """Collects report lines and writes to file and stdout."""

    def __init__(self):
        self.lines = []

    def print(self, text=""):
        print(text, flush=True)
        self.lines.append(text)

    def save(self, path):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            f.write("\n".join(self.lines) + "\n")


def run_script(script, args):
    """Run a Python script, capture output. Returns (returncode, stdout)."""
    cmd = [sys.executable, str(SCRIPT_DIR / script)] + args
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode, result.stdout


def get_lib_mtime():
    """Get the modification time of the Bmm5 plugin library."""
    cmssw_base = os.environ.get("CMSSW_BASE", "")
    scram_arch = os.environ.get("SCRAM_ARCH", "")
    lib = os.path.join(cmssw_base, "lib", scram_arch, "pluginBmm5NanoAODPlugins.so")
    if os.path.exists(lib):
        return os.path.getmtime(lib)
    return None


def output_is_valid(output_dir, nevents=None):
    """Check if existing output matches current library state and event count."""
    info_path = output_dir / "test_info.json"
    root_path = output_dir / "nanoaod_test.root"
    if not info_path.exists() or not root_path.exists():
        return False
    try:
        info = json.load(open(info_path))
        if info.get("lib_mtime") != get_lib_mtime():
            return False
        if nevents is not None and info.get("nevents") != nevents:
            return False
        return True
    except (json.JSONDecodeError, KeyError):
        return False


def get_processed_events(log_path):
    """Return number of input events processed by cmsRun, parsed from the log
    via extract_performance.py. Returns None on failure."""
    log_path = str(log_path)
    if not os.path.exists(log_path):
        return None
    rc, text = run_script("extract_performance.py", [log_path])
    if rc != 0:
        return None
    m = re.search(r"^\s*Events:\s+(\d+)", text, re.MULTILINE)
    return int(m.group(1)) if m else None


def produce_output(sample, level, output_dir, cmssw_base=None, nevents_args=None):
    """Run NanoAOD production. Returns True on success.
    If cmssw_base is given, run in that CMSSW environment."""
    run_args = ["-p", sample, "-l", level, "-o", str(output_dir)]
    if nevents_args:
        run_args += nevents_args
    if cmssw_base:
        # Run run_nanoaod.py inside the given CMSSW environment
        cmd = (
            f"cd {cmssw_base}/src && eval `scramv1 runtime -sh` && "
            f"python3 {SCRIPT_DIR}/run_nanoaod.py " + " ".join(run_args)
        )
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            output = (result.stdout + result.stderr).strip()
            if output:
                # Show last 20 lines of combined output
                lines = output.split("\n")
                for line in lines[-20:]:
                    print(f"  {line}", flush=True)
        return result.returncode == 0
    else:
        rc, _ = run_script("run_nanoaod.py", run_args)
        if rc != 0:
            return False
        # Record library mtime in test_info.json
        info_path = output_dir / "test_info.json"
        if info_path.exists():
            info = json.load(open(info_path))
            info["lib_mtime"] = get_lib_mtime()
            with open(info_path, "w") as f:
                json.dump(info, f, indent=2)
        return True


def build_reference_cmssw(ref_commit, report):
    """Build a fresh CMSSW area in /tmp with Bmm5 at the given commit.
    Returns the CMSSW_BASE path on success, None on failure."""
    cmssw_version = os.environ["CMSSW_VERSION"]
    work_dir = Path(f"/tmp/{os.environ.get('USER', 'nobody')}/bmm_ref_build")
    cmssw_dir = work_dir / cmssw_version

    # Check if already built at the right commit.
    # Also sanity-check a key source file and the plugin library, since
    # the marker has been observed to outlive a partially-deleted tree.
    marker = cmssw_dir / "src" / "Bmm5" / ".ref_commit"
    critical_src = cmssw_dir / "src" / "Bmm5" / "NanoAOD" / "python" / "nano_cff.py"
    plugin_lib = cmssw_dir / "lib" / os.environ.get("SCRAM_ARCH", "") / "pluginBmm5NanoAODPlugins.so"
    if (marker.exists() and marker.read_text().strip() == ref_commit
        and critical_src.exists() and plugin_lib.exists()):
        report.print(f"  Reusing existing build: {cmssw_dir}")
        return str(cmssw_dir)

    # Clean and create
    if cmssw_dir.exists():
        shutil.rmtree(str(cmssw_dir))
    work_dir.mkdir(parents=True, exist_ok=True)

    report.print(f"  Creating {cmssw_version} in {work_dir}...")
    result = subprocess.run(
        f"cd {work_dir} && scram p CMSSW {cmssw_version}",
        shell=True, capture_output=True, text=True
    )
    if result.returncode != 0:
        report.print(f"[FAIL] scram p failed: {result.stderr.strip()}")
        return None

    report.print(f"  Cloning Bmm5 at {ref_commit[:7]}...")
    src_dir = cmssw_dir / "src"
    # Clone from local repo to avoid SSH/network issues
    result = subprocess.run(
        f"cd {src_dir} && git clone {BMM5_DIR} Bmm5 && "
        f"cd Bmm5 && git checkout {ref_commit}",
        shell=True, capture_output=True, text=True
    )
    if result.returncode != 0:
        report.print(f"[FAIL] git clone/checkout failed: {result.stderr.strip()}")
        return None

    report.print(f"  Building...")
    result = subprocess.run(
        f"cd {src_dir} && eval `scramv1 runtime -sh` && scram b -j 8",
        shell=True, capture_output=True, text=True
    )
    if result.returncode != 0:
        report.print(f"[FAIL] scram build failed: {result.stderr[-500:]}")
        return None

    # Mark the commit
    marker.write_text(ref_commit)
    report.print(f"  Build complete: {cmssw_dir}")
    return str(cmssw_dir)


def ensure_reference(sample, level, ref_commit, ref_dir, report, nevents=None,
                     label="reference", force=False, check_nevents=False):
    """Ensure output exists for a given commit. Builds a separate CMSSW area if needed.
    Returns True if output is available.

    By default (check_nevents=False, used for the tag reference) the event-count
    match is NOT validated here; callers compare processed input events post-hoc
    via get_processed_events and re-invoke with force=True if they don't agree.

    For the working-tree/current-commit path, callers pass check_nevents=True so
    a request for different -n invalidates the cached output.

    If force=True, any existing outputs are removed first."""
    ref_root = ref_dir / "nanoaod_test.root"
    ref_info = ref_dir / "test_info.json"
    if not force and ref_root.exists() and ref_info.exists():
        if check_nevents and nevents is not None:
            try:
                info = json.load(open(ref_info))
                if info.get("nevents") != nevents:
                    report.print(f"{label.capitalize()} has {info.get('nevents')} events, "
                                 f"need {nevents}")
                else:
                    report.print(f"Using existing {label}: {ref_dir}")
                    return True
            except (json.JSONDecodeError, KeyError):
                pass
        else:
            report.print(f"Using existing {label}: {ref_dir}")
            return True

    if ref_dir.exists():
        for fname in ("nanoaod_test.root", "nanoaod_test.log", "test_info.json"):
            p = ref_dir / fname
            if p.exists():
                p.unlink()

    report.print(f"Producing {label} for commit {ref_commit[:7]}...")

    # Build a fresh CMSSW area with Bmm5 at the given commit
    cmssw_base = build_reference_cmssw(ref_commit, report)
    if not cmssw_base:
        return False

    # Produce output using the built CMSSW environment
    report.print(f"  Running NanoAOD production for {label}...")
    ref_dir.mkdir(parents=True, exist_ok=True)
    nevents_args = ["-n", str(nevents)] if nevents else []
    success = produce_output(sample, level, ref_dir, cmssw_base=cmssw_base, nevents_args=nevents_args)

    if success:
        report.print(f"  Output produced: {ref_root}")
    else:
        report.print(f"[FAIL] {label.capitalize()} production failed")

    return success



def run_sample(sample, ref_dir, ref_label, ref_commit, level, report, output_dir, save_dir=None, nevents=None, current_commit=None):
    """Run full validation for one sample. Returns (passed, failed, skipped)."""
    passed = 0
    failed = 0
    skipped = 0

    output_root = output_dir / "nanoaod_test.root"
    output_log = output_dir / "nanoaod_test.log"
    ref_root = ref_dir / "nanoaod_test.root"
    ref_log = ref_dir / "nanoaod_test.log"

    nevents_args = ["-n", str(nevents)] if nevents else []

    # Step 1: Produce current output. Always run if no cached output for this
    # library/commit; we need the event count before deciding what to do about
    # the reference.
    current_label = f"commit {current_commit[:7]}" if current_commit else "working tree"
    report.print(f"\n--- Production ({current_label}) ---")
    if current_commit:
        if not ensure_reference(sample, level, current_commit, output_dir, report, nevents,
                                label=current_label, check_nevents=True):
            report.print(f"[FAIL] Production failed for {current_label}")
            return 0, 1, 0
    else:
        if output_is_valid(output_dir, nevents):
            report.print(f"Using cached output (library unchanged)")
        else:
            if not produce_output(sample, level, output_dir, nevents_args=nevents_args):
                report.print(f"[FAIL] NanoAOD production failed")
                return 0, 1, 0
    report.print(f"[PASS] Output produced")
    passed += 1

    # Step 2: Reference. Use cached output if present; otherwise build with the
    # same -n as current. Event-count agreement is verified post-hoc below.
    report.print(f"\n--- Reference ({ref_label}) ---")
    if not ensure_reference(sample, level, ref_commit, ref_dir, report, nevents):
        return 0, 1, 0

    # Step 2b: Verify processed input-event counts agree. If not, rebuild the
    # reference with the count the current run achieved.
    current_events = get_processed_events(output_log)
    ref_events = get_processed_events(ref_log)
    if current_events is not None and ref_events is not None:
        if current_events != ref_events:
            report.print(f"Processed events mismatch: current={current_events}, reference={ref_events}")
            report.print(f"Rebuilding reference at -n {current_events}...")
            if not ensure_reference(sample, level, ref_commit, ref_dir, report,
                                    current_events, force=True):
                return 0, 1, 0
            ref_events = get_processed_events(ref_log)
            if ref_events != current_events:
                report.print(f"[FAIL] Reference still mismatches after rebuild "
                             f"(current={current_events}, reference={ref_events})")
                return 0, 1, 0
            report.print(f"Reference rebuilt at {ref_events} events")
        else:
            report.print(f"Processed events: {current_events} (matches reference)")
    else:
        report.print(f"[WARN] Could not parse processed-event counts "
                     f"(current={current_events}, reference={ref_events})")

    # Step 3: Output validation
    report.print(f"\n--- Validation ---")
    rc, text = run_script("check_output.py", [str(output_root)])
    if rc == 0:
        report.print("[PASS] Output validation passed")
        passed += 1
    else:
        for line in text.strip().split("\n"):
            report.print(f"  {line}")
        report.print("[FAIL] Output validation failed")
        failed += 1

    # Step 4: Content comparison (smoke only - single-threaded for determinism)
    if level == "smoke" and ref_root.exists():
        report.print(f"\n--- Content Comparison ---")
        exclude = processes[sample].get("exclude_branches", "")
        cmp_args = [str(ref_root), str(output_root),
                    "--labels", "reference", "current"]
        if exclude:
            cmp_args += ["--pattern", f"^(?!{exclude})"]
        rc, text = run_script("compare_root_files.py", cmp_args)
        if text.strip():
            for line in text.strip().split("\n"):
                report.print(f"  {line}")
        if rc == 0:
            report.print("[PASS] Content identical to reference")
            passed += 1
        else:
            report.print("[FAIL] Content differs from reference")
            failed += 1

    # Step 5: Selection yields (nominal/release)
    selections = processes[sample].get("selections", {})
    if selections:
        report.print(f"\n--- Selection Yields ---")
        sel_json = json.dumps(selections)
        yield_args = [str(output_root), "--selections", sel_json]
        if ref_root.exists():
            yield_args += ["--ref", str(ref_root)]
        rc, text = run_script("check_yields.py", yield_args)
        if text.strip():
            for line in text.strip().split("\n"):
                report.print(f"  {line}")

    # Step 6: Comparison plots
    plots = processes[sample].get("plots", {})
    if plots and ref_root.exists():
        report.print(f"\n--- Plots ---")
        plot_dir = output_dir / "plots"
        plots_json = json.dumps(plots)
        rc, text = run_script("make_plots.py",
            [str(output_root), str(ref_root), "--plots", plots_json,
             "-o", str(plot_dir)])
        if text.strip():
            for line in text.strip().split("\n"):
                report.print(f"  {line}")
        report.print(f"  Saved to: {plot_dir}")

    # Step 7: Performance
    perf_args = [str(output_log), "--rootfile", str(output_root)]
    report.print(f"\n--- Performance (current) ---")
    rc, text = run_script("extract_performance.py", perf_args)
    if text.strip():
        for line in text.strip().split("\n"):
            report.print(f"  {line}")
    ref_log = ref_dir / "nanoaod_test.log"
    if ref_log.exists():
        ref_perf_args = [str(ref_log), "--rootfile", str(ref_root)]
        report.print(f"\n--- Performance (reference) ---")
        rc, text = run_script("extract_performance.py", ref_perf_args)
        if text.strip():
            for line in text.strip().split("\n"):
                report.print(f"  {line}")

    # Step 8: Save as tagged reference on EOS
    if save_dir:
        report.print(f"\n--- Saving reference ---")
        save_dir.mkdir(parents=True, exist_ok=True)
        for fname in ["nanoaod_test.root", "nanoaod_test.log", "test_info.json"]:
            src = output_dir / fname
            if src.exists():
                shutil.copy2(str(src), str(save_dir / fname))
        report.print(f"Saved to: {save_dir}")

    return passed, failed, skipped


def main():
    parser = argparse.ArgumentParser(description="Bmm5/NanoAOD validation")
    parser.add_argument("-p", "--sample",
                        help="sample name (default: DEFAULT_SAMPLE for smoke/nominal, all for release)")
    parser.add_argument("-r", "--ref",
                        help="reference tag or commit (default: HEAD for smoke/nominal, latest tag for pre-release/release)")
    parser.add_argument("-l", "--level", default="nominal",
                        choices=["smoke", "nominal", "pre-release", "release"])
    parser.add_argument("-o", "--output", help="base output directory")
    parser.add_argument("-n", "--nevents",
                        type=lambda s: -1 if s.lower() == "all" else int(s),
                        help="override number of events (invalidates cached output); 'all' or -1 = all events in input")
    parser.add_argument("--commit",
                        help="build and test a specific commit (instead of current working tree)")
    parser.add_argument("-s", "--save", action="store_true",
                        help="save output as tagged reference on EOS")
    args = parser.parse_args()

    if not os.environ.get("CMSSW_BASE"):
        print("Error: CMSSW environment not set. Run: eval `scramv1 runtime -sh`",
              file=sys.stderr)
        sys.exit(1)

    # Resolve --commit to full hash
    current_commit = None
    if args.commit:
        result = subprocess.run(
            ["git", "rev-parse", args.commit],
            capture_output=True, text=True, cwd=BMM5_DIR
        )
        if result.returncode != 0:
            print(f"Error: cannot resolve commit '{args.commit}'", file=sys.stderr)
            sys.exit(1)
        current_commit = result.stdout.strip()

    # Determine reference tag
    ref_tag = None
    if args.level in ("pre-release", "release") or args.ref:
        ref_tag = args.ref or get_latest_tag()
        if not ref_tag:
            print("Error: no Bmm5 tags found and no -r specified", file=sys.stderr)
            sys.exit(1)

    # Determine samples
    if args.sample:
        samples = [args.sample]
    elif args.level == "release":
        samples = [name for name in sorted(processes) if processes[name].get("input")]
    elif args.level == "pre-release":
        samples = [PRE_RELEASE_SAMPLE]
    else:
        samples = [DEFAULT_SAMPLE]

    # Base output path
    base_output = Path(args.output or LOCAL_OUTPUT)

    # Set up report
    report = Report()
    report.print(f"Bmm5/NanoAOD Validation Report")
    report.print(f"Date:       {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.print(f"Test level: {args.level}")
    if ref_tag:
        report.print(f"Reference:  {ref_tag}")
    if current_commit:
        report.print(f"Target:     commit {current_commit[:7]}")
    report.print(f"Samples:    {', '.join(samples)}")

    total_pass = 0
    total_fail = 0
    total_skip = 0

    for sample in samples:
        report.print(f"\n{'=' * 64}")
        report.print(f"Sample: {sample}")
        report.print(f"{'=' * 64}")

        # Determine directories and reference commit
        if ref_tag:
            # Tag reference: check local first, fall back to EOS
            local_ref = base_output / args.level / ref_tag / sample
            eos_ref = Path(REFERENCE_BASE) / args.level / ref_tag / sample
            if (local_ref / "nanoaod_test.root").exists():
                ref_dir = local_ref
            else:
                try:
                    eos_exists = (eos_ref / "nanoaod_test.root").exists()
                except PermissionError:
                    eos_exists = False
                if eos_exists:
                    ref_dir = eos_ref
                else:
                    # Will be built by ensure_reference
                    ref_dir = local_ref
            ref_label = ref_tag
            # Resolve tag to commit hash
            ref_commit = subprocess.run(
                ["git", "rev-parse", ref_tag],
                capture_output=True, text=True, cwd=BMM5_DIR
            ).stdout.strip()
        else:
            # Default: compare against HEAD (pre-commit workflow)
            commit = get_git_commit()
            ref_dir = base_output / args.level / commit[:7] / sample
            ref_label = f"commit {commit[:7]}"
            ref_commit = commit

        # Prevent comparing a commit against itself
        if current_commit and current_commit == ref_commit:
            report.print(f"[SKIP] Target and reference are the same commit ({ref_commit[:7]})")
            total_skip += 1
            continue

        if current_commit:
            output_dir = base_output / args.level / current_commit[:7] / sample
        else:
            output_dir = base_output / args.level / "current" / sample

        save_dir = None
        if args.save:
            git_version = get_git_version()
            save_dir = Path(REFERENCE_BASE) / "release" / git_version / sample

        p, f, s = run_sample(sample, ref_dir, ref_label, ref_commit, args.level,
                             report, output_dir, save_dir, args.nevents,
                             current_commit=current_commit)

        if f > 0:
            report.print(f"\n>>> {sample}: FAILED ({p} passed, {f} failed)")
            total_fail += 1
        else:
            report.print(f"\n>>> {sample}: PASSED ({p} passed)")
            total_pass += 1
        total_skip += s

    report.print(f"\n{'=' * 64}")
    report.print(f"Overall: {total_pass} passed, {total_fail} failed, {total_skip} skipped")
    report.print(f"{'=' * 64}")

    # Save report
    report_path = base_output / "validation_report.txt"
    report.save(report_path)
    print(f"\nReport saved: {report_path}", flush=True)

    sys.exit(1 if total_fail > 0 else 0)


if __name__ == "__main__":
    main()
