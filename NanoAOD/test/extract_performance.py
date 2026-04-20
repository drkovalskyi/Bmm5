#!/usr/bin/env python3
"""Extract performance metrics from a cmsRun log file.

Extracts: event count, time/event (excluding first event), peak RSS.
Uses "Begin processing" timestamps to compute per-event timing
without first-event initialization overhead.
"""

import sys
import os
import re
import json
import argparse
from datetime import datetime


def parse_event_timestamps(logfile_path):
    """Extract per-event timestamps from 'Begin processing' lines."""
    timestamps = []
    pattern = re.compile(
        r"Begin processing the (\d+)\w+ record\..*at (\d{2}-\w{3}-\d{4} \d{2}:\d{2}:\d{2}\.\d+)"
    )
    with open(logfile_path) as f:
        for line in f:
            m = pattern.search(line)
            if m:
                event_num = int(m.group(1))
                time_str = m.group(2)
                try:
                    ts = datetime.strptime(time_str, "%d-%b-%Y %H:%M:%S.%f")
                    timestamps.append((event_num, ts))
                except ValueError:
                    pass
    return timestamps


def extract_metrics(logfile_path, skip_events=1):
    """Parse a cmsRun log and return performance metrics dict.
    skip_events: number of initial events to exclude from time/event."""
    max_rss = None
    nevents = None
    total_time = None
    loop_time = None

    with open(logfile_path) as f:
        for line in f:
            m = re.search(r"MemoryCheck.*?RSS\s+(\S+)", line)
            if m:
                try:
                    rss = float(m.group(1))
                except ValueError:
                    # Interleaved log output can corrupt the RSS field
                    continue
                if max_rss is None or rss > max_rss:
                    max_rss = rss
                continue

            m = re.search(r"TrigReport Events total = (\d+)", line)
            if m:
                nevents = int(m.group(1))
                continue

            m = re.search(r"TimeReport> Time report complete in\s+(\S+)", line)
            if m:
                try:
                    total_time = float(m.group(1))
                except ValueError:
                    pass
                continue

            m = re.search(r"Total loop:\s+(\S+)", line)
            if m:
                try:
                    loop_time = float(m.group(1))
                except ValueError:
                    pass
                continue

    # Compute time/event from timestamps, excluding first N events
    timestamps = parse_event_timestamps(logfile_path)
    time_per_event = None
    events_measured = None

    if len(timestamps) > skip_events + 1:
        # Time from event (skip_events+1) start to last event start,
        # divided by number of intervals
        first = timestamps[skip_events][1]
        last = timestamps[-1][1]
        n_intervals = len(timestamps) - 1 - skip_events
        if n_intervals > 0:
            elapsed = (last - first).total_seconds()
            time_per_event = elapsed / n_intervals
            events_measured = n_intervals
    elif nevents and loop_time:
        # Fallback: use total loop time
        time_per_event = loop_time / nevents
        events_measured = nevents

    result = {
        "nevents": nevents,
        "max_rss": max_rss,
        "total_time": total_time,
        "loop_time": loop_time,
        "time_per_event": time_per_event,
        "events_measured": events_measured,
        "events_skipped": skip_events if events_measured else None,
    }
    return result


def get_hostname(args):
    """Find hostname from explicit file or companion file next to log."""
    if args.hostname:
        with open(args.hostname) as f:
            return f.read().strip()
    hostname_path = os.path.join(os.path.dirname(os.path.abspath(args.logfile)), "hostname")
    if os.path.exists(hostname_path):
        with open(hostname_path) as f:
            return f.read().strip()
    return None


def extract_module_timing(logfile_path):
    """Parse TimeReport Module Summary for per-module real time per event.
    Returns dict of {module_name: real_sec_per_event}."""
    modules = {}
    in_section = False
    with open(logfile_path) as f:
        for line in f:
            if "Module Summary ---[Real sec]" in line:
                in_section = True
                continue
            if in_section:
                if "per event" in line and "Name" in line:
                    continue
                if "TimeReport ---" in line or (
                    not line.startswith("TimeReport") and line.strip()
                ):
                    in_section = False
                    continue
                parts = line.split()
                if len(parts) >= 5 and parts[0] == "TimeReport":
                    try:
                        per_event = float(parts[1])
                        name = parts[4]
                        modules[name] = per_event
                    except (ValueError, IndexError):
                        continue
    return modules


def extract_event_summary(logfile_path):
    """Parse TimeReport Event Summary for CPU and Real time per event."""
    with open(logfile_path) as f:
        for line in f:
            m = re.search(
                r"CPU/event\s*=\s*(\S+)\s+Real/event\s*=\s*(\S+)", line
            )
            if m:
                return {
                    "cpu_per_event": float(m.group(1)),
                    "real_per_event": float(m.group(2)),
                }
    return None


def extract_timing(logfile_path):
    """Parse DileptonPlusXProducer timing lines and aggregate across streams."""
    totals = {}
    total_events = 0
    pattern = re.compile(r"DileptonPlusXProducer::timing (.+)")
    with open(logfile_path) as f:
        for line in f:
            m = pattern.search(line)
            if m:
                pairs = m.group(1).split()
                for pair in pairs:
                    if "=" not in pair:
                        continue
                    key, val = pair.split("=", 1)
                    if not val:
                        continue
                    try:
                        val = int(val)
                    except ValueError:
                        continue
                    if key == "events":
                        total_events += val
                    else:
                        totals[key] = totals.get(key, 0) + val
    if total_events == 0:
        return None
    # Convert microseconds to ms/evt
    result = {"events": total_events}
    for key, val in totals.items():
        if key.startswith("n_"):
            result[key] = val / total_events  # counts per event
        else:
            result[key] = val / 1000.0 / total_events  # ms per event
    return result


def format_timing(timing):
    """Format timing dict as human-readable report."""
    if not timing:
        return ""
    lines = []
    lines.append(f"  DileptonPlusXProducer timing ({timing['events']} events, ms/evt):")
    order = [
        ("fillDileptonInfo", "  fillDileptonInfo"),
        ("mmGamma", "  mmGamma"),
        ("buildLLXCandidates", "  buildLLXCandidates"),
        ("fillBtoKllInfo", "    fillBtoKllInfo"),
        ("fillBtoLLhhInfo", "    fillBtoLLhhInfo"),
        ("buildKsCandidate", "    buildKsCandidate"),
        ("fit_kin", "      fit_kin"),
        ("fit_jpsikk", "      fit_jpsikk"),
        ("fit_phill", "      fit_phill"),
        ("fit_jpsiks", "      fit_jpsiks"),
        ("displacement", "      displacement"),
        ("dstar", "  dstar"),
        ("kstar", "  kstar"),
        ("mmm", "  mmm"),
        ("isolation", "  isolation"),
        ("tnp", "  tnp"),
        ("ee", "  ee"),
        ("emu", "  emu"),
        ("hh", "  hh"),
    ]
    for key, label in order:
        if key in timing:
            extra = ""
            if key == "fillBtoKllInfo" and "n_fillBtoKllInfo" in timing:
                extra = f"  ({timing['n_fillBtoKllInfo']:.1f} calls/evt)"
            if key == "fillBtoLLhhInfo" and "n_fillBtoLLhhInfo" in timing:
                extra = f"  ({timing['n_fillBtoLLhhInfo']:.1f} calls/evt)"
            lines.append(f"    {label:30s} {timing[key]:8.1f}{extra}")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Extract performance metrics from cmsRun log")
    parser.add_argument("logfile", help="cmsRun log file")
    parser.add_argument("--rootfile", help="ROOT output file (for file size per event)")
    parser.add_argument("--hostname", type=str, help="hostname file path")
    parser.add_argument("--json", action="store_true", help="output as JSON")
    parser.add_argument("--skip", type=int, default=1,
                        help="number of initial events to skip for timing (default: 1)")

    args = parser.parse_args()

    if not os.path.exists(args.logfile):
        print(f"Error: {args.logfile} not found", file=sys.stderr)
        sys.exit(2)

    metrics = extract_metrics(args.logfile, skip_events=args.skip)
    hostname = get_hostname(args)
    if hostname:
        metrics["hostname"] = hostname

    timing = extract_timing(args.logfile)
    module_timing = extract_module_timing(args.logfile)
    event_summary = extract_event_summary(args.logfile)

    # File size per event
    file_size = None
    file_size_per_event = None
    if args.rootfile and os.path.exists(args.rootfile):
        file_size = os.path.getsize(args.rootfile)
        if metrics["nevents"] and metrics["nevents"] > 0:
            file_size_per_event = file_size / metrics["nevents"]

    if args.json:
        if timing:
            metrics["timing"] = timing
        if module_timing:
            metrics["module_timing"] = module_timing
        if event_summary:
            metrics["event_summary"] = event_summary
        if file_size is not None:
            metrics["file_size"] = file_size
        if file_size_per_event is not None:
            metrics["file_size_per_event"] = file_size_per_event
        print(json.dumps(metrics, indent=2))
    else:
        if hostname:
            print(f"Hostname:   {hostname}")
        if metrics["nevents"] is not None:
            print(f"Events:     {metrics['nevents']}")
        if metrics["time_per_event"] is not None:
            skip = metrics.get("events_skipped", 0)
            measured = metrics.get("events_measured", 0)
            print(f"Time/event: {metrics['time_per_event']:.3f} sec"
                  f" (measured on {measured} events, skipped first {skip})")
        if event_summary:
            print(f"Event loop: CPU/event = {event_summary['cpu_per_event']:.3f} sec"
                  f"  Real/event = {event_summary['real_per_event']:.3f} sec")
        if metrics["max_rss"] is not None:
            print(f"Peak RSS:   {metrics['max_rss']:.0f} kB")
        if file_size is not None:
            size_str = f"{file_size / 1024:.1f} kB"
            if file_size_per_event is not None:
                size_str += f" ({file_size_per_event / 1024:.2f} kB/event)"
            print(f"File size:  {size_str}")
        if module_timing:
            sorted_modules = sorted(module_timing.items(), key=lambda x: -x[1])
            print(f"Module timing (real sec/event):")
            for name, t in sorted_modules:
                if t >= 0.001:
                    print(f"  {name:40s} {t:.4f}")
        if timing:
            print(format_timing(timing))


if __name__ == "__main__":
    main()
