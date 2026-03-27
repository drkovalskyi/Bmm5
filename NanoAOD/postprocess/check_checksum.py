import subprocess
import re
import datetime

# Constants
LOCAL_PREFIX = "/eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/532/"
REMOTE_PREFIX = "/store/user/paus/nanoao/532/"
ROOT_LOG_PREFIX = "root://eoscms.cern.ch://eos/cms/store/group/phys_bphys/bmm/bmm6/NanoAOD/532/"
LOCAL_XROOTD = "eoscms.cern.ch"
REMOTE_XROOTD = "xrootd.cmsaf.mit.edu"

def run_checksum(server, path):
    try:
        result = subprocess.run(
            ["xrdfs", server, "query", "checksum", path],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        if result.returncode != 0:
            return None
        match = re.search(r"adler32\s+([a-fA-F0-9]+)", result.stdout)
        return match.group(1) if match else None
    except Exception:
        return None

def extract_relative_path(line):
    if ROOT_LOG_PREFIX in line:
        return line.split(ROOT_LOG_PREFIX)[-1].strip()
    return None

def process_log_file(log_file):
    mismatches = []
    with open(log_file, "r") as f:
        for i, line in enumerate(f, 1):
            rel_path = extract_relative_path(line)
            if not rel_path:
                continue

            print(f"Checking file {i}: {rel_path}")

            local_path = LOCAL_PREFIX + rel_path
            remote_path = REMOTE_PREFIX + rel_path

            local_sum = run_checksum(LOCAL_XROOTD, local_path)
            if local_sum is None:
                print("  Failed to get local checksum.")
                continue

            remote_sum = run_checksum(REMOTE_XROOTD, remote_path)
            if remote_sum is None:
                print("  Failed to get remote checksum.")
                continue

            if local_sum == remote_sum:
                print("  File is good.")
            else:
                print("  File is bad.")
                print(f"    Local checksum : {local_sum}")
                print(f"    Remote checksum: {remote_sum}")
                mismatches.append(local_path)

    return mismatches

def write_report(mismatches):
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    report_filename = f"bad_files_{timestamp}.txt"
    with open(report_filename, "w") as report_file:
        for path in mismatches:
            report_file.write(path + "\n")
    print(f"\nReport saved to: {report_filename}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python check_checksum.py <error_log_file>")
        sys.exit(1)

    log_file = sys.argv[1]
    mismatched_files = process_log_file(log_file)

    if mismatched_files:
        write_report(mismatched_files)
    else:
        print("\nAll files passed checksum comparison.")

