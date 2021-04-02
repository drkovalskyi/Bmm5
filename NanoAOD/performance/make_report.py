#!/usr/bin/env python3
import sys
import re

if len(sys.argv) != 2:
    print("Usage:\n\t%s <log file>\n" % sys.argv[0])
    sys.exit()

nanoaod_block = dict()
total_time = None
loop_time = None
block_time = dict()
nevents = None
max_rss = None
with open(sys.argv[1]) as logfile:
    timereport_block = None
    for line in logfile:
        match = re.search('^MemoryCheck.*?RSS\s+(\S+)', line)
        if match:
            rss = float(match.group(1))
            if not max_rss or max_rss < rss:
                max_rss = rss
            continue
        match = re.search('TrigReport Events total = (\d+)', line)
        if match:
            nevents = int(match.group(1))
            continue
        match = re.search('TimeReport> Time report complete in\s+(\S+)', line)
        if match:
            total_time = float(match.group(1))
            continue
        match = re.search('Total loop:\s+(\S+)', line)
        if match:
            loop_time = float(match.group(1))
            continue
        match = re.search('^TimeReport\s*\-+\s*(\S.*?)\s*\-\-', line)
        if match:
            timereport_block = match.group(1)
            continue
        if timereport_block == 'Modules in Path: nanoAOD_step':
            match = re.search('^TimeReport\s+([\d\.]+)\s+([\d\.]+)\s+(\S+)', line)
            if match:
                nanoaod_block[match.group(3)] = float(match.group(1))
        match = re.search('^TimeReport\s+([\d\.]+)', line)
        if match:
            if timereport_block not in block_time:
                block_time[timereport_block] = 0
            block_time[timereport_block] += float(match.group(1))
            
print("Max RSS: %f" % max_rss)
print("Total time: %0.1f sec" % total_time)
print("Total event loop time: %0.1f sec" % loop_time)
if loop_time/total_time < 0.90:
    print("WARNING: overhead is %0.0f%% of total time" % (100.*(1-loop_time/total_time)))
    print("Consider running more events for more reliable results")

print("Time per event (total time): %0.3f sec" % (total_time/nevents))
print("Time per event (event loop time): %0.3f sec\n" % (loop_time/nevents))

nanoaod_block_time_per_event = 0
for module, time in nanoaod_block.items():
    nanoaod_block_time_per_event += time
print("nanoAOD_step path time per event: %0.3f sec" % nanoaod_block_time_per_event)

bmm_module_patterns = ['BxToMuMu', 'ForMuonFake', 'BmmMuonId']

print("Bmm modules:")
for module in sorted(nanoaod_block, key=nanoaod_block.get, reverse=True):
    report_module = False
    for pattern in bmm_module_patterns:
        if re.search(pattern, module):
            report_module = True
            break
    if not report_module:
        continue
    time = nanoaod_block[module]
    print("\t%-60s\t%0.4f (%4.1f%%)" % (module, time, 100.*time/nanoaod_block_time_per_event))

show_top = 10
print("Top contributors:")
for module in sorted(nanoaod_block, key=nanoaod_block.get, reverse=True):
    time = nanoaod_block[module]
    if show_top > 0:
        print("\t%-60s\t%0.4f (%4.1f%%)" % (module, time, 100.*time/nanoaod_block_time_per_event))
        show_top -= 1

# print("Block time:")
# for block in sorted(block_time, key=block_time.get, reverse=True):
#    time = block_time[block]
#    print("\t%-30s\t%0.4f (%4.1f%%)" % (block, time, 100.*time/total_time*nevents))
