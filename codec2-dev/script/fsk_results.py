# fsk_results.py
# David Rowe Sep 2018
#
# Reads JSON files from HF data system and prints some summary results

import json
import time
import sys

if len(sys.argv) == 1:
    print("\nusage: %s filename.json\n" % (sys.argv[0]))
    sys.exit(0)
    
filepath = sys.argv[1]
EbNodB_sum = 0.0

with open(filepath) as fp:
    line = fp.readline()
    cnt = 1
    while line:
        data = json.loads(line)
        time_str = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(data['secs']))
        print((time_str, data['EbNodB'],data['frames']))
        EbNodB_sum += float(data['EbNodB'])
        line = fp.readline()
        cnt += 1
    print("Average EbNodB: %4.2f\n" % (EbNodB_sum/cnt))
