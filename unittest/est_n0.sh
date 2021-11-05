#!/bin/bash -x
# David Dec 2019
# Tests misc/est_n0

onerun=$(mktemp)
results=$(mktemp)

# generate an impulse with time offset 1
timpulse --n0 1 --f0 100 | c2sim - --modelout - | est_n0 > $results

python3 -c "
import sys; import numpy as np
est_n0 = np.loadtxt(\"$results\")
#ignore first few frames as buffers load up
est_n0 = est_n0[2:]
sys.exit(0) if np.all(est_n0==1) else sys.exit(1)
"
