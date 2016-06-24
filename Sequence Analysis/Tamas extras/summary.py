import glob
import os
import numpy as np
import re
import sys

results = [[] for x in range(7)]

mintrial = -1;
maxtrial = -1;

if len(sys.argv) > 1:
    match = re.match(r'(\d+)-(\d+)',sys.argv[1])
    mintrial = int(match.group(1))
    maxtrial = int(match.group(2))

for file in glob.glob("Out/stdout_*.txt"):
    match = re.search('_(\d+).txt',file)
    trial = int(match.group(1))
    if mintrial > 0 and (trial < mintrial or trial > maxtrial):
        continue
    data = open(file,'r').read()
    # get outputted values
    matches = re.findall(r'Final mutations \d+: (\d+.\d+)',data)
    # convert to numerical values
    vals = [float(x) for x in matches]
    # and save them
    for i,x in enumerate(vals):
        results[i].append(x)

# now go through and display
coverages = [2, 4, 8, 15, 25, 50, 100]
for i,r in enumerate(results):
    print 'Seq %d: %0.2f (%d)' % (coverages[i],np.mean(r),len(r))
