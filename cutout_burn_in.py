import numpy as np
import sys, os
import re

dir = sys.argv[1]
burn_in = int(sys.argv[2])

filename = os.path.join(dir, "chains.npy")
chain = np.load(filename)
chain = chain[:, burn_in:, :]
root, ext = os.path.splitext(filename)
filename = root + ".burnin%d" % burn_in
np.save(filename, chain)
print "file written in: %s.npy" % filename

filename = os.path.join(dir, "lnprobabilities.npy")
lnprobabilities = np.load(filename)
lnprobability = lnprobabilities[:, burn_in:]
root, ext = os.path.splitext(filename)
filename = root + ".burnin%d" % burn_in
np.save(filename, lnprobability)
print "file written in: %s.npy" % filename


