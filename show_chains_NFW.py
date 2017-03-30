import matplotlib.pyplot as plt
import numpy as np
import sys

filename = sys.argv[1]
chain = np.load(filename)

param_names = [r"$M_{200m}$", r"$c_{200m}$"]
plt.figure(figsize = (10,10))
for i in range(chain.shape[2]):
    plt.subplot(2,1,i+1)
    c = chain[:,:,i]
    for j in range(c.shape[0]):
        plt.plot(c[j], "-k", alpha = 0.2)
    plt.ylabel(param_names[i])

plt.show()
