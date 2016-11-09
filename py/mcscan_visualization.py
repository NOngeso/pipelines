#!/usr/bin/python3

import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
import seaborn as sns
sns.set_palette("deep", desat=.6)
sns.set_context(rc={"figure.figsize": (8, 4)})

file_mcscan = sys.argv[1] #'/data2/k821209/KimKH/Gm2Gm.collinearity.kaks'
ks_array = np.array([])
for line in open(file_mcscan):
    if line[0] == '#':
        continue
    cell = line.strip().split('\t')
    fKs = float(cell[5])
    ks_array = np.append(ks_array,fKs)

data = ks_array
density = gaussian_kde(data)
xs = np.linspace(0,3,200)

density.covariance_factor = lambda : .25
density._compute_covariance()
aa = np.argmax(density(xs)) # find the position showing the max value
print 'Ks:',round(xs[aa],3),'density:',round(density(xs)[aa],3)

max_x = xs[aa]
max_y = density(xs)[aa]

plt.plot(xs,density(xs))

plt.annotate('Peak at Ks %f'%round(xs[aa],3), xy=(max_x, max_y), xytext=(max_x+0.3, max_y+0.1),
            arrowprops=dict(facecolor='black', shrink=0.05),
            )

plt.ylabel('Density')
plt.xlabel('Ks value')
plt.savefig('%s.ks_dist.png'%sys.argv[2], dpi=300)
plt.show()
