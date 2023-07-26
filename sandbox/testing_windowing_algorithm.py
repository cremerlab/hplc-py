#%%
import numpy as np 
import matplotlib.pyplot as plt
import scipy.stats
import scipy.signal
import pandas as pd

# Generate a signal
time = np.linspace(0, 10, 1000)
dt = np.diff(time)[1]
sig = np.zeros_like(time)
for i in range(3):
    sig += 100 * scipy.stats.norm(loc=3 + i * 2, scale=0.5).pdf(time)

plt.plot(time, sig)

#%%
locs = [3.1, 5.1, 7.1]
plt.plot(time, sig)
plt.vlines(locs, 0, sig.max(), 'r')

#%%
# With prescribed locations find the windows
indices = np.array([(np.abs(loc - time)).argmin() for loc in locs])
_left = indices - 5/dt 
_right = indices + 5/dt


ranges = []
buffer = 100
for l, r in zip(_left, _right):
    if (l - buffer) <= 0:
        l = 0
    else:
        l -= buffer
    if (r + buffer) <= 0:
        r = len(time)
    else:
        r += buffer
    ranges.append(np.arange(np.round(l), np.round(r), 1))
valid = [True] * len(ranges)
for i, r1 in enumerate(ranges):
    for j, r2 in enumerate(ranges):
        if i != j:
            if set(r2).issubset(r1):
                valid[j] = False
valid

ranges = [r for i, r in enumerate(ranges) if valid[i] is True]

_df = pd.DataFrame(np.array([time, sig]).T, columns=['x', 'y'])
_df['idx'] = 0
_df['time_idx'] = np.arange(len(_df))
for i, r in enumerate(ranges):
    _df.loc[_df['time_idx'].isin(r), 'idx'] = int(i + 1)

_df = _df[_df['idx'] > 0]
for g, d in _df.groupby(['idx']):
    if g == 1:
        plt.plot(d['y'])