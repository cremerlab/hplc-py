#%%
import pandas as pd
import numpy as np 
import hplc.quant
import matplotlib.pyplot as plt 
import seaborn as sns
sns.set()
import imp
imp.reload(hplc.quant)

# Load the simulated data and ground truth
data = pd.read_csv('./sample_chromatogram.txt')
chrom = hplc.quant.Chromatogram(data, cols={'time':'time_min','signal':'intensity_mV'})
chrom.crop([10, 20])
chrom.fit_peaks()
chrom.show()
if len(bg_windows) > 0:
    split_inds = np.nonzero(np.diff(bg_windows['time_idx']) - 1)[0]
    if split_inds[0] != 0:
        n_bg_windows = len(split_inds) + 1
    else:
        n_bg_windows = 1

    # Add the indices for the ranges
    split_inds += 1 
    split_inds = np.insert(split_inds, 0, 0) 
    split_inds = np.append(split_inds, len(tidx)-1)
    bg_ranges = [bg_windows.iloc[np.arange(split_inds[i], split_inds[i+1])]['time_idx'].values for i in range(len(split_inds) - 1)]
#%%

plt.plot(bg_ranges[0], np.ones_like(bg_ranges[0]))
plt.plot(bg_ranges[1], np.ones_like(bg_ranges[1]))
plt.plot(df['y'], '.')
#%%
win_1 = bg_windows.iloc[bg_ranges[0]].time_idx
win_1
#%%
for g, d in chrom.window_df.groupby(['window_id', 'window_type']):
        plt.plot(d['y'], '.', label=g)
plt.legend()