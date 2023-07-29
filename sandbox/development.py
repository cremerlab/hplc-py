#%%
import pandas as pd
import numpy as np 
import hplc.quant
import matplotlib.pyplot as plt 
import seaborn as sns
sns.set()
import imp
imp.reload(hplc.quant)

# Load the simulated data and ground t
# data = pd.read_csv('./sample_chromatogram.txt')
# data = pd.read_csv('./simulated_chromatogram.csv') 
data = pd.read_csv('test_shallow_signal_chrom.csv')
chrom = hplc.quant.Chromatogram(data, cols={'time':'x','signal':'y'})
# chrom.crop([10, 20])
_ = chrom.fit_peaks()

_ = chrom.show()
_ = chrom.assess_fit()
_

#%%
df = chrom.window_df
bg_windows = df[df['window_id']==0]
tidx = bg_windows['time_idx'].values
# plt.plot(bg_windows['time_idx'], 'o')

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
#%%
