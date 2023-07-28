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
# data = pd.read_csv('./sample_chromatogram.txt')
# data = pd.read_csv('./simulated_chromatogram.csv') 
data = pd.read_csv('test_shallow_signal_chrom.csv')
chrom = hplc.quant.Chromatogram(data, cols={'time':'x','signal':'y'})
# chrom.crop([10, 20])
chrom.fit_peaks(prominence=0.5)

chrom.show()
chrom.assess_fit()

#%%
bg_time_id = chrom.window_df[chrom.window_df['window_id']==0].time_id.values
split_inds = np.nonzero(np.diff(bg_time_id) - 1)[0]
chrom.window_df[chrom.window_df['window_id']==0]['x'].values[split_inds+1]

# %%
w2 = chrom.window_df[chrom.window_df['window_id']==2]
plt.plot(w2['time_min'], w2['intensity_mV']+1)
plt.plot(chrom.df['time_min'], np.sum(chrom.unmixed_chromatograms, axis=1))

# %%
recon = np.sum(chrom.unmixed_chromatograms, axis=1)
res = np.sqrt((recon - chrom.df['intensity_mV'])**2)
np.sum(res) / np.sum(chrom.df['intensity_mV'])
# %%

# %%
import scipy.signal
intensity = chrom.df.y
norm_int = (intensity - intensity.min()) / (intensity.max() - intensity.min())
peaks, _ = scipy.signal.find_peaks(norm_int, prominence=0.5)
_widths, _, _, _ = scipy.signal.peak_widths(intensity, peaks, rel_height=0.5)
_, _, left, right = scipy.signal.peak_widths(intensity, peaks, rel_height=1)
print(left, right)
ranges = []
for l, r in zip(left, right):
    _range = np.arange(int(l - 100), int(r + 100), 1)
    _range = _range[(_range >= 0) & (_range <= len(norm_int))]
    ranges.append(_range)
plt.plot(intensity)
plt.vlines(_range[0], 0,intensity.max())
plt.vlines(_range[-1], 0,intensity.max())
# %%



# %%

# %%
chrom.window_df
# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%
np.sqrt(data['y']).sum()

# %%
chrom.peaks
# %%
chrom.window_props