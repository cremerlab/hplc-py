#%%
import pandas as pd
import numpy as np 
import hplc.quant
import matplotlib.pyplot as plt 
import seaborn as sns
sns.set()

# Load the simulated data and ground truth
data = pd.read_csv('./simulated_chromatogram.csv')
peaks = pd.read_csv('./simulated_chromatogram_peaks.csv') 
plt.plot(data['time_min'],  data['intensity_mV'], label='observed chrom')
plt.plot(data['time_min'], data['bg_truth'], label='bg')
plt.legend()



#%%
chrom = hplc.quant.Chromatogram(data, peak_width=2)
chrom.quantify(verbose=True, prominence=0.005, rel_height=1.0)

#%%
chrom.show()

# %%
chrom = hplc.quant.Chromatogram(data, peak_width=1)
df = chrom._bg_subtract(return_df=True, window=10)

#%%
plt.plot(data['time_min'],  data['intensity_mV'], label='observed chrom')
plt.plot(df['time_min'],  df['intensity_mV'], label='subtracted(??)')
plt.plot(df['time_min'], df['bg_truth'], label='bg')
plt.legend()
# %%

# %%

# %%

#%%

# Test the bg identification algorithm. 
intensity = data['intensity_mV'].values.copy()
pos = intensity * (intensity >= 0)
# Log transformation for application of filter 
tform = np.log(np.log(np.sqrt(pos + 1) + 1 ) + 1)

window_t = 3
dt = np.mean(np.diff(data['time_min'].values))
window = int(window_t / dt)
# window = 10
its = []
for i in range(0, window):
    tform_new = np.zeros_like(tform)
    for j in range(i, len(tform) - i):
        tform_new[j] = min(tform[j], 0.5 * (tform[j+i] + tform[j-i])) 
    tform = tform_new

_corr = (np.exp(np.exp(tform)-1)-1)**2 - 1
# plt.plot(pos, 'r-')
# plt.plot(corr, 'b-')
plt.plot(data['time_min'], data['signal_truth'])
plt.plot(data['time_min'], pos - _corr, '--', label='corrected')
plt.legend()


# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

#%%
plt.plot(data['time_min'], data['bg_truth'])
plt.plot(data['time_min'], _corr)
# %%

# %%
