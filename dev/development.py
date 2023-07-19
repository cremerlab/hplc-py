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
plt.plot(data['time_min'],  data['intensity_mV'])
plt.plot(data['time_min'], data['bg_truth'])



#%%
# Test the bg identification algorithm. 
# Heaviside transformation for positivity

intensity = data['intensity_mV'].values.copy()
pos = data['intensity_mV'].values * (data['intensity_mV'].values >= 0)
# Log transformation for application of filter 
tform = np.log(np.log(np.sqrt(pos + 1) + 1 ) + 1)

window_t = 100 
dt = np.mean(np.diff(data['time_min'].values))
window = int(window_t / dt)

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

# %%

# %%

#%%
plt.plot(data['time_min'], data['bg_truth'])
plt.plot(data['time_min'], _corr)
# %%

# %%
