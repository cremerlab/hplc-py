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
data = pd.read_csv('./simulated_chromatogram.csv')
peaks = pd.read_csv('./simulated_chromatogram_peaks.csv') 
plt.plot(data['time_min'],  data['intensity_mV'], label='observed chrom')
plt.plot(data['time_min'], data['bg_truth'], label='bg')
plt.legend()

chrom = hplc.quant.Chromatogram('./simulated_chromatogram.csv', bg_subtract=True)
locs = peaks['retention_time'].values
df = chrom.quantify(locs[:2])


#%%
chrom.show()

# %%

# %%

# %%

# %%
chrom.show(locations=True)

chrom.show()
#%%
x = chrom.quantify()
_ = chrom.show()

#%%
signal = data.intensity_mV.copy()

# Ensure positivity of signal
signal *= np.heaviside(signal, 0)

# Compute the LLS operator
tform = np.log(np.log(np.sqrt(signal + 1) + 1) + 1)

# Compute the number of iterations given the window size.
window = 3/0.015

num_iterations = int((window - 1) / 2)
# num_iterations = 10
# Iteratively filter the signal
for i in reversed(range(0,num_iterations)):
    tform_new = tform.copy()
    for j in range(i, len(tform) - i):
        tform_new[j] = min(tform[j], 0.5 * (tform[j+i] + tform[j-i])) 
    tform = tform_new



# Perform the inverse of the LLS transformation
inv_tform = (np.exp(np.exp(tform) - 1) - 1)**2 - 1 


#%%
plt.plot(data['time_min'], data['bg_truth'], label='true background')
plt.plot(data['time_min'], inv_tform, '--', label='inferred bg')

#%%
plt.plot(data['time_min'], data['signal_truth'], label='true background')
plt.plot(data['time_min'], signal - inv_tform, '--', label='inferred bg')