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
# peaks = pd.read_csv('./simulated_chromatogram_peaks.csv') 
chrom = hplc.quant.Chromatogram(data, cols={'time':'time_min','signal':'intensity_mV'})
chrom.crop([10, 20])
chrom.fit_peaks(buffer=100)
chrom.show()
chrom.assess_fit(tol=1E-3)

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

# %%

# %%

# %%

# %%
np.sqrt(data['y']).sum()

# %%

# %%
