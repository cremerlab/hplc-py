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
df = chrom.quantify()


#%%
chrom.show()