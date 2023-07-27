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


#%%
plt.plot(np.diff(data['intensity_mV']))

#%%
_df = chrom.df
f = np.sum(chrom.deconvolved_peaks, axis=1)
len(chrom.window_props)