#%%
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 
import hplc.quant
import scipy.stats

locs = [10, 12.5]
scale = [ 0.7, 0.85]
amps = [5, 4]
x = np.linspace(0, 20, 1000)
pdf = np.zeros_like(x)
for i in range(len(locs)):
    pdf += amps[i] * scipy.stats.norm(locs[i], scale[i]).pdf(x)
df = pd.DataFrame(np.array([x, pdf]).T, columns=['time_min', 'intensity_mV']) 
chrom = hplc.quant.Chromatogram(df)
_ = chrom.quantify()

mixture = chrom.mix_array
fig, ax = plt.subplots(1, 1, figsize=(2,1))
ax.axis(False)
ax.plot(x, pdf, 'k-')
ax.fill_between(x, mixture[:, 0], zorder=1000)
ax.fill_between(x, mixture[:, 1],zorder=1000)
plt.savefig('./logo.pdf')
# plt.plot(x, pdf)
# plt.xlim([5, 18])
# %%
