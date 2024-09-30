# %%
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import hplc.quant
import hplc.io
import scipy.stats
import importlib
importlib.reload(hplc.quant)
np.random.seed(666)
dt = 0.01
x = np.arange(0, 40, dt)
sig = np.zeros_like(x)
n_peaks = 6
locs = np.linspace(0.2 * x.max(), 0.8 * x.max(), n_peaks) + \
    np.random.normal(size=n_peaks)
scale = np.random.uniform(1, 2, n_peaks)
skews = np.random.normal(size=n_peaks)
amps = np.random.uniform(10, 500, size=n_peaks)

for i in range(n_peaks):
    sig += amps[i] * \
        scipy.stats.skewnorm(skews[i], loc=locs[i], scale=scale[i]).pdf(x)
plt.plot(sig)

_df = pd.DataFrame(np.array([x, sig]).T, columns=['time', 'signal'])
chrom = hplc.quant.Chromatogram(_df)
chrom.fit_peaks(enforced_locations=locs,
                enforcement_tolerance=1,
                correct_baseline=False)
mix = chrom.unmixed_chromatograms
chrom.show()
# %%
pal = sns.color_palette('Blues_r', n_colors=n_peaks+1)
np.random.shuffle(pal[:-1])

fig, ax = plt.subplots(1, 1, figsize=(6, 2))
ax.axis('off')

for i in range(n_peaks):
    ax.fill_between(x, mix[:, i], color=pal[i], linewidth=0, alpha=0.85)

ax.plot(x, sig, color='k', linewidth=3, zorder=1000)
plt.savefig('./page_logo_dists.pdf')
