# %%
import importlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
import hplc.quant
importlib.reload(hplc.quant)

# simulate a dataset with a negative baseline
dt = 0.01
x = np.arange(0, 40, dt)
sig1 = 1000 * scipy.stats.norm(6, 1).pdf(x)
sig2 = -100 * scipy.stats.norm(30, 1).pdf(x)
sig3 = 30 * scipy.stats.norm(20, 2).pdf(x)
sig5 = 100 * scipy.stats.norm(10, 1).pdf(x)
sig = np.round(sig1 + sig2 + sig3 + sig5, decimals=9) - 10

df = pd.DataFrame(np.array([x, sig]).T, columns=['x', 'y'])

chrom = hplc.quant.Chromatogram(df, cols={'time': 'x', 'signal': 'y'})
chrom.correct_baseline(window=5)
chrom.show()
chrom.fit_peaks()
chrom.show()
chrom.assess_fit()

# %%
plt.plot(chrom.df['y'], 'b')
plt.plot(chrom.df['y_corrected'], 'r')
# %%
df = pd.read_csv('./test_SNIP_chrom.csv')
df['y'] -= 20
df['bg'] -= 20
chrom = hplc.quant.Chromatogram(df, cols={'time': 'x', 'signal': 'y'})
chrom.correct_baseline(window=0.5)
chrom.show()


# %%
np.median(df['y'])

# %%
for g, d in chrom.window_df.groupby(['window_type', 'window_id']):
    plt.plot(d['y_corrected'], label=g)
plt.legend()
