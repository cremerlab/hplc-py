# %%

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
from hplc.quant import Chromatogram

locs = [17, 22, 50]
scales = [2, 2, 3]
skews = [-0.5, 1, 0]
amps = [10, 100, 50]
dt = 0.008
x = np.arange(0, 70, dt)
signal = np.zeros_like(x)
for i in range(len(locs)):
    signal += amps[i] * \
        scipy.stats.skewnorm(skews[i], loc=locs[i], scale=scales[i]).pdf(x)


# Generate a background signal
bg = 10 - (1/x.max()) * x
bg *= np.cos(x/20) + 2
bg = np.round(bg, decimals=4)
signal = np.round(signal, decimals=4)
df = pd.DataFrame(np.array([x, signal + bg, signal, bg]).T,
                  columns=['time', 'measured_signal', 'true_signal', 'true_background'])
df.to_csv('../data/strong_baseline_drift.csv', index=False)

bg = 1 - (1/x.max()) * x
bg *= np.cos(x/30) + 2
bg = np.round(bg, decimals=4)
signal = np.round(signal, decimals=4)
df = pd.DataFrame(np.array([x, signal + bg, signal, bg]).T,
                  columns=['time', 'measured_signal', 'true_signal', 'true_background'])
df.to_csv('../data/mild_baseline_drift.csv', index=False)
plt.plot(df['measured_signal'])
# %%
windows = np.arange(0.5, 10.5, 1)
fig, ax = plt.subplots(1, 1)
for w in windows:
    chrom = Chromatogram(
        df, cols={'time': 'time', 'signal': 'measured_signal'})
    df_sub = chrom.correct_baseline(window=w, return_df=True, verbose=False)
    plt.plot(df_sub['time'], df_sub['estimated_background'], label=w)
plt.plot(df['time'], df['true_background'], 'r--', label='true background')
plt.legend()
