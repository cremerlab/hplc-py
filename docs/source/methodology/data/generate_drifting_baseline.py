#%%
# GENERATE DATA FOR SHOWING THE BASELINE SUBTACTION ALGORITHM
import numpy as np 
import pandas as pd 
import scipy.stats

dt = 0.008
time = np.arange(0, 60, dt)
locs = [15, 19, 40]
amps = [100, 50, 100]
scales = [1, 1.5, 2]
skews = [0, 1, -1]
sig = np.zeros_like(time)
for i in range(len(locs)):
    sig += amps[i] * scipy.stats.skewnorm(skews[i], loc=locs[i], scale=scales[i]).pdf(time)
bg = 20 -  (18/ time.max()) * time
bg += np.sin(time/10)
bg = np.round(bg, decimals=3)
sig = np.round(sig, decimals=3)
df = pd.DataFrame(np.array([sig, bg, bg + sig, time]).T, 
                  columns=['true_signal', 'true_background', 'signal', 'time'])
df.to_csv('./sample_baseline.csv', index=False)
