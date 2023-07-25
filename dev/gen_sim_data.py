#%%
import numpy as np 
import scipy.stats 
import matplotlib.pyplot as plt
import pandas as pd
np.random.seed(666)

# Generate synthetic data that looks realistic 
n_peaks = 6
locs = np.linspace(5, 20, n_peaks)
locs[-2:] += 8
locs[-2] -= 1
locs = np.append(locs, [35, 36.4])
stdevs = np.random.uniform(0.3, 0.8, n_peaks)
stdevs[-1] = 0.7
stdevs[-2] = 2 
stdevs = np.append(stdevs, [0.5, 1])
skews = np.random.normal(0, 4, n_peaks)
skews[-2] = 6
skews = np.append(skews, [0, 0])
amps = 10**np.random.normal(3.5, 0.5, n_peaks + 2)
amps[-2] *= 3
amps[-3] *= 2 

dt = 0.015 # In minutes 
time = np.arange(0, 50, dt)
signal = np.zeros_like(time)
areas = []
# Generate the true signals
for i in range(len(locs)):
    norm = scipy.stats.skewnorm.pdf(time, skews[i], loc=locs[i], scale=stdevs[i])
    peak = amps[i] * norm
    areas.append(peak.sum())
    if i == 5:
        signal += peak

# Generate a drifting baseline
bg = np.zeros_like(time)
for i in range(100):
    bg += np.random.uniform(1, 30) * scipy.stats.norm(np.random.uniform(5, 30), np.sqrt(np.random.normal(1, 0.5)**2)).pdf(time)
bg += 200 * np.sin(time/ 20)
# Generate a noise profile
noise = np.random.normal(0, 2, len(time))
plt.plot(signal)
#%%
# Generate the simulated data
shift = 100
chrom_sim = bg + noise + signal + shift
# chrom_sim -= np.min(chrom_sim)
df = pd.DataFrame(np.array([time, chrom_sim]).T, columns=['R.Time (min)', 'Intensity'])
# df = df[df['time_min'] <= 30]
df.to_csv('./example.txt', index=False)

# Generate the truth peak table. 
peak_idx = np.arange(len(locs))
peak_table = pd.DataFrame(np.array([amps, areas, peak_idx, 
                                    locs, stdevs, skews]).T,
                         columns=['amplitude', 'area', 'peak_idx', 'retention_time', 'scale', 'skew'])
peak_table['sample'] = 'simulated data'
peak_table.to_csv('./simulated_chromatogram_peaks.sv', index=False)

# %%
import matplotlib.pyplot as plt
plt.plot(chrom_sim)