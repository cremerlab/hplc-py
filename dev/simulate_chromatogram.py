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
stdevs = np.random.uniform(0.3, 0.8, n_peaks)
# stdevs[-1] = 0.9
skews = np.random.normal(0, 1, n_peaks)
amps = 10**np.random.normal(3.5, 0.5, n_peaks)
dt = 0.015 # In minutes 
time = np.arange(0, 50, dt)
signal = np.zeros_like(time)
areas = []
# Generate the true signals
for i in range(len(locs)):
    norm = scipy.stats.skewnorm.pdf(time, skews[i], loc=locs[i], scale=stdevs[i])
    peak = amps[i] * norm
    areas.append(peak.sum())
    signal += peak

# Generate a drifting baseline
# bg = np.zeros_like(time)
bg = 600 - 20 * time
# for i in range(100):
    # bg += np.random.uniform(1, 30) * scipy.stats.norm(np.random.uniform(5, 30), np.sqrt(np.random.normal(1, 0.5)**2)).pdf(time)
# bg += 200 * np.sin(time/ 50)
# Generate a noise profile
noise = np.random.normal(0, 2, len(time))

# Generate the simulated data
shift = 50
chrom_sim = bg + noise + signal + shift
# chrom_sim -= np.min(chrom_sim)
df = pd.DataFrame(np.array([time, chrom_sim, bg + noise + shift, signal]).T, columns=['time_min', 'intensity_mV', 'bg_truth', 'signal_truth'])
df = df[df['time_min'] <= 30]
df.to_csv('./simulated_chromatogram.csv', index=False)

# Generate the truth peak table. 
peak_idx = np.arange(len(locs))
peak_table = pd.DataFrame(np.array([amps, areas, peak_idx, 
                                    locs, stdevs, skews]).T,
                         columns=['amplitude', 'area', 'peak_idx', 'retention_time', 'scale', 'skew'])
peak_table['sample'] = 'simulated data'
peak_table.to_csv('./simulated_chromatogram_peaks.csv', index=False)

# %%
