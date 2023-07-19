#%%
import numpy as np 
import scipy.stats 
import matplotlib.pyplot as plt
import pandas as pd
np.random.seed(666)
# Generate synthetic data that looks realistic 
n_peaks = 6
locs = np.linspace(5, 18, n_peaks)
stdevs = np.random.uniform(0.1, 0.8, n_peaks)
stdevs[-1] = 0.9
skews = np.random.normal(0, 1, n_peaks)
amps = 10**np.random.normal(3, 0.2, n_peaks)
dt = 0.015 # In minutes 
time = np.arange(0, 30, dt)
signal = np.zeros_like(time)
areas = []
# Generate the true signals
for i in range(len(locs)):
    norm = scipy.stats.skewnorm.pdf(time, skews[i], loc=locs[i], scale=stdevs[i])
    peak = amps[i] * norm
    areas.append(peak.sum())
    signal += peak

# Generate a drifting baseline
# bg = 100 - 10 * time
bg = np.zeros_like(time)
for i in range(20):
    bg += np.random.uniform(10, 150) * scipy.stats.norm(np.random.uniform(-10, 70), np.sqrt(np.random.normal(4, 1)**2)).pdf(time)
# Generate a noise profile
noise = np.random.normal(0, 1, len(time))

# Generate the simulated data
shift = 0
chrom_sim = bg + noise + signal 
chrom_sim -= np.min(chrom_sim)
df = pd.DataFrame(np.array([time, chrom_sim, bg, signal]).T, columns=['time_min', 'intensity_mV', 'bg_truth', 'signal_truth'])
df.to_csv('./simulated_chromatogram.csv', index=False)

# Generate the truth peak table. 
peak_idx = np.arange(len(locs))
peak_table = pd.DataFrame(np.array([amps, areas, peak_idx, 
                                    locs, stdevs, skews]).T,
                         columns=['amplitude', 'area', 'peak_idx', 'retention_time', 'scale', 'skew'])
peak_table['sample'] = 'simulated data'
peak_table.to_csv('./simulated_chromatogram_peaks.csv', index=False)
