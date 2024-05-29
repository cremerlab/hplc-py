# %%
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import pandas as pd

# ##############################################################################
# TEST DATA FOR PEAK FITTING
# ##############################################################################
# Generate test data for peak fitting
n_peaks = 6
t_bounds = [0, 160]
x_bounds = [0.1 * t_bounds[1], 0.9 * t_bounds[1]]
dt = 0.01
locs = np.linspace(x_bounds[0], x_bounds[1], n_peaks)
scales = np.logspace(-1, np.log10(3), n_peaks)
skews = np.linspace(-5, 5, n_peaks)
amps = np.logspace(2, 4, n_peaks)
x = np.arange(0, t_bounds[1], dt)
chroms = pd.DataFrame([])
peaks = pd.DataFrame([])
_iter = 0
for i, sig in enumerate(scales):
    for j, skew in enumerate(skews):
        signal = np.zeros_like(x)
        _areas = []
        _scales = []
        _skews = []
        _amps = []
        _maxes = []
        for ell, loc in enumerate(locs):
            _pdf = amps[ell] * \
                scipy.stats.skewnorm.pdf(x, skew, loc=loc, scale=sig)
            signal += _pdf
            _areas.append(_pdf.sum())
            _scales.append(sig)
            _skews.append(skew)
            _amps.append(amps[ell])
            _maxes.append(_pdf.max())
        # Generate the peak table
        _peaks = pd.DataFrame(np.array([locs, _areas, _scales, _skews, _amps, _maxes, np.arange(len(locs)) + 1]).T,
                              columns=['retention_time', 'area', 'scale', 'skew', 'amplitude',
                                       'signal_maximum', 'peak_idx'])
        _peaks['iter'] = _iter
        _peaks['peak_idx'] = np.int_(_peaks['peak_idx'])
        _chrom = pd.DataFrame(np.array([x, signal]).T, columns=['x', 'y'])
        _chrom['iter'] = _iter
        peaks = pd.concat([peaks, _peaks], sort=False)
        chroms = pd.concat([chroms, _chrom], sort=False)
        _iter += 1
peaks.to_csv('./test_data/test_fitting_peaks.csv', index=False)
chroms.to_csv('./test_data/test_fitting_chrom.csv', index=False)


# %%
# ##############################################################################
# TEST DATA FOR AUTOMATED PEAK UNMIXING
# ##############################################################################
# Generate test data for peak unmixing
x = np.arange(0, 25, dt)
n_mixes = 20
nudge = 0.2
peak1 = 100 * scipy.stats.norm(8, 1).pdf(x)
chroms = pd.DataFrame([])
peaks = pd.DataFrame([])
for n in range(n_mixes):
    sig = peak1.copy()
    peak2 = 70 * scipy.stats.norm(10.5 + n * 0.2, 1).pdf(x)
    sig += peak2

    # Save chromatogram
    _df = pd.DataFrame(np.array([x, sig]).T, columns=['x', 'y'])
    _df['iter'] = n
    chroms = pd.concat([chroms, _df])

    # Save the peak info
    _df = pd.DataFrame(np.array([[8, 10.5 + n * 0.2], [1, 1], [0, 0],
                                [100, 70], [peak1.sum(), peak2.sum()],
                                [peak1.max(), peak2.max()],
                                [1, 2], [n, n]]).T,
                       columns=['retention_time', 'scale', 'skew',
                                'amplitude', 'area', 'signal_maximum',
                                'peak_idx', 'iter'])
    peaks = pd.concat([peaks, _df])
chroms.to_csv('./test_data/test_unmix_chrom.csv', index=False)
peaks.to_csv('./test_data/test_unmix_peaks.csv', index=False)


# %%
# ##############################################################################
# TEST DATA FOR BACKGROUND SUBTRACTION
# ##############################################################################
# Generate a noisy background to test the background subtraction algorithm
np.random.seed(666)
n_peaks = 10
x = np.arange(0, 40, dt)
bg = np.zeros_like(x)
amps = np.abs(np.random.normal(200, 30, n_peaks))
loc = np.abs(np.random.uniform(-10, 40, n_peaks))
scale = np.abs(np.random.normal(10, 2, n_peaks))
for i in range(n_peaks):
    bg += amps[i] * scipy.stats.norm(loc[i], scale[i]).pdf(x)


for i in range(5):
    if i == 0:
        sig = 100 * scipy.stats.norm(5, 0.1).pdf(x)
    else:
        sig += 100 * scipy.stats.norm(5 + 7 * i, 0.1).pdf(x)

# Add strong candidate signal
noise = np.random.normal(size=len(x))
df = pd.DataFrame(np.array([x, bg, sig + bg]).T, columns=['x', 'bg', 'y'])
df.to_csv('./test_data/test_SNIP_chrom.csv', index=False)

# %%
# ##############################################################################
# TEST DATA FOR MANUAL LOCATION PEAK UNMIXING
# ##############################################################################
x = np.arange(0, 25, dt)
n_mixes = 20
nudge = 0.2
peak1 = 100 * scipy.stats.norm(8, 1).pdf(x)
chroms = pd.DataFrame([])
peaks = pd.DataFrame([])
amps = np.linspace(10, 100, n_mixes)
for n in range(n_mixes):
    sig = peak1.copy()
    loc = 11
    scale = 2
    amp = 1.5 * amps[n]
    peak2 = amp * scipy.stats.norm(loc, scale).pdf(x)
    sig += peak2
    # Save chromatogram
    _df = pd.DataFrame(np.array([x, sig]).T, columns=['x', 'y'])
    _df['iter'] = n
    chroms = pd.concat([chroms, _df])

    # Save the peak info
    _df = pd.DataFrame(np.array([[8, loc], [1, scale], [0, 0],
                                [100, amp], [peak1.sum(), peak2.sum()],
                                [peak1.max(), peak2.max()],
                                [1, 2], [n, n]]).T,
                       columns=['retention_time', 'scale', 'skew',
                                'amplitude', 'area', 'signal_maximum', 'peak_idx', 'iter'])
    peaks = pd.concat([peaks, _df])
chroms.to_csv('./test_data/test_manual_unmix_chrom.csv', index=False)
peaks.to_csv('./test_data/test_manual_unmix_peaks.csv', index=False)

# %%
# ##############################################################################
# TEST DATA FOR MANUAL LOCATION PEAK ASSIGNMENT
# ##############################################################################
# Generate data with a very shallow peak that would not normally be detected
# but can be identified if the manual position is included.
x = np.arange(0, 80, dt)
sig1 = 100 * scipy.stats.norm(10, 1).pdf(x)
sig2 = 10 * scipy.stats.norm(50, 3).pdf(x)
sig = sig1 + sig2
df = pd.DataFrame(np.array([x, sig]).T, columns=['x', 'y'])
df.to_csv('./test_data/test_shallow_signal_chrom.csv', index=False)
peak_df = pd.DataFrame(np.array([[10, 50], [1, 3], [0, 0], [100, 10],
                                 [sig1.sum(), sig2.sum()],
                                 [sig1.max(), sig2.max()], [1, 2]]).T,
                       columns=['retention_time', 'scale', 'skew',
                                'amplitude', 'area', 'signal_maximum', 'peak_idx'])
peak_df.to_csv('./test_data/test_shallow_signal_peaks.csv', index=False)

# %%
# ##############################################################################
# TEST DATA FOR CHROMATOGRAM RECONSTRUCTION REPORTING
# ##############################################################################
# Generate data with an overlapping peak pair, a very shallow peak, and two
# background windows with noise

x = np.arange(0, 150, dt)
sig1 = 100 * scipy.stats.norm(15, 1).pdf(x)
sig2 = 100 * scipy.stats.norm(19, 2).pdf(x)
sig3 = 20 * scipy.stats.norm(80, 3).pdf(x)
sig4 = 100 * scipy.stats.norm(110, 1).pdf(x)


sig = sig1 + sig2 + sig3 + sig4
sig[0:int(7/dt)] += np.random.normal(0, 0.01, len(sig[0:int(7/dt)]))
df = pd.DataFrame(np.array([x, sig]).T, columns=['x', 'y'])
peak_df = pd.DataFrame(np.array([[10, 25], [1, 3], [0, 0], [100, 10],
                                 [sig1.sum(), sig2.sum()],
                                 [sig1.max(), sig2.max()],
                                 [1, 2]]).T,
                       columns=['retention_time', 'scale', 'skew',
                                'amplitude', 'area', 'signal_maximum', 'peak_idx'])
df.to_csv('./test_data/test_assessment_chrom.csv', index=False)
peak_df.to_csv('./test_data/test_assessment_peaks.csv', index=False)

# Generate the scoring table
score_df = pd.DataFrame(np.array([1, 2, 3, 1, 2]).T, columns=['window_id'])
score_df['window_type'] = ['interpeak',
                           'interpeak', 'interpeak', 'peak', 'peak']
score_df['status'] = ['needs review', 'invalid', 'valid', 'invalid', 'valid']
score_df.to_csv('./test_data/test_assessment_scores.csv', index=False)

# %%
# ##############################################################################
# TEST DATA FOR  MANY PEAKS IN A WINDOW
# ##############################################################################
x = np.arange(0, 100, dt)
_locs = np.linspace(20, 55, 10)
locs = [10]
for _l in _locs:
    locs.append(_l)
locs.append(70)
sig = np.zeros_like(x)
for i in range(len(locs)):
    sig += 100 * scipy.stats.norm(locs[i], 1).pdf(x)
df = pd.DataFrame(np.array([x, sig]).T, columns=['x', 'y'])
df.to_csv('./test_data/test_many_peaks.csv', index=False)

# %%

# ##############################################################################
# TEST DATA FOR SPECIFIC BOUND SETTING
# ##############################################################################
x = np.arange(0, 20, dt)
amps = np.linspace(100, 1000, 8)
sig1_props = [1000, -5, 12, 2]
sig1 = sig1_props[0] * scipy.stats.skewnorm(sig1_props[1], loc=sig1_props[2],
                                            scale=sig1_props[3]).pdf(x)
chrom_dfs = pd.DataFrame([])
peak_dfs = pd.DataFrame([])
for i, a in enumerate(amps):
    sig2 = a * scipy.stats.norm(11, 1).pdf(x)
    sig = sig1 + sig2
    _df = pd.DataFrame(np.array([x, sig]).T, columns=['x', 'y'])
    _df['iter'] = i + 1
    chrom_dfs = pd.concat([chrom_dfs, _df], sort=False)
    _df = pd.DataFrame(np.array([[sig1_props[0], a],
                                 [sig1_props[1], 0],
                                 [sig1_props[2], 11],
                                 [sig1_props[3], 1],
                                 [sig1.sum(), sig2.sum()],
                                 [sig1.max(), sig2.max()],
                                 [2, 1]]).T,
                       columns=['amplitude', 'skew',
                                'retention_time', 'scale',
                                'area', 'signal_maximum',
                                'peak_id'])
    _df['iter'] = i + 1
    peak_dfs = pd.concat([peak_dfs, _df], sort=False)
chrom_dfs.to_csv('./test_data/test_bounding_chroms.csv', index=False)
peak_dfs.to_csv('./test_data/test_bounding_peaks.csv', index=False)


# %%
# ##############################################################################
# TEST DATA FOR INTEGRATION WINDOWING
# ##############################################################################
x = np.arange(0, 20, dt)
sig = 500 * scipy.stats.skewnorm(0, loc=10, scale=1.5).pdf(x)
sig = np.round(sig, decimals=3)


# Save the chromatogram
chrom_df = pd.DataFrame(np.array([x, sig]).T, columns=['time', 'signal'])
chrom_df.to_csv('./test_data/test_integration_window_chrom.csv', index=False)

# Integrate the area over a range of windows
windows = [[x[0], x[-1]], [2.5, 17.5], [5, 15], [7.5, 12.5], [9, 11]]
df = pd.DataFrame([])
for i, w in enumerate(windows):
    area = chrom_df[(chrom_df['time'] >= w[0]) &
                    (chrom_df['time'] <= w[1])]['signal'].sum()
    _df = pd.DataFrame(
        {'t_start': w[0], 't_end': w[1], 'area': area, 'window': i+1}, index=[0])
    df = pd.concat([df, _df], sort=False)
df.to_csv('./test_data/test_integration_window_areas.csv', index=False)
