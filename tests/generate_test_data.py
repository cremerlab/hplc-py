#%%
import numpy as np 
import scipy.stats 
import pandas as pd

# Generate test data for peak fitting
n_peaks = 5
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
        for ell, loc in enumerate(locs):
            _pdf = amps[ell] * scipy.stats.skewnorm.pdf(x, skew, loc=loc, scale=sig)
            signal += _pdf
            _areas.append(_pdf.sum())
            _scales.append(sig)
            _skews.append(skew)
            _amps.append(amps[ell])
        # Generate the peak table
        _peaks = pd.DataFrame(np.array([locs, _areas, _scales, _skews, _amps, np.arange(len(locs)) + 1]).T,
                              columns=['retention_time', 'area', 'scale', 'skew', 'amplitude', 'peak_idx'])
        _peaks['iter'] = _iter
        _peaks['peak_idx'] = np.int_(_peaks['peak_idx'])
        _chrom = pd.DataFrame(np.array([x, signal]).T, columns=['x', 'y'])
        _chrom['iter'] = _iter
        peaks = pd.concat([peaks, _peaks], sort=False)
        chroms = pd.concat([chroms, _chrom], sort=False)
        _iter += 1  
peaks.to_csv('../tests/test_fitting_peaks.csv', index=False) 
chroms.to_csv('../tests/test_fitting_chrom.csv', index=False) 

#%%
# Generate test data for peak unmixing
x = np.arange(0, 25, dt)
n_mixes = 20
nudge=0.2
peak1 = 100 * scipy.stats.norm(8, 1).pdf(x)
chroms = pd.DataFrame([])
peaks = pd.DataFrame([])
for n in range(n_mixes):
    sig = peak1.copy()
    peak2 = 70 * scipy.stats.norm(10.5+ n * 0.2, 1).pdf(x)
    sig += peak2

    # Save chromatogram
    _df = pd.DataFrame(np.array([x, sig]).T, columns=['x', 'y'])
    _df['iter'] = n
    chroms = pd.concat([chroms, _df])


    # Save the peak info
    _df = pd.DataFrame(np.array([[8, 10.5 + n * 0.2], [1, 1], [0, 0], 
                                [100, 70], [peak1.sum(), peak2.sum()],
                                [1, 2], [n, n]]).T, 
                                columns=['retention_time', 'scale', 'skew', 
                                         'amplitude', 'area', 'peak_idx', 'iter'])
    peaks = pd.concat([peaks, _df])

chroms.to_csv('../tests/test_unmix_chrom.csv', index=False)
peaks.to_csv('../tests/test_unmix_peaks.csv', index=False)
