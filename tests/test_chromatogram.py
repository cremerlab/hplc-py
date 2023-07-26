#%%
import hplc.quant
import pandas as pd
import numpy as np

def compare(a, b, tol):
    """
    Compares all elements in a and b and assigns equality within a tolerance, 
    accounting for annoying values near zero.
    """
    val = np.round(a, decimals=np.abs(int(np.log10(tol))))
    tru = np.round(b, decimals=np.abs(int(np.log10(tol))))
    _zeros = tru == 0 
    tru[_zeros] = np.sign(val[_zeros]) * tol 
    assert np.isclose(val, tru, rtol=tol).all()

def fit_peaks(test_data, truth, colnames={'time':'x', 'signal':'y'}, tol=1E-2):
    """
    Uses the `hplc.quant.Chromatogram.quantify` method to fit peaks in a chromatogram
    and compares teh value with the ground truth.
    """
    # Define constants
    props = ['retention_time', 'amplitude', 'area', 'scale', 'skew']

    # Execute analysis
    chrom = hplc.quant.Chromatogram(test_data, cols=colnames)
    peaks = chrom.fit_peaks(correct_baseline=False, prominence=1E-3)

    # Enforce that the correct number of peaks have been identified    
    assert len(peaks) == truth['peak_idx'].max()

    # Enforce that the inferred peak parameters are within a tolerance of 1%
    for p in props:  
        compare(peaks[p].values, truth[p].values, tol)

def test_peak_fitting():
    """
    Tests that peaks with known parameters can be faithfully estimated within 
    1% of the true value. If true parameter values are close to zero, victory is declared 
    if the estimated parameter is within 0.01.
    """
    # Load test data
    chrom_df = pd.read_csv('../tests/test_fitting_chrom.csv')
    peak_df = pd.read_csv('../tests/test_fitting_peaks.csv')
    for g, d in chrom_df.groupby('iter'):
        truth = peak_df[peak_df['iter']==g]
        fit_peaks(d, truth)

def test_peak_unmixing():
    """
    Tests that peaks can be properly unmixed and parameters estimates lie within 
    1% of the  ground truth.
    """
    # Load test data
    chrom_df = pd.read_csv('../tests/test_unmix_chrom.csv')
    peak_df = pd.read_csv('../tests/test_unmix_peaks.csv')
    for g, d in chrom_df.groupby('iter'):
        # Select true peak info
        truth = peak_df[peak_df['iter']==g]
        fit_peaks(d, truth)


def test_bg_estimation():
    """
    Tests that a background signal with a known profile can be correctly estimated 
    within 1% of the ground truth with a fixed window size.
    """
    tol = 1E-2
    data = pd.read_csv('./test_SNIP_chrom.csv')
    chrom = hplc.quant.Chromatogram(data, cols={'time':'x', 'signal':'y'})
    chrom.correct_baseline(window=0.5)
    window = int(0.5 / np.mean(np.diff(data['x'].values)))
    assert np.isclose(chrom.df['estimated_background'].values[window:-window],
                      data['bg'].values[window:-window], rtol=tol).all()

def test_add_peak():
    """
    Tests that manually applied peaks can be properly deconvolved to within 1% 
    of the known parameter values.
    """
    tol = 1E-2
    data = pd.read_csv('./test_manual_unmix_chrom.csv')
    peak_df = pd.read_csv('./test_manual_unmix_peaks.csv')
    props = ['retention_time', 'amplitude', 'area', 'scale', 'skew']
    for g, d in data.groupby('iter'):
        truth = peak_df[peak_df['iter']==g]
        chrom = hplc.quant.Chromatogram(d, cols={'time':'x', 'signal':'y'})
        peaks = chrom.fit_peaks(enforced_locations=truth['retention_time'].values,
                                correct_baseline=False,
                                # enforced_widths = truth['scale'].values/2,
                                enforcement_tolerance=0.25)
        print(peaks)
        # assert len(peaks) == len(truth)
        # for p in props:
            # compare(peaks[p].values, truth[p].values, tol)

# %%

# %%

# %%

import importlib
importlib.reload(hplc.quant)
# test_add_peak()
# %%

# %%

# %%

# %%

tol = 1E-2
data = pd.read_csv('./test_manual_unmix_chrom.csv')
peak_df = pd.read_csv('./test_manual_unmix_peaks.csv')
props = ['retention_time', 'amplitude', 'area', 'scale', 'skew']
for g, d in data.groupby('iter'):
    truth = peak_df[peak_df['iter']==g]
    chrom = hplc.quant.Chromatogram(d, cols={'time':'x', 'signal':'y'})
    peaks = chrom.fit_peaks(enforced_locations=truth['retention_time'].values,
                            correct_baseline=False,
                            enforced_widths = truth['scale'].values/2,
                            enforcement_tolerance=0.5)
    # print(peaks)
    break
    # assert len(

#%%

import matplotlib.pyplot as plt
chrom.show()

#%%
for p in props:
   compare(peaks[p].values, truth[p].values, 1E-2)
#%%
peaks[p].values
truth[p].values
# %%

# truth[p].values

peaks
#%%
peak_df