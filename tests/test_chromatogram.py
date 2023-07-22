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
    chrom = hplc.quant.Chromatogram(test_data, bg_subtract=False, cols=colnames)
    peaks = chrom.detect_peaks(prominence=1E-3, verbose=False)

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
    chrom = hplc.quant.Chromatogram(data, cols={'time':'x', 'signal':'y'},
                                    peak_width=0.5)
    window = int(0.5 / np.mean(np.diff(data['x'].values)))
    assert np.isclose(chrom.df['estimated_background'].values[window:-window],
                      data['bg'].values[window:-window], rtol=tol).all()

test_bg_estimation()