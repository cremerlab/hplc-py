#%%
import hplc.quant
import pandas as pd
import numpy as np

def compare(a, b, tol):
    """
    Compares all elements in a and b and assigns equality within a tolerance, 
    accounting for annoying values near zero.
    """
    a = np.round(a, decimals=int(np.abs(np.round(np.log10(tol)))))
    b = np.round(b, decimals=int(np.abs(np.round(np.log10(tol))))) 
    _zeros = b == 0 
    b[_zeros] = np.sign(a[_zeros]) * tol
    assert np.isclose(a, b, rtol=tol).all()

def fit_peaks(test_data, truth, colnames={'time':'x', 'signal':'y'}, tol=1.5E-2):
    """
    Uses the `hplc.quant.Chromatogram.quantify` method to fit peaks in a chromatogram
    and compares the value with the ground truth.
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
    chrom_df = pd.read_csv('./tests/test_data/test_fitting_chrom.csv')
    peak_df = pd.read_csv('./tests/test_data/test_fitting_peaks.csv')
    for g, d in chrom_df.groupby('iter'):
        truth = peak_df[peak_df['iter']==g]
        fit_peaks(d, truth)

def test_peak_unmixing():
    """
    Tests that peaks can be properly unmixed and parameters estimates lie within 
    1% of the  ground truth.
    """
    # Load test data
    chrom_df = pd.read_csv('./tests/test_data/test_unmix_chrom.csv')
    peak_df = pd.read_csv('./tests/test_data/test_unmix_peaks.csv')
    for g, d in chrom_df.groupby('iter'):
        # Select true peak info
        truth = peak_df[peak_df['iter']==g]
        fit_peaks(d, truth)


def test_bg_estimation():
    """
    Tests that a background signal with a known profile can be correctly estimated 
    within 1.5% of the ground truth with a fixed window size.
    """
    tol = 1.5E-2
    data = pd.read_csv('./tests/test_data/test_SNIP_chrom.csv')
    chrom = hplc.quant.Chromatogram(data, cols={'time':'x', 'signal':'y'})
    chrom.correct_baseline(window=0.5)
    window = int(0.5 / np.mean(np.diff(data['x'].values)))
    assert np.isclose(chrom.df['estimated_background'].values[window:-window],
                      data['bg'].values[window:-window], rtol=tol).all()

def test_shouldered_peaks():
    """
    Tests that manually applied peaks can be properly deconvolved to within 1.5% 
    of the known parameter values.
    """
    tol = 1.5E-2
    data = pd.read_csv('./tests/test_data/test_manual_unmix_chrom.csv')
    peak_df = pd.read_csv('./tests/test_data/test_manual_unmix_peaks.csv')
    props = ['retention_time', 'amplitude', 'area', 'scale', 'skew']
    for g, d in data.groupby('iter'):
        truth = peak_df[peak_df['iter']==g]
        chrom = hplc.quant.Chromatogram(d, cols={'time':'x', 'signal':'y'})
        peaks = chrom.fit_peaks(enforced_locations=[11],# truth['retention_time'].values,
                                correct_baseline=False,
                                enforcement_tolerance=0.5)

        assert len(peaks) == len(truth)
        for p in props:
            compare(peaks[p].values, truth[p].values, tol)

def test_add_peak():
    """
    Tests that a peak that is not automatically detected that is not within 
    an extant peak window can be identified and deconvolved to within 1.5% of
    the known parameter values.    
    """
    data = pd.read_csv('./tests/test_data/test_shallow_signal_chrom.csv')
    props = ['retention_time', 'amplitude', 'area', 'scale', 'skew']
    peak_df = pd.read_csv('./tests/test_data/test_shallow_signal_peaks.csv')
    chrom = hplc.quant.Chromatogram(data, cols={'time':'x', 'signal':'y'})
    peaks = chrom.fit_peaks(enforced_locations=[25.0], correct_baseline=False, prominence=0.5)
    for p in props:
        compare(peaks[p].values, peak_df[p].values, 1.5E-2)

def test_score_reconstruction():
    """
    Tests that a known peak mixture is accurately reconstructed using R-scores 
    and Fano ratios. 
    """
    data = pd.read_csv('./tests/test_data/test_assessment_chrom.csv')
    peaks = pd.read_csv('./tests/test_data/test_assessment_peaks.csv')
    scores = pd.read_csv('./tests/test_data/test_assessment_scores.csv')
    loc = peaks[peaks['peak_idx']==2]['retention_time'].values[0]
    chrom = hplc.quant.Chromatogram(data, cols={'time':'x', 'signal':'y'})
    _ = chrom.fit_peaks(prominence=0.5)
    fit_scores =  chrom.assess_fit(verbose=False, tol=1E-3)
    for g, d in scores.groupby(['window_id', 'window_type']):
        _d = fit_scores[(fit_scores['window_id']==g[0]) & (fit_scores['window_type']==g[1])]['status'].values
        assert (_d == d['status'].values).all()
