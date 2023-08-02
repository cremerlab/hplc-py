# %%
import hplc.quant
import pandas as pd
import numpy as np
import pytest

def compare(a, b, tol):
    """
    Compares all elements in a and b and assigns equality within a tolerance, 
    accounting for annoying values near zero.
    """
    dec = int(np.abs(np.round(np.log10(tol))))
    a = np.round(a, decimals=dec)
    b = np.round(b, decimals=dec)
    _zeros = b == 0
    b[_zeros] = np.sign(a[_zeros]) * tol
    assert np.isclose(a, b, rtol=tol).all()


def fit_peaks(test_data, truth, colnames={'time': 'x', 'signal': 'y'}, tol=1.5E-2):
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
    # Make sure it fails if anything other than a dataframe is given. 
    try:
        chrom = hplc.quant.Chromatogram('test', cols={'time':'x', 'signal':'y'})
        assert False
    except RuntimeError:
        assert True

    # Load test data
    chrom_df = pd.read_csv('./tests/test_data/test_fitting_chrom.csv')
    chrom = hplc.quant.Chromatogram(chrom_df, cols={'time':'x', 'signal':'y'})


    try:
        chrom._assign_windows(prominence=-1)
    except ValueError:
        assert True
    try:
        chrom._assign_windows(prominence=2)
    except ValueError:
        assert True 
    try: 
        chrom._assign_windows(rel_height=-1)
    except ValueError:
        assert True
    try: 
        chrom._assign_windows(rel_height=2)
    except ValueError:
        assert True

    # Make sure a warning is thrown if a given buffer is < 10
    chrom_df = pd.read_csv('./tests/test_data/test_fitting_chrom.csv') 
    chrom = hplc.quant.Chromatogram(chrom_df[chrom_df['iter']==1], cols={'time':'x', 'signal':'y'})
    with pytest.warns():
        chrom.fit_peaks(buffer=9)

    peak_df = pd.read_csv('./tests/test_data/test_fitting_peaks.csv')
    for g, d in chrom_df.groupby('iter'):
        truth = peak_df[peak_df['iter'] == g]
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
        truth = peak_df[peak_df['iter'] == g]
        fit_peaks(d, truth)


def test_bg_estimation():
    """
    Tests that a background signal with a known profile can be correctly estimated 
    within 1.5% of the ground truth with a fixed window size.
    """
    tol = 1.5E-2
    data = pd.read_csv('./tests/test_data/test_SNIP_chrom.csv')
    chrom = hplc.quant.Chromatogram(data, cols={'time': 'x', 'signal': 'y'})
    chrom.correct_baseline(window=0.5)
    window = int(0.5 / np.mean(np.diff(data['x'].values)))
    assert np.isclose(chrom.df['estimated_background'].values[window:-window],
                      data['bg'].values[window:-window], rtol=tol).all()

    with pytest.warns():
        chrom.correct_baseline(window=0.5)

    data['y'] -= 100 
    chrom = hplc.quant.Chromatogram(data, cols={'time': 'x', 'signal': 'y'})
    with pytest.warns():
        chrom.correct_baseline(window=0.5)


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
        truth = peak_df[peak_df['iter'] == g]
        chrom = hplc.quant.Chromatogram(d, cols={'time': 'x', 'signal': 'y'})
        peaks = chrom.fit_peaks(enforced_locations=[11],
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
    chrom = hplc.quant.Chromatogram(data, cols={'time': 'x', 'signal': 'y'})
    peaks = chrom.fit_peaks(
        enforced_locations=[50.0], enforced_widths=[3], prominence=0.5, correct_baseline=False)
    for p in props:
        compare(peaks[p].values, peak_df[p].values, 1.5E-2)


def test_score_reconstruction():
    """
    Tests that a known peak mixture is accurately reconstructed using R-scores 
    and Fano ratios. 
    """
    data = pd.read_csv('./tests/test_data/test_assessment_chrom.csv')
    scores = pd.read_csv('./tests/test_data/test_assessment_scores.csv')
    chrom = hplc.quant.Chromatogram(data, cols={'time': 'x', 'signal': 'y'})
    _ = chrom.fit_peaks(prominence=0.9, rel_height=0.99)
    fit_scores = chrom.assess_fit(tol=1E-3, verbose=False)
    for g, d in scores.groupby(['window_id', 'window_type']):
        _d = fit_scores[(fit_scores['window_id'] == g[0]) & (
            fit_scores['window_type'] == g[1])]['status'].values
        assert (_d == d['status'].values).all()

def test_crop():
    """
    Tests that the crop function works as expected and throws exceptions when 
    improper time windows are given. 
    """
    data = pd.read_csv('./tests/test_data/test_assessment_chrom.csv')
    chrom = hplc.quant.Chromatogram(data, cols={'time':'x', 'signal':'y'})
    try:
        chrom.crop([1, 2, 3])
        assert False
    except ValueError:
        assert True
    try:
        chrom.crop([2, 1])
        assert False
    except RuntimeError:
        assert True
    chrom.crop([10, 20])
    assert (chrom.df.x.values[0] >= 10) & (chrom.df.x.values[-1] <= 20)
    _ = chrom.fit_peaks()
    try:
        chrom.crop([1, 2])
        assert False
    except RuntimeError:
        assert True

    # Test that cropping happens if a time window is provided.
    chrom = hplc.quant.Chromatogram(data, time_window=[10, 20], 
                                    cols={'time':'x', 'signal':'y'})
    assert (chrom.df.x.values[0] >= 10) & (chrom.df.x.values[-1] <= 20)

def test_deconvolve_peaks():
    """
    Tests that exception is properly thrown if peak fitting hasn't been performed.
    """
    data = pd.read_csv('./tests/test_data/test_assessment_chrom.csv')
    chrom = hplc.quant.Chromatogram(data, cols={'time':'x', 'signal':'y'})
    try:
        chrom.deconvolve_peaks()
    except RuntimeError:
        assert True

def test_map_peaks():
    """
    Tests that the peakmapping function correctly identifies peaks given retention
    times and tolerance and makes sure a linear calibration curve is used correctly. 
    """
    data = pd.read_csv('./tests/test_data/test_assessment_chrom.csv')
    chrom = hplc.quant.Chromatogram(data, cols={'time':'x', 'signal':'y'})

    # Check that peak mapping and calculation of the concentration works. 
    peaks = chrom.fit_peaks()
    orig_peaks = peaks.copy()
    params = {g:{'retention_time':d['retention_time'].values[0], 'intercept':0, 'slope':2} for g, d in orig_peaks.groupby('peak_id')}
    peaks = chrom.map_peaks(params)
    
    assert (orig_peaks['peak_id'].values == peaks['compound'].values).all() 
    assert (peaks['area'].values == 2 * peaks['concentration'].values).all() 

    # Check that mapping works if retention times are within tolerance
    chrom.fit_peaks() 
    params = {g:{'retention_time':d['retention_time'].values[0] + 0.1} for g, d in orig_peaks.groupby('peak_id')}
    peaks = chrom.map_peaks(params)
    assert (orig_peaks['peak_id'].values == peaks['compound'].values).all()

    # Check that it fails if no peaks can be found
    chrom.fit_peaks() 
    params = {g:{'retention_time':d['retention_time'].values[0] + 0.1} for g, d in orig_peaks.groupby('peak_id')}
    try:
        peaks = chrom.map_peaks(params, loc_tolerance=0.05)
        assert False
    except ValueError:
        assert True

    # Check that it fails if multiple peaks within the tolerance are found
    chrom.fit_peaks() 
    params = {g:{'retention_time':d['retention_time'].values[0]} for g, d in orig_peaks.groupby('peak_id')}
    try:
        peaks = chrom.map_peaks(params, loc_tolerance=5)
        assert False
    except ValueError:
        assert True
    with pytest.warns():
        params['test'] = {'retention_time': 1000}
        peaks = chrom.map_peaks(params)

def test_many_peaks():
    """
    Ensures that a warning is raised if there are 10 or more peaks in a given window. 
    """
    data = pd.read_csv('./tests/test_data/test_many_peaks.csv')
    chrom = hplc.quant.Chromatogram(data, cols={'time':'x', 'signal':'y'})
    with pytest.warns():  
        chrom.fit_peaks()
