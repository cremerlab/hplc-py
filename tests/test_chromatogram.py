# %%

import importlib
import hplc.quant
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
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
    assert chrom._fitting_progress_state == 1

    # Ensure that proper representation is applied.
    assert 'Peak(s) Detected' in chrom.__repr__()
    assert 'Baseline Subtracted' not in chrom.__repr__()
    assert 'Enforced Peak Location(s)' not in chrom.__repr__()
    assert 'Compoun(s) Assigned' not in chrom.__repr__()
    assert 'Cropped' not in chrom.__repr__()

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
        chrom = hplc.quant.Chromatogram(
            'test', cols={'time': 'x', 'signal': 'y'})
        assert False
    except RuntimeError:
        assert True

    # Load test data
    chrom_df = pd.read_csv('./tests/test_data/test_fitting_chrom.csv')
    chrom = hplc.quant.Chromatogram(
        chrom_df, cols={'time': 'x', 'signal': 'y'})
    try:
        chrom._assign_windows(rel_height=-1)
    except ValueError:
        assert True
    try:
        chrom._assign_windows(rel_height=2)
    except ValueError:
        assert True

    chrom_df = pd.read_csv('./tests/test_data/test_fitting_chrom.csv')
    chrom = hplc.quant.Chromatogram(chrom_df[chrom_df['iter'] == 1], cols={
                                    'time': 'x', 'signal': 'y'})
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
    _df = chrom.correct_baseline(window=0.5, return_df=True)
    with pytest.warns():
        __df = chrom.correct_baseline(window=0.5, return_df=False)

    # Ensure that dataframe returning is working.
    assert _df is not None
    assert __df is None

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
        peaks = chrom.fit_peaks(known_peaks=[11],
                                correct_baseline=False,
                                tolerance=0.5)

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
        known_peaks={50.0: {'width': 3}}, prominence=0.5, correct_baseline=False)
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

    # Test that proper RuntimeException is thrown if no reconstruction is present.
    try:
        chrom.assess_fit()
        assert False
    except RuntimeError:
        assert True

    _ = chrom.fit_peaks(prominence=0.9, rel_height=0.99, buffer=100)
    fit_scores = chrom.assess_fit(rtol=1E-3, verbose=False)
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
    chrom = hplc.quant.Chromatogram(data, cols={'time': 'x', 'signal': 'y'})
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

    # Test that a dataframe is returned only if specified.
    no_returned_df = chrom.crop([10, 20], return_df=False)
    returned_df = chrom.crop([10, 20], return_df=True)
    assert no_returned_df is None
    assert returned_df is not None

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
                                    cols={'time': 'x', 'signal': 'y'})
    assert 'Cropped' in chrom.__repr__()
    assert 'Baseline Subtracted' not in chrom.__repr__()
    assert 'Peak(s) Detected' not in chrom.__repr__()
    assert 'Enforced Peak Location(s)' not in chrom.__repr__()
    assert (chrom.df.x.values[0] >= 10) & (chrom.df.x.values[-1] <= 20)


def test_deconvolve_peaks():
    """
    Tests that exception is properly thrown if peak fitting hasn't been performed.
    """
    data = pd.read_csv('./tests/test_data/test_assessment_chrom.csv')
    chrom = hplc.quant.Chromatogram(data, cols={'time': 'x', 'signal': 'y'})
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
    chrom = hplc.quant.Chromatogram(data, cols={'time': 'x', 'signal': 'y'})

    # Check that peak mapping and calculation of the concentration works.
    peaks = chrom.fit_peaks()
    orig_peaks = peaks.copy()
    params = {g: {'retention_time': d['retention_time'].values[0],
                  'intercept': 0, 'slope': 2} for g, d in orig_peaks.groupby('peak_id')}
    peaks = chrom.map_peaks(params)
    assert 'Compound(s) Assigned' in chrom.__repr__()
    assert 'Baseline Subtracted' in chrom.__repr__()
    assert 'Cropped' not in chrom.__repr__()
    assert 'Enforced Peak Location(s)' not in chrom.__repr__()
    assert (orig_peaks['peak_id'].values == peaks['compound'].values).all()
    assert (peaks['area'].values == 2 * peaks['concentration'].values).all()

    # Check that mapping works if retention times are within tolerance
    chrom.fit_peaks()
    params = {g: {'retention_time': d['retention_time'].values[0] + 0.1}
              for g, d in orig_peaks.groupby('peak_id')}
    peaks = chrom.map_peaks(params)
    assert (orig_peaks['peak_id'].values == peaks['compound'].values).all()

    # Check that it fails if no peaks can be found
    chrom.fit_peaks()
    params = {g: {'retention_time': d['retention_time'].values[0] + 0.1}
              for g, d in orig_peaks.groupby('peak_id')}
    try:
        peaks = chrom.map_peaks(params, loc_tolerance=0.05)
        assert False
    except ValueError:
        assert True

    # Check that it fails if multiple peaks within the tolerance are found
    chrom.fit_peaks()
    params = {g: {'retention_time': d['retention_time'].values[0]}
              for g, d in orig_peaks.groupby('peak_id')}
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
    chrom = hplc.quant.Chromatogram(data, cols={'time': 'x', 'signal': 'y'})
    with pytest.warns():
        chrom.fit_peaks()


def test_bounding():
    """
    Ensures that custom bounding of parameters can resolve heavily overlapping 
    peaks within a tolerance of 5%. This higher tolerance is due to the 
    difficulty of resolving heavily overlapping peaks. 
    """
    tol = 5E-2
    data = pd.read_csv('./tests/test_data/test_bounding_chroms.csv')
    true_peaks = pd.read_csv('./tests/test_data/test_bounding_peaks.csv')
    bounding_factors = np.array([0.9, 1.1])
    for g, d in data.groupby('iter'):
        truth = true_peaks[true_peaks['iter'] == g].copy()
        truth.sort_values('peak_id', inplace=True)
        peak2 = truth[truth['peak_id'] == 2]
        peak1 = truth[truth['peak_id'] == 1]
        bounds = {peak2['retention_time'].values[0]: {
            'amplitude': peak2['amplitude'].values[0] * bounding_factors,
            'scale': peak2['scale'].values[0] * bounding_factors,
            'skew': peak2['skew'].values[0] * bounding_factors,
            'location': peak2['retention_time'].values[0] * bounding_factors},
            peak1['retention_time'].values[0]: {}}

        # Assert that it fails without providing locations.
        chrom = hplc.quant.Chromatogram(d,  cols={'time': 'x', 'signal': 'y'})
        bad_peaks = chrom.fit_peaks(correct_baseline=False)
        assert len(truth) != len(bad_peaks)

        # Fit with provided bounding
        peaks = chrom.fit_peaks(
            known_peaks=bounds, correct_baseline=False, buffer=100)
        assert 'Peak(s) Detected' in chrom.__repr__()
        assert 'Enforced Peak Location(s)' in chrom.__repr__()
        assert 'Baseline Subtracted' not in chrom.__repr__()
        assert 'Compound(s) Assigned' not in chrom.__repr__()
        assert len(truth) == len(peaks)

        # Ensure that it's close to within tolerance
        peaks.sort_values('peak_id', inplace=True)
        props = ['retention_time', 'amplitude', 'area', 'scale', 'skew']
        for p in props:
            compare(peaks[p].values, truth[p].values, tol)


def test_variable_integration_area():
    """
    Tests that the integration window is adjusted correctly and measured peak 
    area agrees with the ground truth to within a tolerance of 1.5%
    """
    df = pd.read_csv('./tests/test_data/test_integration_window_chrom.csv')
    chrom = hplc.quant.Chromatogram(df)

    # Ensure that the test fails if a nonsense integration window is supplied.
    win = [1]
    try:
        _ = chrom.fit_peaks(integration_window=win)
        assert False
    except RuntimeError:
        assert True

    # Load the window area dataframe
    areas = pd.read_csv('./tests/test_data/test_integration_window_areas.csv')
    for g, d in areas.groupby(['t_start', 't_end', 'window']):
        _area = d['area'].values[0]
        if g[2] == 1:
            peaks = chrom.fit_peaks()
        else:
            win = [g[0], g[1]]
            peaks = chrom.fit_peaks(integration_window=win)
        assert np.isclose(peaks['area'].values[0], _area, rtol=1.5E-2)


def test_verbosity():
    """
    Ensures that verbosity flags are respective 
    """
    df = pd.read_csv('./tests/test_data/test_integration_window_chrom.csv')

    # Peak Fitting and Baseline Subtraction
    chrom = hplc.quant.Chromatogram(df)
    _ = chrom.fit_peaks(verbose=True)
    assert chrom._fitting_progress_state == 1
    assert chrom._bg_correction_progress_state == 1

    chrom = hplc.quant.Chromatogram(df)
    _ = chrom.fit_peaks(verbose=False)
    assert chrom._fitting_progress_state == 0
    assert chrom._bg_correction_progress_state == 0

    # Reconstruction Reporting
    _ = chrom.assess_fit(verbose=False)
    assert chrom._report_card_progress_state == 0
    _ = chrom.assess_fit(verbose=True)
    assert chrom._report_card_progress_state == 1


def test_show():
    """
    Ensures that chromatogram visualization is showing features as expected. 
    """
    df = pd.read_csv('./tests/test_data/test_integration_window_chrom.csv')
    chrom = hplc.quant.Chromatogram(df)
    _ = chrom.show()
    plt.close()
    assert chrom._viz_ylabel_subtraction_indication == False
    assert chrom._viz_subtracted_baseline == False
    assert chrom._viz_peak_reconstruction == False
    assert chrom._viz_adjusted_xlim == False

    _ = chrom.show(time_range=[2, 3])
    plt.close()
    assert chrom._viz_ylabel_subtraction_indication == False
    assert chrom._viz_subtracted_baseline == False
    assert chrom._viz_peak_reconstruction == False
    assert chrom._viz_adjusted_xlim == True

    chrom.correct_baseline()
    _ = chrom.show()
    plt.close()
    assert chrom._viz_ylabel_subtraction_indication == True
    assert chrom._viz_peak_reconstruction == False
    assert chrom._viz_subtracted_baseline == True
    assert chrom._viz_adjusted_xlim == False

    _ = chrom.fit_peaks(verbose=False)
    _ = chrom.show()
    plt.close()
    assert chrom._viz_peak_reconstruction == True
    assert chrom._viz_ylabel_subtraction_indication == True
    assert chrom._viz_subtracted_baseline == True
    assert chrom._viz_adjusted_xlim == False
    assert chrom._viz_mapped_peaks == False

    _ = chrom.map_peaks({'test': {'retention_time': 10}},)
    _ = chrom.show()
    plt.close()
    assert chrom._viz_peak_reconstruction == True
    assert chrom._viz_ylabel_subtraction_indication == True
    assert chrom._viz_subtracted_baseline == True
    assert chrom._viz_adjusted_xlim == False
    assert chrom._viz_mapped_peaks == True
    assert chrom._viz_min_one_concentration == False

    _ = chrom.map_peaks(
        {'test': {'retention_time': 10, 'slope': 1, 'intercept': 1}})
    _ = chrom.show()
    plt.close()
    assert chrom._viz_peak_reconstruction == True
    assert chrom._viz_ylabel_subtraction_indication == True
    assert chrom._viz_subtracted_baseline == True
    assert chrom._viz_adjusted_xlim == False
    assert chrom._viz_mapped_peaks == True
    assert chrom._viz_min_one_concentration == True
    assert chrom._viz_unit_display == False

    _ = chrom.map_peaks(
        {'test': {'retention_time': 10, 'slope': 1, 'intercept': 1, 'unit': 'test'}})
    _ = chrom.show()
    plt.close()
    assert chrom._viz_peak_reconstruction == True
    assert chrom._viz_ylabel_subtraction_indication == True
    assert chrom._viz_subtracted_baseline == True
    assert chrom._viz_adjusted_xlim == False
    assert chrom._viz_mapped_peaks == True
    assert chrom._viz_min_one_concentration == True
    assert chrom._viz_unit_display == True


def test_generic_param_bounding():
    """
    Tests that global parameter bounds can be maniupulated. 
    """
    df = pd.read_csv('./tests/test_data/test_integration_window_chrom.csv')
    chrom = hplc.quant.Chromatogram(df)
    # Adjust the parameters
    adjustments = {'amplitude': [0.9, 1.1],
                   'location': [-1, 1],
                   'scale': [1, 3],
                   'skew': [-10, 10]}
    _ = chrom.fit_peaks(param_bounds=adjustments, verbose=False)
    adj_pars = chrom._param_bounds[0]

    # Make sure the adjustments match
    _loc = chrom.df['time'].values[chrom._peak_indices]
    _amp = chrom.df['signal_corrected'].values[chrom._peak_indices]
    truth = {'amplitude': _amp * adjustments['amplitude'],
             'location': _loc + adjustments['location'],
             'scale': [1, 3],
             'skew': [-10, 10]}
    for k, v in adj_pars.items():
        assert np.array(truth[k] == v).all()
