# %%
import hplc.quant
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pytest
import scipy.stats


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
    props = ['retention_time', 'amplitude',
             'area', 'scale', 'skew', 'signal_maximum']

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
    with pytest.raises(ValueError):
        chrom._assign_windows(rel_height=-1)
    with pytest.raises(ValueError):
        chrom._assign_windows(rel_height=2)

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
    _df = chrom.correct_baseline(window=0.75, return_df=True)
    # return (chrom, _df, data)
    with pytest.warns():
        __df = chrom.correct_baseline(window=0.75, return_df=False)

    # Ensure that dataframe returning is working.
    assert _df is not None
    assert __df is None

    window = int(0.75 / np.mean(np.diff(data['x'].values)))
    assert np.isclose(chrom.df['estimated_background'].values[window:-window],
                      data['bg'].values[window:-window], rtol=tol).all()

    with pytest.warns():
        chrom.correct_baseline(window=0.75)

    data['y'] -= 100
    chrom = hplc.quant.Chromatogram(data, cols={'time': 'x', 'signal': 'y'})
    with pytest.warns():
        chrom.correct_baseline(window=0.75)

    # Test that the correct exception is thrown if the window is too small.
    try:
        with pytest.warns():
            chrom.correct_baseline(window=np.diff(_df['x'].values)[0])
            assert False
    except ValueError:
        assert True


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
    with pytest.raises(ValueError):
        chrom.crop([1, 2, 3])
    with pytest.raises(RuntimeError):
        chrom.crop([2, 1])
    # A missing/None time window should give a clear ValueError, not a TypeError.
    with pytest.raises(ValueError, match='must be provided as a list'):
        chrom.crop(None)
    with pytest.raises(ValueError, match='must be provided as a list'):
        chrom.crop()

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
    with pytest.raises(RuntimeError):
        chrom.deconvolve_peaks()


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


def test_map_peaks_missing_compound_not_last():
    """
    Regression test for GitHub issue #23: a compound with no matching peak must
    not abort mapping of the remaining (present) compounds. Previously the loop
    `break`ed on the first unmapped compound, so any compound listed after it was
    silently skipped, and a leading miss raised "No peaks could be properly
    mapped!" even though every other peak was present.
    """
    data = pd.read_csv('./tests/test_data/test_assessment_chrom.csv')
    chrom = hplc.quant.Chromatogram(data, cols={'time': 'x', 'signal': 'y'})
    peaks = chrom.fit_peaks()

    # Build a params dict whose FIRST entry has no matching peak, followed by the
    # real, present peaks. Insertion order is preserved by Python dicts.
    params = {'absent': {'retention_time': peaks['retention_time'].min() - 5}}
    for g, d in peaks.groupby('peak_id'):
        params[f'compound_{g}'] = {'retention_time': d['retention_time'].values[0]}

    with pytest.warns(match='No peak found for absent'):
        mapped = chrom.map_peaks(params)

    # All present compounds are mapped despite the leading miss.
    assert len(mapped) == len(peaks)
    assert set(mapped['compound'].values) == {
        f'compound_{g}' for g in peaks['peak_id'].values}
    assert 'absent' not in mapped['compound'].values


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

    # Make sure that incorrect bounds throw a value error.
    adjustments = {'scale': [10, 11]}
    try:
        _ = chrom.fit_peaks(param_bounds=adjustments,
                            correct_baseline=False,
                            verbose=False)
        assert False
    except ValueError:
        assert True


def _skewnorm_signal(t, params):
    """Build a signal as a sum of amplitude-weighted (skew)normal peaks."""
    sig = np.zeros_like(t, dtype=float)
    for amp, loc, scale, alpha in params:
        sig += amp * scipy.stats.skewnorm(alpha, loc, scale).pdf(t)
    return sig


def test_known_peaks_nonzero_start_time():
    """
    Regression test for bug E: enforced (`known_peaks`) locations are mapped to
    array indices relative to the chromatogram's start time, not by assuming the
    time axis begins at t = 0. On a chromatogram whose time starts at t != 0, a
    shallow peak that is not auto-detected must still be enforced at the correct
    retention time.
    """
    t = np.arange(10, 30, 0.01)
    # A tall auto-detected peak and a shallow peak that prominence filtering skips.
    sig = _skewnorm_signal(t, [(1000, 16.0, 0.3, 0), (80, 22.0, 0.3, 0)])
    df = pd.DataFrame({'time': t, 'signal': sig})
    chrom = hplc.quant.Chromatogram(df)
    peaks = chrom.fit_peaks(known_peaks={22.0: {'width': 1}},
                            prominence=0.5, correct_baseline=False, verbose=False)
    # The enforced peak is placed at ~22 (would land out-of-range on the t0=0 bug).
    assert np.any(np.abs(peaks['retention_time'].values - 22.0) < 0.5)
    assert np.any(np.abs(peaks['retention_time'].values - 16.0) < 0.5)


def test_unmixed_columns_match_peak_id():
    """
    Regression test for bug F: `unmixed_chromatograms[:, peak_id - 1]` must hold
    the trace for that `peak_id`. When detection order differs from
    retention-time order (here forced by enforcing an early peak in a window that
    also contains a later auto-detected peak), the columns were previously left
    in detection order and disagreed with the sorted `peak_id`.
    """
    t = np.arange(0, 30, 0.01)
    # Auto-detected tall peak at 16.0; shallow enforced peak earlier at 15.0.
    sig = _skewnorm_signal(t, [(1000, 16.0, 0.2, 0), (90, 15.0, 0.2, 0)])
    df = pd.DataFrame({'time': t, 'signal': sig})
    chrom = hplc.quant.Chromatogram(df)
    peaks = chrom.fit_peaks(known_peaks={15.0: {'width': 1}}, prominence=0.5,
                            correct_baseline=False, buffer=200, verbose=False)

    time = chrom.df['time'].values
    for _, row in peaks.iterrows():
        col = chrom.unmixed_chromatograms[:, int(row['peak_id']) - 1]
        t_at_max = time[np.argmax(col)]
        assert np.abs(t_at_max - row['retention_time']) < 0.5, (
            f"column for peak_id {row['peak_id']} peaks at {t_at_max}, "
            f"expected near {row['retention_time']}")


def test_infeasible_location_guess_does_not_crash():
    """
    Regression test for bug H: the location initial guess is the peak time
    *rounded* to `_time_precision`, while its bounds are the window's *raw*
    (unrounded) time range. When `_time_precision = |ceil(log10(dt))|` is coarser
    than the actual sample spacing, rounding can push the guess just past the
    window's max, which scipy rejects as "`x0` is infeasible". The clamp must pull
    the guess back inside its bounds so fitting proceeds.

    This is purely a function of `dt` vs `_time_precision`, so it is reproduced
    synthetically rather than from recorded data: a step of `dt = 0.011` gives
    `_time_precision = |ceil(log10(0.011))| = 1`, so a peak near the end of the
    record (here ~2.985) rounds *up* to 3.0 — past its window's raw max (~2.98).
    On the unfixed code this raises before any optimizer iteration; with the clamp
    the fit completes. (Verified: removing the clamp makes this signal raise
    "Initial guess for 'location' ... lies outside its bounds".)
    """
    t = np.arange(0, 3.0, 0.011)
    sig = _skewnorm_signal(t, [(6000, 1.0, 0.03, 0), (6000, 2.985, 0.03, 0)])
    df = pd.DataFrame({'time': t, 'signal': sig})
    assert int(np.abs(np.ceil(np.log10(np.mean(np.diff(t)))))) == 1, (
        'fixture must keep `_time_precision == 1` for the rounding excursion')

    chrom = hplc.quant.Chromatogram(df)
    try:
        peaks = chrom.fit_peaks(correct_baseline=False, verbose=False,
                                max_iter=5000)
    except ValueError as e:
        # The clamp removes the bounds/feasibility crash entirely; if a ValueError
        # still surfaces it must not be the infeasible-bounds one under test.
        assert 'infeasible' not in str(e).lower()
        assert 'outside its bounds' not in str(e).lower()
        assert 'lower bound' not in str(e).lower()
        raise
    # The fit proceeds past the bounds check and recovers both peaks.
    assert len(peaks) == 2


def test_peak_adjacent_to_start_assigns_interpeak():
    """
    Smoke test for bug G: a chromatogram with a peak right at the start and a
    long trailing background must run end-to-end and assign interpeak windows
    over the background. The G fix corrects background-window splitting when a
    gap lands at the first sample (`split_inds[0] == 0`). That exact branch is
    hard to trigger deterministically without depending on `scipy.peak_widths`
    internals, so this guards the surrounding behavior rather than isolating the
    branch.
    """
    t = np.arange(0, 30, 0.01)
    sig = _skewnorm_signal(t, [(400, 0.4, 0.15, 0), (600, 18.0, 0.4, 0)])
    df = pd.DataFrame({'time': t, 'signal': sig})
    chrom = hplc.quant.Chromatogram(df)
    chrom.fit_peaks(correct_baseline=False, verbose=False, max_iter=5000)
    scores = chrom.assess_fit(verbose=False)
    # The trailing background between the two peaks is captured as interpeak.
    assert (chrom.window_df['window_type'] == 'interpeak').any()
    assert (scores['window_type'] == 'interpeak').any()


def test_degenerate_bounds_raise_clear_error():
    """
    The pre-fit bounds validation (the safety net added for bug H) must surface
    degenerate or infeasible parameter bounds as a clear, actionable ValueError
    instead of scipy's opaque message. Issue #22 is exactly the degenerate
    amplitude-bound case.
    """
    t = np.arange(0, 20, 0.01)
    sig = _skewnorm_signal(t, [(1000, 10.0, 0.3, 0)])
    df = pd.DataFrame({'time': t, 'signal': sig})

    # Degenerate amplitude bound [0, 0]: lower bound not strictly less than upper.
    chrom = hplc.quant.Chromatogram(df)
    with pytest.raises(ValueError, match='Invalid bounds'):
        chrom.fit_peaks(param_bounds={'amplitude': [0, 0]},
                        correct_baseline=False, verbose=False)

    # Amplitude bound [2x, 3x] the peak value excludes the initial guess (1x).
    chrom = hplc.quant.Chromatogram(df)
    with pytest.raises(ValueError, match='lies outside its bounds'):
        chrom.fit_peaks(param_bounds={'amplitude': [2, 3]},
                        correct_baseline=False, verbose=False)


def test_multipeak_param_bounds_validated_per_peak():
    """
    Regression test for bug B: when a global `param_bounds` is applied to a
    window containing more than one peak, the initial-guess-vs-bounds check must
    use *each* peak's own guess. Previously it always inspected the first peak's
    guess, so a bound that excluded a later peak's guess was silently accepted
    (and later surfaced as an opaque scipy error rather than the informative one).
    """
    # Two overlapping peaks with clearly different widths that share one window.
    t = np.arange(0, 30, 0.01)
    sig = (200 * scipy.stats.norm(14.0, 0.2).pdf(t)
           + 200 * scipy.stats.norm(15.0, 0.6).pdf(t))
    df = pd.DataFrame({'time': t, 'signal': sig})

    # Inspect the per-peak scale initial guesses (scale guess = width / 2).
    probe = hplc.quant.Chromatogram(df)
    probe._assign_windows(buffer=200)
    multi = [v for v in probe.window_props.values() if v['num_peaks'] == 2]
    assert len(multi) == 1, "expected the two peaks to share a single window"
    guesses = sorted(w / 2 for w in multi[0]['width'])
    # Bound contains the smaller guess (first peak) but excludes the larger one.
    bound = [0, 0.5 * (guesses[0] + guesses[1])]

    chrom = hplc.quant.Chromatogram(df)
    with pytest.raises(ValueError, match='exclusive of initial guess'):
        chrom.fit_peaks(correct_baseline=False, buffer=200,
                        param_bounds={'scale': bound}, verbose=False)


def test_narrow_peak_scale_guess_clamped():
    """
    Regression test: a very narrow peak (e.g. a noise spike passing a low
    `prominence`) collapses its fitting window to roughly its own footprint, so
    the scale initial guess (half-width) can slightly exceed the upper scale
    bound (half the window's time range). Like the location guess, the scale
    guess must be clamped into its bounds so the window is deconvolved instead of
    raising "Initial guess for 'scale' ... lies outside its bounds".
    """
    t = np.arange(60) * 0.1
    sig = np.full(60, 0.01)
    sig[30] = 1.0  # single-sample spike -> window collapses to the peak footprint
    chrom = hplc.quant.Chromatogram(pd.DataFrame({'time': t, 'signal': sig}))
    # Must not raise, and must recover the single peak near t = 3.0.
    peaks = chrom.fit_peaks(prominence=0.05, buffer=0,
                            correct_baseline=False, verbose=False)
    assert np.isclose(peaks['retention_time'].values, 3.0, atol=0.3).any()
