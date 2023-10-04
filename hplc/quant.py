import pandas as pd
import numpy as np
import scipy.signal
import scipy.optimize
import scipy.special
import tqdm
import matplotlib.pyplot as plt
import warnings
import seaborn as sns
import termcolor


class Chromatogram:
    """
    Base class for the processing and quantification of an HPLC chromatogram.

    Attributes
    ----------
    df : `pandas.core.frame.DataFrame`    
        A Pandas DataFrame containing the chromatogram, minimally with columns 
        of time and signal intensity. 
    window_props : `dict`
       A dictionary of each peak window, labeled as increasing integers in 
       linear order. Each key has its own dictionary with the following keys:
    peaks : `pandas.core.frame.DataFrame` 
        A Pandas DataFrame containing the inferred properties of each peak 
        including the retention time, scale, skew, amplitude, and total
        area under the peak across the entire chromatogram.
    unmixed_chromatograms : `numpy.ndarray`
        A matrix where each row corresponds to a time point and each column corresponds
        to the value of the probability density for each individual peak. This 
        is used primarily for plotting in the `show` method. 
    quantified_peaks : `pandas.core.frame.DataFrame`
        A Pandas DataFrame with peak areas converted to 
    scores : `pandas.core.frame.DataFrame`
        A Pandas DataFrame containing the reconstruction scores and Fano factor 
        ratios for each peak and interpeak region. This is generated only afer 
        `assess_fit()` is called.
    """

    def __init__(self,
                 file,
                 time_window=None,
                 cols={'time': 'time',
                       'signal': 'signal'}):
        """
        Instantiates a chromatogram object on which peak detection and quantification
        is performed.

        Parameters
        ----------
        file: `str` or pandas.core.frame.DataFrame`
            The path to the csv file of the chromatogram to analyze or 
            the pandas DataFrame of the chromatogram. If None, a pandas DataFrame 
            of the chromatogram must be passed.
        dataframe : `pandas.core.frame.DataFrame`
            a Pandas DataFrame of the chromatogram to analyze. If None, 
            a path to the csv file must be passed
        time_window: `list` [start, end], optional
            The retention time window of the chromatogram to consider for analysis.
            If None, the entire time range of the chromatogram will be considered.
        cols: `dict`, keys of 'time', and 'signal', optional
            A dictionary of the retention time and intensity measurements 
            of the chromatogram. Default is 
            `{'time':'time', 'signal':'signal'}`. 
       """

        # Peform type checks and throw exceptions where necessary.
        if (type(file) is not pd.core.frame.DataFrame):
            raise RuntimeError(
                f'Argument must be a Pandas DataFrame. Argument is of type {type(file)}')

        # Assign class variables
        self.time_col = cols['time']
        self.int_col = cols['signal']

        # Load the chromatogram and necessary components to self.
        dataframe = file.copy()
        self.df = dataframe

        # Define the average timestep in the chromatogram. This computes a mean
        # but values will typically be identical.
        self._dt = np.mean(np.diff(dataframe[self.time_col].values))
        self._time_precision = int(np.abs(np.ceil(np.log10(self._dt))))

       # Blank out vars that are used elsewhere
        self.window_props = None
        self.scores = None
        self._peak_indices = None
        self.peaks = None
        self._guesses = None
        self._bg_corrected = False
        self._mapped_peaks = None
        self._added_peaks = None
        self.unmixed_chromatograms = None
        self._crop_offset = 0

        # Prune to time window
        if time_window is not None:
            self.crop(time_window)
        else:
            self.df = dataframe

    def __repr__(self):
        trange = f'(t: {self.df[self.time_col].values[0]} - {self.df[self.time_col].values[-1]})'
        rep = f"""Chromatogram:"""
        if self._crop_offset > 0:
            rep += f'\n\t\u2713 Cropped {trange}'
        if self._bg_corrected:
            rep += f'\n\t\u2713 Baseline Subtracted'
        if self.peaks is not None:
            rep += f"\n\t\u2713 {self.peaks.peak_id.max()} Peak(s) Detected"
            if self._added_peaks is not None:
                rep += f'\n\t\u2713 Enforced Peak Location(s)'
        if self._mapped_peaks is not None:
            rep += f'\n\t\u2713 Compound(s) Assigned'
        return rep

    def crop(self,
             time_window=None,
             return_df=False):
        R"""
        Restricts the time dimension of the DataFrame in place.

        Parameters
        ----------
        time_window : `list` [start, end], optional
            The retention time window of the chromatogram to consider for analysis.
            If None, the entire time range of the chromatogram will be considered.
        return_df : `bool`
            If `True`, the cropped DataFrame is 

        Returns
        -------
        cropped_df : pandas DataFrame
            If `return_df = True`, then the cropped dataframe is returned.
        """
        if self.peaks is not None:
            raise RuntimeError("""
You are trying to crop a chromatogram after it has been fit. Make sure that you 
do this before calling `fit_peaks()` or provide the argument `time_window` to the `fit_peaks()`.""")
        if len(time_window) != 2:
            raise ValueError(
                f'`time_window` must be of length 2 (corresponding to start and end points). Provided list is of length {len(time_window)}.')
        if time_window[0] >= time_window[1]:
            raise RuntimeError(
                f'First index in `time_window` must be ≤ second index.')

        # Apply the crop and return
        self.df = self.df[(self.df[self.time_col] >= time_window[0]) &
                          (self.df[self.time_col] <= time_window[1])]
        # set the offset
        self._crop_offset = int(time_window[0] / self._dt)
        if return_df:
            return self.df

    def _assign_windows(self,
                        known_peaks=[],
                        tolerance=0.5,
                        prominence=0.01,
                        rel_height=1,
                        buffer=0,
                        peak_kwargs={}):
        R"""
        Breaks the provided chromatogram down to windows of likely peaks. 

        Parameters
        ----------
        known_peaks : `list`
            The approximate locations of the peaks. If this is not provided, 
            peak locations will be automatically detected. 
        tolerance: `float`, optional
            If an enforced peak location is within tolerance of an automatically 
            identified peak, the automatically identified peak will be preferred. 
            This parameter is in units of time. Default is one-half time unit.
        prominence : `float`
            The promimence threshold for identifying peaks. Prominence is the 
            relative height of signal relative to the local background. 
        rel_height : `float`, [0, 1]
            The relative height of the peak where the baseline is determined. 
            Default is 1.
        buffer : positive `int`
            The padding of peak windows in units of number of time steps. Default 
            is 100 points on each side of the identified peak window.

        Returns
        ------- 
        window_df : `pandas.core.frame.DataFrame`
            A Pandas DataFrame with each measurement assigned to an identified 
            peak or overlapping peak set. This returns a copy of the chromatogram
            DataFrame with  a column  for the local baseline and one column for 
            the window IDs. Window ID of -1 corresponds to area not assigned to 
            any peaks
        """
        if (rel_height < 0) | (rel_height > 1):
            raise ValueError(f' `rel_height` must be [0, 1].')

        # Correct for a negative baseline
        df = self.df
        intensity = self.df[self.int_col].values
        int_sign = np.sign(intensity)
        norm_int = (intensity - intensity.min()) / \
            (intensity.max() - intensity.min())
        self.normint = int_sign * norm_int

        # Preform automated peak detection and set window ranges
        peaks, _ = scipy.signal.find_peaks(
            int_sign * norm_int, prominence=prominence, **peak_kwargs)
        self._peak_indices = peaks

        # Get the amplitudes and the indices of each peak
        amps = np.sign(intensity[peaks])
        pos_inds = amps > 0
        neg_inds = amps < 0

        # Set up storage vectors for peak quantities
        _widths = np.zeros_like(amps)
        _left = np.zeros(len(amps)).astype(int)
        _right = np.zeros(len(amps)).astype(int)

        # Get the peak properties
        with warnings.catch_warnings():
            warnings.simplefilter(
                "ignore", category=scipy.signal._peak_finding_utils.PeakPropertyWarning)
            for i, inds in enumerate([pos_inds, neg_inds]):
                if i == 0:
                    _intensity = intensity
                else:
                    _intensity = -intensity
                if len(inds) > 0:
                    __widths, _, _, _ = scipy.signal.peak_widths(_intensity,
                                                                 self._peak_indices,
                                                                 rel_height=0.5)
                    _widths[inds] = __widths[inds]

                    _, _, __left, __right = scipy.signal.peak_widths(_intensity,
                                                                     self._peak_indices,
                                                                     rel_height=rel_height)
                    _left[inds] = __left[inds]
                    _right[inds] = __right[inds]

        # Determine if peaks should be added.
        if len(known_peaks) > 0:
            # Get the enforced peak positions
            if type(known_peaks) == dict:
                _known_peaks = list(known_peaks.keys())
            else:
                _known_peaks = known_peaks

            # Find the nearest location in the time array given the user-specified time
            enforced_location_inds = np.int_(
                np.array(_known_peaks) / self._dt) - self._crop_offset

            # Update the user specified times with the nearest location
            updated_loc = np.round(
                self._dt * (enforced_location_inds + self._crop_offset), decimals=self._time_precision)
            if type(known_peaks) == dict:
                updated_known_peaks = known_peaks.copy()
                for _new, _old in zip(updated_loc, _known_peaks):
                    updated_known_peaks[_new] = updated_known_peaks.pop(_old)
                _known_peaks = list(updated_known_peaks)
            else:
                updated_known_peaks = list(np.copy(known_peaks))
                for _old, _new in enumerate(updated_loc):
                    updated_known_peaks[_old] = _new
                _known_peaks = updated_known_peaks
            self._known_peaks = updated_known_peaks

            # Clear peacks that are within a tolerance of the provided locations
            # of known peaks.
            for i, loc in enumerate(enforced_location_inds):
                autodetected = np.nonzero(
                    np.abs(self._peak_indices - loc) <= (tolerance / self._dt))[0]
                if len(autodetected) > 0:
                    # Remove the autodetected peak
                    self._peak_indices = np.delete(
                        self._peak_indices, autodetected[0])
                    _widths = np.delete(_widths, autodetected[0])
                    _left = np.delete(_left, autodetected[0])
                    _right = np.delete(_right, autodetected[0])

            # Add the provided locations of known peaks and adjust parameters as necessary.
            for i, loc in enumerate(enforced_location_inds):
                self._peak_indices = np.append(self._peak_indices, loc)
                if self._added_peaks is None:
                    self._added_peaks = []
                self._added_peaks.append((loc + self._crop_offset) * self._dt)
                if type(known_peaks) == dict:
                    _sel_loc = updated_known_peaks[_known_peaks[i]]
                    if 'width' in _sel_loc.keys():
                        _widths = np.append(
                            _widths, _sel_loc['width'] / self._dt)
                    else:
                        _widths = np.append(_widths, 1 / self._dt)
                else:
                    _widths = np.append(_widths, 1 / self._dt)
                _left = np.append(_left, loc - _widths[-1] - buffer)
                _right = np.append(_right, loc + _widths[-1] + buffer)
        else:
            self._known_peaks = known_peaks

        # Set window ranges
        ranges = []
        for l, r in zip(_left, _right):
            _range = np.arange(int(l - buffer), int(r + buffer), 1)
            _range = _range[(_range >= 0) & (_range <= len(norm_int))]
            ranges.append(_range)

        # Identiy subset ranges and remove
        valid = [True] * len(ranges)
        if len(ranges) > 1:
            for i, r1 in enumerate(ranges):
                for j, r2 in enumerate(ranges):
                    if i != j:
                        if set(r2).issubset(r1):
                            valid[j] = False

        # Keep only valid ranges and baselines
        ranges = [r for i, r in enumerate(ranges) if valid[i] is True]

        # Copy the dataframe and return the windows
        window_df = df.copy(deep=True)
        window_df.sort_values(by=self.time_col, inplace=True)
        window_df['time_idx'] = np.arange(len(window_df))
        window_df['window_id'] = 0
        window_df['window_type'] = 'peak'
        self.ranges = ranges
        for i, r in enumerate(ranges):
            window_df.loc[window_df['time_idx'].isin(r),
                          'window_id'] = int(i + 1)

        # Determine the windows for the background (nonpeak) areas.
        bg_windows = window_df[window_df['window_id'] == 0]
        tidx = bg_windows['time_idx'].values

        if len(bg_windows) > 0:
            split_inds = np.nonzero(
                np.diff(bg_windows['time_idx'].values) - 1)[0]

            # If there is only one background window
            if (len(split_inds) == 0):
                window_df.loc[window_df['time_idx'].isin(
                    bg_windows['time_idx'].values), 'window_id'] = 1
                window_df.loc[window_df['time_idx'].isin(
                    bg_windows['time_idx'].values), 'window_type'] = 'interpeak'

            # If more than one split ind, set up all ranges.
            elif split_inds[0] != 0:
                split_inds += 1
                split_inds = np.insert(split_inds, 0, 0)
                split_inds = np.append(split_inds, len(tidx))

            bg_ranges = [bg_windows.iloc[np.arange(
                split_inds[i], split_inds[i+1], 1)]['time_idx'].values for i in range(len(split_inds)-1)]
            win_idx = 1
            for i, rng in enumerate(bg_ranges):
                if len(rng) >= 10:
                    window_df.loc[window_df['time_idx'].isin(
                        rng), 'window_id'] = win_idx
                    window_df.loc[window_df['time_idx'].isin(
                        rng), 'window_type'] = 'interpeak'
                    win_idx += 1
        window_df = window_df[window_df['window_id'] > 0]

        # Convert this to a dictionary for easy parsing
        window_dict = {}
        for g, d in window_df[window_df['window_type'] == 'peak'].groupby('window_id'):
            if g > 0:
                _peaks = [
                    p for p in self._peak_indices if p in d['time_idx'].values]
                peak_inds = [x for _p in _peaks for x in np.where(
                    self._peak_indices == _p)[0]]
                _dict = {'time_range': d[self.time_col].values,
                         'signal': d[self.int_col].values,
                         'signal_area': d[self.int_col].values.sum(),
                         'num_peaks': len(_peaks),
                         'amplitude': [d[d['time_idx'] == p][self.int_col].values[0] for p in _peaks],
                         'location': [np.round(d[d['time_idx'] == p][self.time_col].values[0], decimals=self._time_precision) for p in _peaks],
                         'width':  [_widths[ind] * self._dt for ind in peak_inds]}
                window_dict[int(g)] = _dict

        self.window_df = window_df
        self.window_props = window_dict
        return window_df

    def _compute_skewnorm(self,
                          x,
                          *params):
        R"""
        Computes the lineshape of a skew-normal distribution given the shape,
        location, and scale parameters

        Parameters
        ----------
        x : `float` or `numpy.ndarray`
            The time dimension of the skewnorm 
        params : `list`, [`amplitude`, `loc`, `scale`, `alpha`]
            Parameters for the shape and scale parameters of the skewnorm 
            distribution.
                `amplitude` : positive `float`
                    Height of the peak.
                `loc` : positive `float`
                    The location parameter of the distribution.
                `scale` : positive `float`
                    The scale parameter of the distribution.
                `alpha` : positive `float`
                    The skew shape parameter of the distribution.

        Returns
        -------
        scaled_pdf : `float or numpy array, same shape as `x`
            The PDF of the skew-normal distribution scaled with the supplied 
            amplitude.

        Notes
        -----
        This function infers the parameters defining skew-normal distributions 
        for each peak in the chromatogram. The fitted distribution has the form 

        .. math:: 
            I = 2S_\text{max} \left(\frac{1}{\sqrt{2\pi\sigma^2}}\right)e^{-\frac{(t - r_t)^2}{2\sigma^2}}\left[1 + \text{erf}\frac{\alpha(t - r_t)}{\sqrt{2\sigma^2}}\right]

        where :math:`S_\text{max}` is the maximum signal of the peak, 
        :math:`t` is the time, :math:`r_t` is the retention time, :math:`\sigma`
        is the scale parameter, and :math:`\alpha` is the skew parameter.

        """
        amp, loc, scale, alpha = params
        _x = alpha * (x - loc) / scale
        norm = np.sqrt(2 * np.pi * scale**2)**-1 * \
            np.exp(-(x - loc)**2 / (2 * scale**2))
        cdf = 0.5 * (1 + scipy.special.erf(_x / np.sqrt(2)))
        return amp * 2 * norm * cdf

    def _fit_skewnorms(self,
                       x,
                       *params):
        R"""
        Estimates the parameters of the distributions which consititute the 
        peaks in the chromatogram. 

        Parameters
        ----------
        x : `float`
            The time dimension of the skewnorm 
        params : list of length 4 x number of peaks, [amplitude, loc, scale, alpha]
            Parameters for the shape and scale parameters of the skewnorm 
            distribution. Must be provided in following order, repeating
            for each distribution.
                `amplitude` : float; > 0
                    Height of the peak.
                `loc` : float; > 0
                    The location parameter of the distribution.
                `scale` : float; > 0
                    The scale parameter of the distribution.
                `alpha` : float; > 
                    The skew shape parater of the distribution.

        Returns
        -------
        out : `float`
            The evaluated distribution at the given time x. This is the summed
            value for all distributions modeled to construct the peak in the 
            chromatogram.
        """
        # Get the number of peaks and reshape for easy indexing
        n_peaks = int(len(params) / 4)
        params = np.reshape(params, (n_peaks, 4))
        out = 0

        # Evaluate each distribution
        for i in range(n_peaks):
            out += self._compute_skewnorm(x, *params[i])
        return out

    def deconvolve_peaks(self,
                         verbose=True,
                         known_peaks=[],
                         param_bounds={},
                         integration_window=[],
                         max_iter=1000000,
                         optimizer_kwargs={}):
        R"""
        .. note::
           In most cases, this function should not be called directly. Instead, 
           it should called through the :func:`~hplc.quant.Chromatogram.fit_peaks`

        For each peak window, estimate the parameters of skew-normal distributions 
        which makeup the peak(s) in the window. See "Notes" for information on
        default parameter bounds.

        Parameters
        ----------
        verbose : `bool`
            If `True`, a progress bar will be printed during the inference.

        param_bounds : `dict`
            Modifications to the default parameter bounds (see Notes below) as 
            a dictionary for each parameter. A dict entry should be of the 
            form `parameter: [lower, upper]`. Modifications have the following effects:
                * Modifications to `amplitude` bounds are multiplicative of the 
                  observed magnitude at the peak position. 
                * Modifications to `location` are values that are subtracted or 
                  added from the peak position for lower and upper bounds, respectively.
                * Modifications to `scale` replace the default values. 
                * Modifications to `skew` replace the default values. 
        integration_window: `list`
            The time window over which the integrated peak areas should be computed. 
            If empty, the area will be integrated over the entire duration of the 
            cropped chromatogram.
        max_iter : `int`
            The maximum number of iterations the optimization protocol should 
            take before erroring out. Default value is 10^6.

        optimizer_kwargs : `dict`
            Keyword arguments to be passed to `scipy.optimize.curve_fit`.

        Returns 
        --------
        peak_props: `dict`
            A dataframe containing properties of the peak fitting procedure. 

        Notes
        -----
        The parameter boundaries are set automatically to prevent run-away estimation 
        into non-realistic regimes that can seriously slow down the inference. The 
        default parameter boundaries for each peak are as follows.

            * `amplitude`: The lower and upper peak amplitude boundaries correspond to one-tenth and ten-times the value of the peak at the peak location in the chromatogram.

            * `location`: The lower and upper location bounds correspond to the minimum and maximum time values of the chromatogram.

            * `scale`: The lower and upper bounds of the peak standard deviation defaults to the chromatogram time-step and one-half of the chromatogram duration, respectively.  

            * `skew`: The skew parameter by default is allowed to take any value between (-`inf`, `inf`).
        """
        if self.window_props is None:
            raise RuntimeError(
                'Function `_assign_windows` must be run first. Go do that.')
        if verbose:
            self._fitting_progress_state = 1
            iterator = tqdm.tqdm(self.window_props.items(),
                                 desc='Deconvolving mixture')
        else:
            self._fitting_progress_state = 0
            iterator = self.window_props.items()
        parorder = ['amplitude', 'location', 'scale', 'skew']
        paridx = {k: -(i+1) for i, k in enumerate(reversed(parorder))}
        peak_props = {}
        self._bounds = []

        # Determine the areas over which to integrate the window
        if len(integration_window) == 0:
            t_range = self.df[self.time_col].values
        elif (type(integration_window) == list):
            if len(integration_window) == 2:
                t_range = np.arange(
                    integration_window[0], integration_window[1], self._dt)
            else:
                raise RuntimeError(
                    'Provided integration bounds has wrong dimensions. Should have a length of 2.')

        # Instantiate a state variable to ensure that parameter bounds are being adjusted.
        self._param_bounds = []
        for k, v in iterator:
            window_dict = {}

            # Set up the initial guess
            p0 = []
            bounds = [[],  []]

            # If there are more than 5 peaks in a mixture, throw a warning
            if v['num_peaks'] >= 10:
                warnings.warn(f"""
-------------------------- Hey! Yo! Heads up! ----------------------------------
| This time window (from {np.round(v['time_range'].min(), decimals=4)} to {np.round(v['time_range'].max(), decimals=3)}) has {v['num_peaks']} candidate peaks.
| This is a complex mixture and may take a long time to properly fit depending 
| on how well resolved the peaks are. Reduce `buffer` if the peaks in this      
| window should be separable by eye. Or maybe just go get something to drink.
--------------------------------------------------------------------------------
""")
            if v['num_peaks'] == 0:
                continue
            for i in range(v['num_peaks']):
                p0.append(v['amplitude'][i])
                p0.append(v['location'][i]),
                p0.append(v['width'][i] / 2)  # scale parameter
                p0.append(0)  # Skew parameter, starts with assuming Gaussian

                # Set default parameter bounds
                _param_bounds = {'amplitude': np.sort([0.1 * v['amplitude'][i], 10 * v['amplitude'][i]]),
                                 'location': [v['time_range'].min(), v['time_range'].max()],
                                 'scale': [self._dt, (v['time_range'].max() - v['time_range'].min())/2],
                                 'skew': [-np.inf, np.inf]}

                # Modify the parameter bounds given arguments
                if len(param_bounds) != 0:
                    for p in parorder:
                        if p in param_bounds.keys():
                            if p == 'amplitude':
                                _param_bounds[p] = v['amplitude'][i] * \
                                    np.sort(param_bounds[p])
                            elif p == 'location':
                                _param_bounds[p] = [v['location']
                                                    [i] + p for p in param_bounds[p]]
                            else:
                                _param_bounds[p] = param_bounds[p]

                # Add peak-specific bounds if provided
                if (type(known_peaks) == dict) & (len(known_peaks) != 0):
                    if v['location'][i] in known_peaks.keys():
                        newbounds = known_peaks[v['location'][i]]
                        tweaked = False
                        if len(newbounds) > 0:
                            for p in parorder:
                                if p in newbounds.keys():
                                    _param_bounds[p] = np.sort(newbounds[p])
                                    tweaked = True
                                    if p != 'location':
                                        p0[paridx[p]] = np.mean(newbounds[p])

                           # Check if width is the only key
                            if (len(newbounds) >= 1) & ('width' not in newbounds.keys()):
                                if tweaked == False:
                                    raise ValueError(
                                        f"Could not adjust bounds for peak at {v['location'][i]} because bound keys do not contain at least one of the following: `location`, `amplitude`, `scale`, `skew`. ")

                for _, val in _param_bounds.items():
                    bounds[0].append(val[0])
                    bounds[1].append(val[1])
            self._param_bounds.append(_param_bounds)

            # Perform the inference
            popt, _ = scipy.optimize.curve_fit(self._fit_skewnorms, v['time_range'],
                                               v['signal'], p0=p0, bounds=bounds, maxfev=max_iter,
                                               **optimizer_kwargs)

            # Assemble the dictionary of output
            if v['num_peaks'] > 1:
                popt = np.reshape(popt, (v['num_peaks'], 4))
            else:
                popt = [popt]
            for i, p in enumerate(popt):
                window_dict[f'peak_{i + 1}'] = {
                    'amplitude': p[0],
                    'retention_time': np.round(p[1], decimals=self._time_precision),
                    'scale': p[2],
                    'alpha': p[3],
                    'area': self._compute_skewnorm(t_range, *p).sum(),
                    'reconstructed_signal': self._compute_skewnorm(v['time_range'], *p)}
            peak_props[k] = window_dict

        self._peak_props = peak_props
        return peak_props

    def fit_peaks(self,
                  known_peaks=[],
                  tolerance=0.5,
                  prominence=1E-2,
                  rel_height=1,
                  approx_peak_width=5,
                  buffer=0,
                  param_bounds={},
                  integration_window=[],
                  verbose=True,
                  return_peaks=True,
                  correct_baseline=True,
                  max_iter=1000000,
                  precision=9,
                  peak_kwargs={},
                  optimizer_kwargs={}):
        R"""
        Detects and fits peaks present in the chromatogram

        Parameters
        ----------
        known_peaks : `list` or `dict`
            The approximate locations of peaks whose position is known. If 
            provided as a list, only the locations wil be used as initial guesses. 
            If provided as a dictionary, locations and parameter bounds will be 
            set.
        tolerance: `float`, optional
            If an enforced peak location is within tolerance of an automatically 
            identified peak, the automatically identified peak will be preferred. 
            This parameter is in units of time. Default is one-half time unit.
        prominence : `float`,  [0, 1]
            The promimence threshold for identifying peaks. Prominence is the 
            relative height of the normalized signal relative to the local
            background. Default is 1%. If `locations` is provided, this is 
            not used.
        rel_height : `float`, [0, 1]
            The relative height of the peak where the baseline is determined. This
            is used to split into windows and is *not* used for peak detection.
            Default is 100%. 
        approx_peak_width: `float`, optional
            The approximate width of the signal you want to quantify. This is 
            used as filtering window for automatic baseline correction. If `correct_baseline==False`,
            this has no effect. 
        buffer : positive `int`
            The padding of peak windows in units of number of time steps. Default 
            is 100 points on each side of the identified peak window. Must have a value 
            of at least 10.
        verbose : `bool`
            If True, a progress bar will be printed during the inference. 
        param_bounds: `dict`, optional
            Parameter boundary modifications to be used to constrain fitting of 
            all peaks. 
            See docstring of :func:`~hplc.quant.Chromatogram.deconvolve_peaks`
            for more information.
        integration_window: `list`
            The time window over which the integrated peak areas should be computed. 
            If empty, the area will be integrated over the entire duration of the 
            cropped chromatogram.
        correct_baseline : `bool`, optional
            If True, the baseline of the chromatogram will be automatically 
            corrected using the SNIP algorithm. See :func:`~hplc.quant.Chromatogram.correct_baseline`
            for more information.
        return_peaks : `bool`, optional
            If True, a dataframe containing the peaks will be returned. Default
            is True.
        max_iter : `int`
            The maximum number of iterations the optimization protocol should 
            take before erroring out. Default value is 10^6.
        precision : `int`
            The number of decimals to round the reconstructed signal to. Default
            is 9.
        peak_kwargs : `dict`
            Additional arguments to be passed to `scipy.signal.find_peaks`.
        optimizer_kwargs : `dict`
            Additional arguments to be passed to `scipy.optimize.curve_fit`.

        Returns
        -------
        peak_df : `pandas.core.frame.DataFrame`
            A dataframe containing information for each detected peak. This is
            only returned if `return_peaks == True`. The peaks are always 
            stored as an attribute `peak_df`.

        Notes
        -----
        This function infers the parameters defining skew-norma distributions 
        for each peak in the chromatogram. The fitted distribution has the form 

        .. math:: 
            I = 2S_\text{max} \left(\frac{1}{\sqrt{2\pi\sigma^2}}\right)e^{-\frac{(t - r_t)^2}{2\sigma^2}}\left[1 + \text{erf}\frac{\alpha(t - r_t)}{\sqrt{2\sigma^2}}\right]

        where :math:`S_\text{max}` is the maximum signal of the peak, 
        :math:`t` is the time, :math:`r_t` is the retention time, :math:`\sigma`
        is the scale parameter, and :math:`\alpha` is the skew parameter.

        """
        if correct_baseline and not self._bg_corrected:
            self.correct_baseline(window=approx_peak_width,
                                  verbose=verbose, return_df=False)

        # Assign the window bounds
        _ = self._assign_windows(known_peaks=known_peaks,
                                 tolerance=tolerance,
                                 prominence=prominence, rel_height=rel_height,
                                 buffer=buffer, peak_kwargs=peak_kwargs)

        # Infer the distributions for the peaks
        peak_props = self.deconvolve_peaks(verbose=verbose,
                                           known_peaks=self._known_peaks,
                                           param_bounds=param_bounds,
                                           max_iter=max_iter,
                                           integration_window=integration_window,
                                           **optimizer_kwargs)

        # Set up a dataframe of the peak properties
        peak_df = pd.DataFrame([])
        iter = 0
        for _, peaks in peak_props.items():
            for _, params in peaks.items():
                _dict = {'retention_time': params['retention_time'],
                         'scale': params['scale'],
                         'skew': params['alpha'],
                         'amplitude': params['amplitude'],
                         'area': params['area']}
                iter += 1
                peak_df = pd.concat([peak_df, pd.DataFrame(_dict, index=[0])])

        peak_df.sort_values(by='retention_time', inplace=True)
        peak_df['peak_id'] = np.arange(len(peak_df)) + 1
        peak_df['peak_id'] = peak_df['peak_id'].astype(int)
        self.peaks = peak_df

        # Compute the mixture
        time = self.df[self.time_col].values
        out = np.zeros((len(time), len(peak_df)))
        iter = 0
        for _, _v in self._peak_props.items():
            for _, v in _v.items():
                params = [v['amplitude'], v['retention_time'],
                          v['scale'], v['alpha']]
                out[:, iter] = self._compute_skewnorm(time, *params)
                iter += 1
        self.unmixed_chromatograms = np.round(out, decimals=precision)
        if return_peaks:
            return peak_df

    def correct_baseline(self,
                         window=5,
                         return_df=False,
                         verbose=True,
                         precision=9):
        R"""
        Performs Sensitive Nonlinear Iterative Peak (SNIP) clipping to estimate 
        and subtract background in chromatogram.

        Parameters
        ----------
        window : `int`
            The approximate size of signal objects in the chromatogram in dimensions
            of time. This is related to the number of iterations undertaken by 
            the SNIP algorithm.
        return_df : `bool`
            If `True`, then chromatograms (before and after background correction) are returned
        verbose: `bool`
            If `True`, progress will be printed to screen as a progress bar. 
        precision: `int`
            The number of decimals to round the subtracted signal to. Default is 9.

        Returns
        -------
        corrected_df : `pandas.core.frame.DataFrame`
            If `return_df = True`, then the original and the corrected chromatogram are returned.

        Notes
        -----
        This implements the SNIP algorithm as presented and summarized in `Morhác
        and Matousek 2008 <https://doi.org/10.1366/000370208783412762>`_. The 
        implementation here also rounds to 9 decimal places in the subtracted signal
        to avoid small values very near zero.
        """
        if self._bg_corrected == True:
            warnings.warn(
                'Baseline has already been corrected. Rerunning on original signal...')
            self.int_col = self.int_col.split('_corrected')[0]

        # Unpack and copy dataframe and intensity profile
        df = self.df
        signal = df[self.int_col].copy()

        # Look at the relative magnitudes of the maximum and minimum values
        # And throw a warning if there are appreciable negative peaks.
        min_val = np.min(signal)
        max_val = np.max(signal)
        if min_val < 0:
            if (np.abs(min_val) / max_val) >= 0.1:
                warnings.warn("""
\x1b[30m\x1b[43m\x1b[1m
The chromatogram appears to have appreciable negative signal . Automated background 
subtraction may not work as expected. Proceed with caution and visually 
check if the subtraction is acceptable!
\x1b[0m""")

        # Clip the signal if the median value is negative
        if (signal < 0).any():
            shift = np.median(signal[signal < 0])
        else:
            shift = 0
        signal -= shift
        signal *= np.heaviside(signal, 0)

        # Compute the LLS operator
        tform = np.log(np.log(np.sqrt(signal.values + 1) + 1) + 1)

        # Compute the number of iterations given the window size.
        n_iter = int(((window / self._dt) - 1) / 2)

        # Iteratively filter the signal
        if verbose:
            self._bg_correction_progress_state = 1
            iter = tqdm.tqdm(range(1, n_iter + 1),
                             desc="Performing baseline correction")
        else:
            self._bg_correction_progress_state = 0
            iter = range(1, n_iter + 1)

        for i in iter:
            tform_new = tform.copy()
            for j in range(i, len(tform) - i):
                tform_new[j] = min(tform_new[j], 0.5 *
                                   (tform_new[j+i] + tform_new[j-i]))
            tform = tform_new

        # Perform the inverse of the LLS transformation and subtract
        inv_tform = ((np.exp(np.exp(tform) - 1) - 1)**2 - 1)
        df[f'{self.int_col}_corrected'] = np.round(
            (df[self.int_col].values - shift - inv_tform), decimals=precision)
        df[f'estimated_background'] = inv_tform + shift
        self.df = df
        self._bg_corrected = True
        self.int_col = f'{self.int_col}_corrected'
        if return_df:
            return df

    def map_peaks(self,
                  params,
                  loc_tolerance=0.5,
                  include_unmapped=False):
        R"""
        Maps user-provided mappings to arbitrarily labeled peaks. If a linear 
        calibration curve is also provided, the concentration will be computed.

        Parameters
        ----------
        params : `dict` 
            A dictionary mapping each peak to a slope and intercept used for 
            converting peak areas to units of concentraions. Each peak should
            have a key that is the compound name (e.g. "glucose"). Each key
            should have another dict as the key with `retention_time` , `slope` ,
            and `intercept` as keys. If only `retention_time` is given,
            concentration will not be computed. The key `retention_time` will be
            used to map the compound to the `peak_id`. If `unit` are provided,
            this will be added as a column
       loc_tolerance : `float`
           The tolerance for mapping the compounds to the retention time. The 
           default is 0.5 time units.
       include_unmapped : `bool`
            If True, unmapped compounds will remain in the returned peak dataframe,
            but will be populated with Nan. Default is False.

       Returns
       -------
       peaks : `pandas.core.frame.DataFrame`
            A modified peak table with the compound name and concentration 
            added as columns.

        Notes
        -----
        .. note::
            As of `v0.1.0`, this function can only accommodate linear calibration 
            functions.
        """
        # Create a mapper for peak id to compound
        mapper = {}
        unmapped = {}
        peak_df = self.peaks.copy()
        for k, v in params.items():
            ret_time = v['retention_time']
            peak_id = np.abs(
                peak_df['retention_time'].values - ret_time) < loc_tolerance

            if np.sum(peak_id) > 1:
                raise ValueError(
                    f"Multiple compounds found within tolerance of retention time for {k}. Reduce the tolerance or correct the provided value.")

            if np.sum(peak_id) == 0:
                unmapped[k] = v['retention_time']
                break
            peak_id = peak_df.peak_id.values[np.argmax(peak_id)]
            peak_df.loc[peak_df['peak_id'] == peak_id, 'compound'] = k
            mapper[peak_id] = k
        if len(mapper) == 0:
            raise ValueError(
                "No peaks could be properly mapped! Check your provided retention times.")
        if len(unmapped) > 0:
            for k, v in unmapped.items():
                warnings.warn(
                    f"\nNo peak found for {k} (retention time {v}) within the provided tolerance.")

        # Iterate through the compounds and calculate the concentration.
        for g, d in peak_df.groupby('compound'):
            if (g in params.keys()):
                if 'slope' in params[g].keys():
                    conc = (d['area'] - params[g]['intercept']) / \
                        params[g]['slope']
                    peak_df.loc[peak_df['compound']
                                == g, 'concentration'] = conc
                    if 'unit' in params[g].keys():
                        peak_df.loc[peak_df['compound'] ==
                                    g, 'unit'] = params[g]['unit']
        if include_unmapped == False:
            peak_df.dropna(inplace=True)
        self.quantified_peaks = peak_df
        self._mapped_peaks = mapper
        return peak_df

    def _score_reconstruction(self):
        R"""
        Computes the reconstruction score on a per-window and total chromatogram
        basis. 

        Returns
        -------
        score_df : `pandas.core.frame.DataFrame`
            A DataFrame reporting the scoring statistic for each window as well 
            as for the entire chromatogram. A window value of `0` corresponds 
            to the chromatogram regions which don't have peaks. A window 
            value of `-1` corresponds to the chromatogram as a whole

        Notes
        -----
        The reconstruction score is defined as
        ..math:: 

            R = \frac{\text{area of inferred mixture in window} + 1}{\text{area of observed signal in window} + 1} = \frac{\frac{\sum\limits_{i\in t}^t \sum\limits_{j \in N_\text{peaks}}^{N_\text{peaks}}2A_j \text{SkewNormal}(\alpha_j, r_{t_j}, \sigma_j) + 1}{\sum\limits_{i \in t}^t S_i + 1}

        where :math:`t` is the total time of the region, :math:`A`is the inferred 
        peak amplitude, :math:`\alpha` is the inferred skew paramter, :math:`r_t` is
        the inferred peak retention time, :math:`\sigma` is the inferred scale 
        parameter and :math:`S_i` is the observed signal intensity at time point
        :math:`i`. Note that the signal and reconstruction is cast to be positive
        to compute the score.  

        """
        columns = ['window_id', 'time_start', 'time_end', 'signal_area',
                   'inferred_area', 'signal_variance', 'signal_mean', 'signal_fano_factor', 'reconstruction_score']
        score_df = pd.DataFrame([])
        # Compute the per-window reconstruction

        for g, d in self.window_df[self.window_df['window_type'] == 'peak'].groupby('window_id'):
            # Compute the non-peak windows separately.
            window_area = np.abs(d[self.int_col].values).sum() + 1
            window_peaks = self._peak_props[g]
            window_peak_area = np.array(
                [np.abs(v['reconstructed_signal']) for v in window_peaks.values()]).sum() + 1
            score = np.array(window_peak_area / window_area).astype(float)
            signal_var = np.var(np.abs(d[self.int_col].values))
            signal_mean = np.mean(np.abs(d[self.int_col].values))
            # Account for an edge case to avoid dividing by zero
            if signal_mean == 0:
                signal_mean += 1E-9
            signal_fano = signal_var / signal_mean
            x = [g, d[self.time_col].min(),
                 d[self.time_col].max(), window_area, window_peak_area,
                 signal_var, signal_mean, signal_fano, score]
            _df = pd.DataFrame(
                {_c: _x for _c, _x in zip(columns, x)}, index=[g - 1])
            _df['window_type'] = 'peak'
            score_df = pd.concat([score_df, _df])

        # Compute the score for the non-peak regions
        nonpeak = self.window_df[self.window_df['window_type'] == 'interpeak']
        if len(nonpeak) > 0:
            for g, d in nonpeak.groupby('window_id'):
                total_area = np.abs(d[self.int_col].values).sum() + 1
                recon_area = np.sum(np.abs(self.unmixed_chromatograms), axis=1)[
                    d['time_idx'].values].sum() + 1
                nonpeak_score = recon_area / total_area
                signal_var = np.var(np.abs(d[self.int_col].values))
                signal_mean = np.mean(np.abs(d[self.int_col].values))
                # Account for an edge case to avoide dividing by zero
                if signal_mean == 0:
                    signal_mean += 1E-9
                signal_fano = signal_var / signal_mean

                # Add to score dataframe
                x = [g, d[self.time_col].min(),
                     d[self.time_col].max(),
                     total_area, recon_area, signal_var, signal_mean, signal_fano, nonpeak_score]
                _df = pd.DataFrame(
                    {c: xi for c, xi in zip(columns, x)}, index=[g - 1])
                _df['window_type'] = 'interpeak'
                score_df = pd.concat([score_df, _df])
        score_df['signal_area'] = score_df['signal_area'].astype(float)
        score_df['inferred_area'] = score_df['inferred_area'].astype(float)
        self.scores = score_df
        return score_df

    def assess_fit(self,
                   rtol=1E-2,
                   fano_tol=1E-2,
                   verbose=True):
        R"""
        Assesses whether the computed reconstruction score is adequate, given a tolerance.

        Parameters
        ----------
        rtol : `float`
            The tolerance for a reconstruction to be valid. This is the tolerated 
            deviation from a score of 1 which indicates a perfectly reconstructed
            chromatogram. 
        fano_tol : `float`
            The tolerance away from zero for evaluating the Fano factor of 
            inerpeak windows. See note below.
        verbose : `bool`
            If True, a summary of the fit will be printed to screen indicating 
            problematic regions if detected.

        Returns
        -------
        score_df : `pandas.core.frame.DataFrame`  
            A DataFrame reporting the scoring statistic for each window as well 
            as for the entire chromatogram. A window value of `0` corresponds 
            to the entire chromatogram. A column `accepted` with a boolean 
            value represents whether the reconstruction is within tolerance (`True`)
            or (`False`).

        Notes
        -----
        The reconstruction score is defined as

        .. math:: 
            R = \frac{\text{area of inferred mixture in window} + 1}{\text{area of observed signal in window} + 1} 

        where :math:`t` is the total time of the region, :math:`A` is the inferred 
        peak amplitude, :math:`\alpha` is the inferred skew paramter, :math:`r_t` is
        the inferred peak retention time, :math:`\sigma` is the inferred scale 
        parameter and :math:`S_i` is the observed signal intensity at time point
        :math:`i`. Note that the signal and reconstruction is cast to be positive
        to compute the score.  

        A reconstruction score of :math:`R = 1` indicates a perfect 
        reconstruction of the chromatogram. For practical purposes, a chromatogram
        is deemed to be adequately reconstructed if :math:`R` is within a tolerance
        :math:`\epsilon` of 1 such that

        .. math::
            \left| R - 1 \right| \leq \epsilon \Rightarrow \text{Valid Reconstruction}

        Interpeak regions may have a poor reconstruction score due to noise or
        short durations. To determine if this poor reconstruction score is due 
        to a missed peak, the signal Fano factor of the region is computed as 

        .. math::
            F = \frac{\sigma^2_{S}}{\langle S \rangle}.

        This is compared with the average Fano factor of :math:`N` peak windows such 
        that the Fano factor ratio is 

        .. math::
            \frac{F}{\langle F_{peak} \rangle} = \frac{\sigma^2_{S} / \langle S \rangle}{\frac{1}{N} \sum\limits_{i}^N \frac{\sigma_{S,i}^2}{\langle S_i \rangle}}.

        If the Fano factor ratio is below a tolerance `fano_tol`, then that 
        window is deemed to be noisy and peak-free. 
        """

        if self.unmixed_chromatograms is None:
            raise RuntimeError(
                "No reconstruction found! `.fit_peaks()` must be called first. Go do that.")

        # Compute the reconstruction score
        _score_df = self._score_reconstruction()

        # Apply the tolerance parameter
        _score_df['applied_tolerance'] = rtol
        score_df = pd.DataFrame([])
        mean_fano = _score_df[_score_df['window_type']
                              == 'peak']['signal_fano_factor'].mean()
        for g, d in _score_df.groupby(['window_type', 'window_id']):

            tolpass = np.round(np.abs(d['reconstruction_score'].values[0] - 1),
                               decimals=int(np.abs(np.ceil(np.log10(rtol))))) <= rtol

            d = d.copy()
            if g[0] == 'peak':
                if tolpass:
                    d['status'] = 'valid'
                else:
                    d['status'] = 'invalid'

            else:
                fanopass = (d['signal_fano_factor'].values[0] /
                            mean_fano) <= fano_tol
                if tolpass:
                    d['status'] = 'valid'
                elif fanopass:
                    d['status'] = 'needs review'
                else:
                    d['status'] = 'invalid'
            score_df = pd.concat([score_df, d], sort=False)

        # Define colors printing parameters to avoid retyping everything.
        print_colors = {'valid': ('A+, Success: ', ('black', 'on_green')),
                        'invalid': ('F, Failed: ', ('black', 'on_red')),
                        'needs review': ('C-, Needs Review: ', ('black', 'on_yellow'))}
        if verbose:
            self._report_card_progress_state = 1
            print("""
-------------------Chromatogram Reconstruction Report Card----------------------

Reconstruction of Peaks
======================= 
""")
        else:
            self._report_card_progress_state = 0

        for g, d in score_df[score_df['window_type'] == 'peak'].groupby('window_id'):
            status = d['status'].values[0]
            if status == 'valid':
                warning = ''
            else:
                warning = """
Peak mixture poorly reconstructs signal. You many need to adjust parameter bounds 
or add manual peak positions (if you have a shouldered pair, for example). If 
you have a very noisy signal, you may need to increase the reconstruction 
tolerance `rtol`."""
            if (d['reconstruction_score'].values[0] >= 10) | \
                    (d['reconstruction_score'].values[0] <= 0.1):
                if d['reconstruction_score'].values[0] == 0:
                    r_score = f'0'
                else:
                    r_score = f"10^{int(np.log10(d['reconstruction_score'].values[0]))}"
            else:
                r_score = f"{d['reconstruction_score'].values[0]:0.4f}"
            if verbose:
                termcolor.cprint(f"{print_colors[status][0]} Peak Window {int(g)} (t: {d['time_start'].values[0]:0.3f} - {d['time_end'].values[0]:0.3f}) R-Score = {r_score}",
                                 *print_colors[status][1], attrs=['bold'], end='')
                print(warning)

        if len(score_df[score_df['window_type'] == 'interpeak']) > 0:
            if verbose:
                print("""
Signal Reconstruction of Interpeak Windows
==========================================
                  """)
            for g, d in score_df[score_df['window_type'] == 'interpeak'].groupby('window_id'):
                status = d['status'].values[0]
                if status == 'valid':
                    warning = ''
                elif status == 'needs review':
                    warning = f"""
Interpeak window {g} is not well reconstructed by mixture, but has a small Fano factor  
compared to peak region(s). This is likely acceptable, but visually check this region.\n"""
                elif status == 'invalid':
                    warning = f"""
Interpeak window {g} is not well reconstructed by mixture and has an appreciable Fano 
factor compared to peak region(s). This suggests you have missed a peak in this 
region. Consider adding manual peak positioning by passing `known_peaks` 
to `fit_peaks()`."""

                if (d['reconstruction_score'].values[0] >= 10) | \
                        (d['reconstruction_score'].values[0] <= 0.1):
                    if d['reconstruction_score'].values[0] == 0:
                        r_score = f'0'
                    else:
                        r_score = f"10^{int(np.log10(d['reconstruction_score'].values[0]))}"
                else:
                    r_score = f"{d['reconstruction_score'].values[0]:0.4f}"
                if ((d['signal_fano_factor'].values[0] / mean_fano) > 10) | \
                   ((d['signal_fano_factor'].values[0] / mean_fano) <= 1E-5):
                    if (d['signal_fano_factor'].values[0] / mean_fano) == 0:
                        fano = '0'
                    else:
                        fano = f"10^{int(np.log10(d['signal_fano_factor'].values[0] / mean_fano))}"
                else:
                    fano = f"{d['signal_fano_factor'].values[0] / mean_fano:0.4f}"

                if verbose:
                    termcolor.cprint(f"{print_colors[status][0]} Interpeak Window {int(g)} (t: {d['time_start'].values[0]:0.3f} - {d['time_end'].values[0]:0.3f}) R-Score = {r_score} & Fano Ratio = {fano}",
                                     *print_colors[status][1], attrs=['bold'], end='')
                    print(warning)
        if verbose:
            print("""
--------------------------------------------------------------------------------""")
        return score_df

    def show(self,
             time_range=[]):
        """
        Displays the chromatogram with mapped peaks if available.

        Parameters
        ----------
        time_range : `List`
            Adjust the limits to show a restricted time range. Should 
            be provided as two floats in the range of [`lower`, `upper`]. Note
            that this does not affect the chromatogram directly as in `crop`. 


        Returns
        -------
        fig : `matplotlib.figure.Figure`
            The matplotlib figure object.
        ax : `matplotlib.axes._axes.Axes`
            The matplotlib axis object.
        """
        sns.set()

        # Set up the figure
        fig, ax = plt.subplots(1, 1)
        ax.set_xlabel(self.time_col)
        if self._bg_corrected:
            self._viz_ylabel_subtraction_indication = True
            ylabel = f"{self.int_col.split('_corrected')[0]} (baseline corrected)"
        else:
            self._viz_ylabel_subtraction_indication = False
            ylabel = self.int_col
        ax.set_ylabel(ylabel)

        # Plot the raw chromatogram
        ax.plot(self.df[self.time_col], self.df[self.int_col], 'k-',
                label='raw chromatogram')

        # Compute the skewnorm mix
        self._viz_min_one_concentration = None
        if self.peaks is not None:
            self._viz_peak_reconstruction = True
            time = self.df[self.time_col].values
            # Plot the mix
            convolved = np.sum(self.unmixed_chromatograms, axis=1)
            ax.plot(time, convolved, 'r--', label='inferred mixture')
            for g, d in self.peaks.groupby('peak_id'):
                label = f'peak {int(g)}'
                if self._mapped_peaks is not None:
                    self._viz_mapped_peaks = True
                    if g in self._mapped_peaks.keys():
                        d = self.quantified_peaks[self.quantified_peaks['compound']
                                                  == self._mapped_peaks[g]]
                        if 'concentration' in d.keys():
                            self._viz_min_one_concentration = True
                            label = f"{self._mapped_peaks[g]}\n[{d.concentration.values[0]:0.3g}"
                        else:
                            if self._viz_min_one_concentration is None:
                                self._viz_min_one_concentration = False
                            label = f"{self._mapped_peaks[g]}"
                        if 'unit' in d.keys():
                            self._viz_unit_display = True
                            label += f" {d['unit'].values[0]}]"
                        else:
                            self._viz_unit_display = False
                            label += ']'
                    else:
                        self._viz_mapped_peaks = False
                        label = f'peak {int(g)}'
                else:
                    self._viz_mapped_peaks = False
                ax.fill_between(time, self.unmixed_chromatograms[:, int(g) - 1], label=label,
                                alpha=0.5)
        else:
            self._viz_peak_reconstruction = False
        if 'estimated_background' in self.df.keys():
            self._viz_subtracted_baseline = True
            ax.plot(self.df[self.time_col], self.df['estimated_background'],
                    color='dodgerblue', label='estimated background', zorder=1)
        else:
            self._viz_subtracted_baseline = False

        if self._added_peaks is not None:
            ymax = ax.get_ylim()[1]
            for i, loc in enumerate(self._added_peaks):
                if i == 0:
                    label = 'suggested peak location'
                else:
                    label = '__nolegend__'
                ax.vlines(loc, 0, ymax, linestyle='--',
                          color='dodgerblue', label=label)
        ax.legend(bbox_to_anchor=(1.5, 1))
        fig.patch.set_facecolor((0, 0, 0, 0))
        if len(time_range) == 2:
            self._viz_adjusted_xlim = True
            ax.set_xlim(time_range)
            # Determine the max min and max value of the chromatogram within range.
            _y = self.df[(self.df[self.time_col] >= time_range[0]) & (
                self.df[self.time_col] <= time_range[1])][self.int_col].values
            ax.set_ylim([ax.get_ylim()[0], 1.1 * _y.max()])
        else:
            self._viz_adjusted_xlim = False
        return [fig, ax]
