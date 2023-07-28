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

class Chromatogram(object):
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
        A Pandas Dataframe with peak areas converted to 
    param_opt : `numpy.ndarray`
        An array of the parameter estimates in order of amplitude, location, scale, 
        and skew for each peak in temporal order. 
    param_pcov: 2-D `numpy.ndarray`
        The estimated approximate covariance matrix of the parameters. Uncertainty
        for each parameter can be calculated as `numpy.sqrt(numpy.diag(param_pcov))`,
        with the following big caveat:

        .. attention::
            `param.pcov` is only an *estimate* of the *approximate* covariance
            matrix and computation of the error is only valid if the linear 
            approximation to the model about the optimum is valid. Use this 
            attribute with caution.
 
    """
    def __init__(self, file, time_window=None, 
                 cols={'time':'time', 'signal':'signal'}):
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
            a Pandas Dataframe of the chromatogram to analyze. If None, 
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
        if (type(file) is not str) & (type(file) is not pd.core.frame.DataFrame):
            raise RuntimeError(f'Argument must be either a filepath or pandas DataFrame. Argument is of type {type(file)}')
        if (time_window is not None):
            if type(time_window) != list:
                raise TypeError(f'`time_window` must be of type `list`. Type {type(time_window)} was proivided')
            if len(time_window) != 2:
                raise ValueError(f'`time_window` must be of length 2 (corresponding to start and end points). Provided list is of length {len(time_window)}.')

        # Assign class variables 
        self.time_col = cols['time']
        self.int_col = cols['signal']

        # Load the chromatogram and necessary components to self. 
        dataframe = file.copy()
        self.df = dataframe

        # Define the average timestep in the chromatogram. This computes a mean
        # but values will typically be identical.
        self._dt = np.mean(np.diff(dataframe[self.time_col].values))

        # Prune to time window
        if time_window is not None:
            self.crop(time_window)
        else: 
            self.df = dataframe

        # Blank out vars that are used elsewhere
        self.window_props = None
        self._peak_indices = None
        self.peaks = None
        self._guesses = None
        self._bg_corrected = False
        self._mapped_peaks = None
        self._added_peaks = None
        self.unmixed_chromatograms = None

    def crop(self, time_window=None, return_df=False):
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
        if type(time_window) != list:
                raise TypeError(f'`time_window` must be of type `list`. Type {type(time_window)} was proivided')
        if len(time_window) != 2:
                raise ValueError(f'`time_window` must be of length 2 (corresponding to start and end points). Provided list is of length {len(time_window)}.')
        self.df = self.df[(self.df[self.time_col] >= time_window[0]) & 
                          (self.df[self.time_col] <= time_window[1])]
        if return_df:
            return self.df

    def _assign_windows(self, enforced_locations=[], enforced_widths=[], 
                             enforcement_tolerance=0.5, 
                             prominence=0.01, rel_height=0.95, buffer=100):
        R"""
        Breaks the provided chromatogram down to windows of likely peaks. 

        Parameters
        ----------
        enforced_locations : `list`
            The approximate locations of the peaks. If this is not provided, 
            peak locations will be automatically detected. 
        enforced_widths : `list`
            The approximate widths of the peaks. If this is not provided  but
            `locations` is, approximate widths of one time unit (0.5 / dt)
            will be assumed.
        enforce_tolerance: `float`, optional
            If an enforced peak location is within tolerance of an automatically 
            identified peak, the automatically identified peak will be preferred. 
            This parameter is in units of time. Default is one-half time unit.
        prominence : `float`,  [0, 1]
            The promimence threshold for identifying peaks. Prominence is the 
            relative height of the normalized signal relative to the local
            background. Default is 1%.
        rel_height : `float`, [0, 1]
            The relative height of the peak where the baseline is determined. 
            Default is 95%.
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
        for param, param_name, param_type in zip([prominence, rel_height, buffer], 
                                     ['prominence', 'rel_height',  'buffer'],
                                     [float, float, int]):
            if type(param) is not param_type:
                raise TypeError(f'Parameter {param_name} must be of type `{param_type}`. Type `{type(param)}` was supplied.') 
        if (prominence < 0) | (prominence > 1):
            raise ValueError(f'Parameter `prominence` must be [0, 1].')
        if (rel_height < 0) | (rel_height > 1):  
            raise ValueError(f' `rel_height` must be [0, 1].')
        if (buffer < 0):
            raise ValueError('Parameter `buffer` cannot be less than 0.')

        # Correct for a negative baseline 
        df = self.df
        intensity = self.df[self.int_col].values
        norm_int = (intensity - intensity.min()) / (intensity.max() - intensity.min())

        # Preform automated peak detection and set window ranges
        peaks, _ = scipy.signal.find_peaks(norm_int, prominence=prominence)
        self._peak_indices = peaks
        _widths, _, _, _ = scipy.signal.peak_widths(intensity, self._peak_indices, 
                                       rel_height=0.5)
        # Compute the peak widths  
        _, _, _left, _right = scipy.signal.peak_widths(intensity, self._peak_indices, 
                                       rel_height=1) 

        # Set window ranges
        ranges = []
        for l, r in zip(_left, _right):
            _range = np.arange(int(l - buffer), int(r + buffer), 1)
            _range = _range[(_range >= 0) & (_range <= len(norm_int))]
            ranges.append(_range)
        self._ranges = ranges
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
        self.ranges = ranges
        # If manual locations are provided, ensure that they are identified        
        if len(enforced_locations) != 0:
            enforced_location_inds = np.int_(np.array(enforced_locations) / self._dt)

            # Keep track of what locations and widths need to be added.
            added_peaks = []
            added_peak_inds = []
            for i, loc in enumerate(enforced_location_inds):
                if np.sum(np.abs(self._peak_indices - loc) < enforcement_tolerance/self._dt) == 0:
                    added_peaks.append(loc) 
                    added_peak_inds.append(i)

            # Consider the edge case where all enforced locations have been automatically detected
            if len(added_peaks) > 0:
                self._peak_indices = np.append(self._peak_indices, added_peaks)                
                self._added_peaks = self.df[self.time_col].values[added_peaks]
                if len(enforced_widths) == 0:
                    _enforced_widths = np.ones_like(added_peaks) / self._dt
                else:
                    if type(enforced_widths) == list:
                        enforced_widths = np.array(enforced_widths)
                    _enforced_widths = enforced_widths[added_peak_inds]  / self._dt
                _widths = np.append(_widths, _enforced_widths)

                # Ensure that the newly added peak is within one of the identified 
                # ranges
                orphan_locs = []
                orphan_widths = [] 
                for i, p in enumerate(added_peaks):
                    located = False 
                    for j, r in enumerate(ranges):
                        if p in r:
                            located = True
                            break
                    if not located:
                        orphan_locs.append(p)
                        orphan_widths.append(_enforced_widths[i])
                # If there are orphan peaks, create ranges                    
                if len(orphan_locs) > 0:
                    for i, o in enumerate(orphan_locs):
                        _added_range = np.arange(o - orphan_widths[i] - buffer, o + orphan_widths[i] + buffer, 1)
                        _added_range = _added_range[(_added_range >= 0) & (_added_range <= len(norm_int))]
                        ranges.append(_added_range)
            else:    
                warnings.warn(f'All manually provided peaks are within {enforcement_tolerance} of an automatically identified peak. If this location is desired, decrease value of `enforcement_tolerance`.')

        # Copy the dataframe and return the windows
        window_df = df.copy(deep=True)
        window_df.sort_values(by=self.time_col, inplace=True)
        window_df['time_idx'] = np.arange(len(window_df))
        window_df['window_id'] = 0
        window_df['window_type'] = 'peak'
        for i, r in enumerate(ranges):
            window_df.loc[window_df['time_idx'].isin(r), 
                                    'window_id'] = int(i + 1)

        # Determine the windows for the background (nonpeak) areas.
        bg_windows = window_df[window_df['window_id']==0]
        tidx = bg_windows['time_idx'].values

        if len(bg_windows) > 0:
            split_inds = np.nonzero(np.diff(bg_windows['time_idx'].values) - 1)[0]

            # If there is only one background window
            if split_inds[0] == 0:
                window_df.loc[window_df['time_idx'].isin(bg_windows['time_idx'].values), 'window_id'] = 1 
                window_df.loc[window_df['time_idx'].isin(bg_windows['time_idx'].values), 'window_type'] = 'background'

            # If more than one split ind, set up all ranges.
            if split_inds[0] != 0:
                split_inds += 1
                split_inds = np.insert(split_inds, 0, 0)
                split_inds = np.append(split_inds, len(tidx)) 

            bg_ranges = [bg_windows.iloc[np.arange(split_inds[i], split_inds[i+1], 1)]['time_idx'].values for i in range(len(split_inds)-1)]
            for i, rng in enumerate(bg_ranges):
                window_df.loc[window_df['time_idx'].isin(rng), 'window_id'] = i + 1
                window_df.loc[window_df['time_idx'].isin(rng), 'window_type'] = 'background' 
        
        # Convert this to a dictionary for easy parsing
        window_dict = {}
        for g, d in window_df[window_df['window_type']=='peak'].groupby('window_id'):
            if g > 0:
                _peaks = [p for p in self._peak_indices if p in d['time_idx'].values]
                peak_inds = [x for _p in _peaks for x in np.where(self._peak_indices == _p)[0]]
                _dict = {'time_range':d[self.time_col].values,
                         'signal': d[self.int_col].values,
                         'signal_area': d[self.int_col].values.sum(),
                         'num_peaks': len(_peaks),
                         'amplitude': [d[d['time_idx']==p][self.int_col].values[0] for p in _peaks],
                         'location' : [d[d['time_idx']==p][self.time_col].values[0] for p in _peaks],
                         'width' :  [_widths[ind] * self._dt for ind in peak_inds]
                             }
                window_dict[int(g)] = _dict

        self.window_df = window_df  
        self.window_props = window_dict

        return window_df  

    def _compute_skewnorm(self, x, *params):
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
                    The skew shape parater of the distribution.

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
        norm = np.sqrt(2 * np.pi * scale**2)**-1 * np.exp(-(x - loc)**2 / (2 * scale**2))
        cdf = 0.5 * (1 + scipy.special.erf(_x / np.sqrt(2))) 
        return amp * 2 * norm * cdf

    def _fit_skewnorms(self, x, *params):
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


    def deconvolve_peaks(self, verbose=True, param_bounds={}, max_iter=1000000, **optimizer_kwargs):
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
        param_bounds : `dict`, optional
            Modifications to the default parameter bounds (see Notes below) as 
            a dictionary for each parameter. A dict entry should be of the 
            form `parameter: [lower, upper]`. Modifications have the following effects:
            + Modifications to `amplitude` bounds are multiplicative of the 
              observed magnitude at the peak position. 
            + Modifications to `location` are values that are subtracted or 
              added from the peak position for lower and upper bounds, respectively.
            + Modifications to `scale` replace the default values. 
            + Modifications to `skew` replace the default values. 
        max_iter : int
            The maximum number of iterations the optimization protocol should 
            take before erroring out. Default value is 10^6.
        optimizer_kwargs : dict
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

            + `amplitude`: The lower and upper peak amplitude boundaries correspond 
            to one-tenth and ten-times the value of the peak at the peak location 
            in the chromatogram.

            + `location`: The lower and upper location bounds correspond to the 
            minimum and maximum time values of the chromatogram.

            + `scale`: The lower and upper bounds of the peak standard deviation
            defaults to the chromatogram time-step and one-half of the chromatogram
            duration, respectively.  

            + `skew`: The skew parameter by default is allowed to take any value
            between (-5, 5).
        """ 
        if self.window_props is None:
            raise RuntimeError('Function `_assign_windows` must be run first. Go do that.')
        if verbose:
            iterator = tqdm.tqdm(self.window_props.items(), desc='Deconvolving mixture')  
        else:
            iterator = self.window_props.items()
        if (len(param_bounds)) > 0 & (param_bounds.keys() not in ['amplitude', 'location', 'scale', 'skew']):
            raise ValueError(f"`param_bounds` must have keys of `amplitude`, `location`, `scale`, and `skew`. Provided keys are {param_bounds.keys()}")
        peak_props = {}
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

            for i in range(v['num_peaks']):
                p0.append(v['amplitude'][i])
                p0.append(v['location'][i]),
                p0.append(v['width'][i] / 2) # scale parameter
                p0.append(0) # Skew parameter, starts with assuming Gaussian

                if len(param_bounds) == 0:
                    # Lower bounds
                    bounds[0].append(0.1 * v['amplitude'][i]) 
                    bounds[0].append(v['time_range'].min()) 
                    bounds[0].append(self._dt) 
                    bounds[0].append(-np.inf) 
                    # Upper bounds
                    bounds[1].append(10 * v['amplitude'][i])
                    bounds[1].append(v['time_range'].max())
                    bounds[1].append((v['time_range'].max() - v['time_range'].min())/2)
                    bounds[1].append(np.inf)
                else:
                    bounds[0].append(param_bounds['amplitude'][0] * v['amplitude'][i]) 
                    bounds[0].append(v['location'] - param_bounds['location'][0]) 
                    bounds[0].append(param_bounds['scale'][0]) 
                    bounds[0].append(param_bounds['skew'][0]) 
                    # Upper bounds
                    bounds[1].append(param_bounds['amplitude'][1] * v['amplitude'][i])
                    bounds[1].append(v['location'] + param_bounds['location'][1])
                    bounds[1].append(param_bounds['scale'][1])
                    bounds[1].append(param_bounds['skew'][1]) 
            
            # Perform the inference
            popt, pcov = scipy.optimize.curve_fit(self._fit_skewnorms, v['time_range'],
                                               v['signal'], p0=p0, bounds=bounds, maxfev=max_iter,
                                               **optimizer_kwargs)
            self.param_opts = popt
            self.param_cov = pcov

            # Assemble the dictionary of output 
            if v['num_peaks'] > 1:
                popt = np.reshape(popt, (v['num_peaks'], 4)) 
            else:
                popt = [popt]
            for i, p in enumerate(popt):
                window_dict[f'peak_{i + 1}'] = {
                            'amplitude': p[0],
                            'retention_time': p[1],
                            'scale': p[2],
                            'alpha': p[3],
                            'area':self._compute_skewnorm(v['time_range'], *p).sum(),
                            'reconstructed_signal':self._compute_skewnorm(v['time_range'], *p)}
            peak_props[k] = window_dict
         
        self._peak_props = peak_props
        return peak_props


    def fit_peaks(self, enforced_locations=[], enforced_widths=[],
                  enforcement_tolerance=0.5, prominence=1E-2, rel_height=1.0, 
                  approx_peak_width=3, buffer=100, param_bounds={}, verbose=True, return_peaks=True, 
                 correct_baseline=True, max_iter=1000000, **optimizer_kwargs):
        R"""
        Detects and fits peaks present in the chromatogram

        Parameters
        ----------
        enforced_locations : `list`
            The approximate locations of the peaks. If this is not provided, 
            peak locations will be automatically detected. 
        enforced_widths : `list`
            The approximate widths of the peaks. If this is not provided  but
            `locations` is, approximate widths of one time unit (1 / dt)
            will be assumed.
        enforce_tolerance: `float`, optional
            If an enforced peak location is within tolerance of an automatically 
            identified peak, the automatically identified peak will be preferred. 
            This parameter is in units of time. Default is one-half time unit.
        prominence : `float`,  [0, 1]
            The promimence threshold for identifying peaks. Prominence is the 
            relative height of the normalized signal relative to the local
            background. Default is 1%. If `locations` is provided, this is 
            not used.
        rel_height : `float`, [0, 1]
            The relative height of the peak where the baseline is determined. 
            Default is 100%. If `locations` is provided, this is not used.
        approx_peak_width: `float`, optional
            The approximate width of the signal you want to quantify. This is 
            used as filtering window for automatic baseline correction. If `correct_baseline==False`,
            this has no effect. 
        buffer : positive `int`
            The padding of peak windows in units of number of time steps. Default 
            is 100 points on each side of the identified peak window. If `locations` 
            is provided, this is not used.
        verbose : `bool`
            If True, a progress bar will be printed during the inference. 
        param_bounds: `dict`, optional
            Parameter boundary modifications to be used to constrain fitting. 
            See docstring of :func:`~hplc.quant.Chromatogram.deconvolve_peaks`
            for more information.
        correct_baseline : `bool`, optional
            If True, the baseline of the chromatogram will be automatically 
            corrected using the SNIP algorithm. See :func:`~hplc.quant.Chromatogram.correct_baseline`
            for more information.
        return_peaks : `bool`, optional
            If True, a dataframe containing the peaks will be returned. Default
            is True.
        max_iter : int
            The maximum number of iterations the optimization protocol should 
            take before erroring out. Default value is 10^6.
        **optimizer_kwargs : `dict`
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
            self.correct_baseline(window=approx_peak_width, verbose=verbose, return_df=False)

        # Assign the window bounds
        _ = self._assign_windows(enforced_locations=enforced_locations, 
                                      enforced_widths=enforced_widths,
                                      enforcement_tolerance=enforcement_tolerance,
                                      prominence=prominence, rel_height=rel_height, 
                                      buffer=buffer)

        # Infer the distributions for the peaks
        peak_props = self.deconvolve_peaks(verbose=verbose, param_bounds=param_bounds, max_iter=max_iter, 
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
        for _ , _v in self._peak_props.items():
            for _, v in _v.items():
                params = [v['amplitude'], v['retention_time'], 
                          v['scale'], v['alpha']]
                out[:, iter] = self._compute_skewnorm(time, *params)
                iter += 1
        self.unmixed_chromatograms = out

        # Assess the quality of the fit and throw a cased warning 
        if return_peaks:
            return peak_df
    
    def correct_baseline(self, window=3, return_df=False, verbose=True):
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

        Returns
        -------
        corrected_df : `pandas.core.frame.DataFrame`
            If `return_df = True`, then the original and the corrected chromatogram are returned.

        Notes
        -----
        This implements the SNIP algorithm as presented and summarized in `Morhác
        and Matousek 2008 <https://doi.org/10.1366/000370208783412762>`_.
        """

        # Unpack and copy dataframe and intensity profile
        df = self.df
        signal = df[self.int_col].copy()

        # Ensure positivity of signal
        signal *= np.heaviside(signal, 0)

        # Compute the LLS operator
        tform = np.log(np.log(np.sqrt(signal.values + 1) + 1) + 1)

        # Compute the number of iterations given the window size.
        n_iter = int(((window / self._dt) - 1) / 2)

        # Iteratively filter the signal
        if verbose:
            iter = tqdm.tqdm(range(1, n_iter), desc="Performing baseline correction")
        else:
            iter = range(1, n_iter)
        for i in iter:
            tform_new = tform.copy()
            for j in range(i, len(tform) - i):
                tform_new[j] = min(tform_new[j], 0.5 * (tform_new[j+i] + tform_new[j-i])) 
            tform = tform_new

        # Perform the inverse of the LLS transformation and subtract
        inv_tform = (np.exp(np.exp(tform) - 1) - 1)**2 - 1 
        df[self.int_col] = signal - inv_tform
        df[f'estimated_background'] = inv_tform 
        self.df = df  
        if return_df:
            return df

    def map_peaks(self, params, loc_tolerance=0.5, include_unmapped=False):
        """
        Maps user-provided mappings to arbitrarily labeled peaks. If a linear 
        calibration curve is also provided, the concentration will be computed.

        .. note::
            As of `v0.1.0`, this function can only accommodate linear calibration 
            functions.

        Parameters
        ----------
        params : `dict` of `dict`s
            A dictionary mapping each peak to a slope and intercept used for 
            converting peak areas to units of concentraions. Each peak 
            should have a key that is the compound name (e.g. "glucose"). Each 
            key should have another dict as the key with `retention_time`, `slope`, and `intercept`
            as keys. If only `retention_time` is given, concentration will 
            not be computed. The key `retention_time` will be used to map the compound to the 
            `peak_id`. If `unit` are provided, this will be added as a column
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
        """
        # Create a mapper for peak id to compound
        mapper = {}
        peak_df = self.peaks.copy()
        for k, v in params.items():
            ret_time = v['retention_time']
            peak_id = np.abs(peak_df['retention_time'].values - ret_time) < loc_tolerance

            if np.sum(peak_id) > 1:
                raise ValueError(f"Multiple compounds found within tolerance of retention time for {k}. Reduce the tolerance or correct the provided value.")

            if np.sum(peak_id) == 0:
                warnings.warn(f"\nNo peak found for {k} (retention time {v['retention_time']}) within the provided tolerance.")
                break
            peak_id = peak_df.peak_id.values[np.argmax(peak_id)] 
            peak_df.loc[peak_df['peak_id']==peak_id, 'compound'] = k
            mapper[peak_id] = k
        if len(mapper) == 0:
            raise ValueError("No peaks could be properly mapped! Check your provided retention times.")

        # Iterate through the compounds and calculate the concentration. 
        for g, d in peak_df.groupby('compound'):
            if (g in params.keys()):
                if 'slope' in params[g].keys():
                    conc = (d['area'] - params[g]['intercept']) / params[g]['slope']
                    peak_df.loc[peak_df['compound']==g, 'concentration'] = conc
                    if 'unit' in params[g].keys():
                        peak_df.loc[peak_df['compound']==g, 'unit'] = params[g]['unit']
        if include_unmapped == False:
            peak_df.dropna(inplace=True)
        self.quantified_peaks = peak_df
        self._mapped_peaks = mapper
        return peak_df

    def score_reconstruction(self):
        """
        Computes the reconstruction score on a per-window and total chromatogram
        basis.

        Parameters
        ----------
        tol : `float`
            The tolerance for a reconstruction to be valid. This is the tolerated 
            deviation from a score of 1 which indicates a perfectly reconstructed
            chromatogram. 

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
            R = \frac{\log\left[\sum\limits_{i\in t}^t \sum\limits_{j \in N_\text{peaks}}^{N_\text{peaks}}2A_j \text{SkewNormal}(\alpha_j, r_{t_j}, \sigma_j)\right]}{\log\left[\sum\limits_{i \in t}^t S_i\right]}
           
        where :math:`t` is the total time of the region, :math:`A`is the inferred 
        peak amplitude, :math:`\alpha` is the inferred skew paramter, :math:`r_t` is
        the inferred peak retention time, :math:`\sigma` is the inferred scale 
        parameter and :math:`S_i` is the observed signal intensity at time point
        :math:`i`. 

        """

        # tform = np.log(np.log(np.sqrt(signal.values + 1) + 1) + 1)
        columns = ['window_id', 'time_start', 'time_end', 'signal_area', 
                   'inferred_area', 'signal_fano_factor', 'reconstruction_score']
        score_df = pd.DataFrame([])  
        # Compute the per-window reconstruction

        for g, d in self.window_df[self.window_df['window_type']=='peak'].groupby('window_id'):
            # Compute the non-peak windows separately.
            window_area = d[self.int_col].values.sum()
            window_peaks = self._peak_props[g]
            window_peak_area = np.array([v['reconstructed_signal'] for v in window_peaks.values()]).sum()
            score = np.array(window_peak_area / window_area).astype(float)
            signal_fano = np.array(np.var(d[self.int_col].values) / np.mean(d[self.int_col].values)).astype(float)
            x = np.array([g, d[self.time_col].min(),  
                            d[self.time_col].max(), window_area, 
                            window_peak_area, signal_fano, score]) 
            _df = pd.DataFrame({_c:_x for _c, _x in zip(columns, x)}, index=[g - 1])
            _df['window_type'] = 'peak'
            score_df = pd.concat([score_df, _df]) 
         
        # Compute the score for the non-peak regions
        nonpeak = self.window_df[self.window_df['window_type'] == 'background']
        if len(nonpeak) > 0:
            for g, d in nonpeak.groupby('window_id'):
                total_area = d[self.int_col].values.sum()
                recon_area = np.sum(self.unmixed_chromatograms, axis=1)[d['time_idx'].values].sum()
                nonpeak_score = recon_area / total_area
                signal_fano = np.var(self.df[self.int_col].values[d['time_idx'].values]) / np.mean(self.df[self.int_col].values[d['time_idx'].values])
                # Add to score dataframe
                x = [g, d[self.time_col].min(),
                    d[self.time_col].max(),
                    total_area, recon_area, signal_fano, nonpeak_score]
                _df = pd.DataFrame({c:xi for c, xi in zip(columns, x)}, index=[g - 1])
                _df['window_type'] = 'background'
                score_df = pd.concat([score_df, _df])
        score_df['signal_area'] = score_df['signal_area'].astype(float)
        score_df['inferred_area'] = score_df['inferred_area'].astype(float)
        self.scores = score_df
        return score_df

    def assess_fit(self, tol=1E-2, fano_tol=1E-3):
        """
        Assesses whether the computed reconstruction score is adequate, given a tolerance.

        Parameters
        ----------
        tol : `float`
            The tolerance for a reconstruction to be valid. This is the tolerated 
            deviation from a score of 1 which indicates a perfectly reconstructed
            chromatogram. 
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
        A reconstruction score of :math:`R = 1` indicates a perfect 
        reconstruction of the chromatogram. For practical purposes, a chromatogram
        is deemed to be adequately reconstructed if :math:`R` is within a tolerance
        :math:`\epsilon` of 1 such that

        .. math::
            \left| R - 1 \right| \leq \epsilon \Rightarrow \text{Valid Reconstruction}
        """

        if self.unmixed_chromatograms is None:
            raise RuntimeError("No reconstruction found! `.fit_peaks()` must be called first. Go do that.")

        # Compute the reconstruction score
        score_df = self.score_reconstruction() 

        # Apply the tolerance parameter
        score_df = self.scores.copy()
        score_df['applied_tolerance'] = tol
        score_df['accepted'] = np.abs(np.round(score_df['reconstruction_score'].values - 1, 
                                          decimals=int(np.abs(np.ceil(np.log10(tol)))))) <= tol 
        print("""
-------------------Chromatogram Reconstruction Report Card----------------------""")
        for g, d in score_df.groupby(['window_type','window_id']): 
            a = 1
        else:
            if score_df['window_id'].min() == 0:
                if (score_df[score_df['window_id'] > 0]['accepted']==True).all() \
                    & (score_df[score_df['window_id']==0]['accepted'].values[0]==False):
                    fano_ratio = score_df[score_df['window_id'] == 0]['signal_fano_factor'].values / score_df[score_df['window_id'] != 0]['signal_fano_factor'].values.mean()
                    if fano_ratio < fano_tol:
                        self._assessment = 1
                    else:
                        self._assessment = 2
                elif (score_df[score_df['window_id'] > 0]['accepted']==False).all():
                        self._assessment = 3
            else:
                self._assessment = 3 
 
#         # Print the summary to the screen
#         if verbose:
#             print("""
# -------------------Chromatogram Reconstruction Report Card----------------------""")
#             if self._assessment == 0:
#                 termcolor.cprint("""
# Success: Inferred chromatogram is a faithful reconstruction of the observed 
# signal.""", "black", "on_green", attrs=["bold"])
#             elif self._assessment == 1:
#                 termcolor.cprint("""
# Success (Mostly): Peak window chromatograms are faithfully reconstructed, but 
# poorly reconstructed in peak-free regions.""", "black", "on_yellow", attrs=["bold"])
#                 print("""
# This means that there are peak-free regions that are not well described by the inferred
# mixtures, but have a small Fano factor compared to detected peak regions. This 
# is probably fine, but call `Chromatogram.show()` to ensure you're not missing 
# a peak.""")
 
#             elif self._assessment == 2:
#                     termcolor.cprint("""
# Warning: Chromatogram only faithfully reconstructed within peak windows.
# """, "black", "on_magenta", attrs=["bold"])
#                     print("""
# This means that there are peak-free regions that are not well described by the
# inferred mixtures, but have a Fano factor above the noise. This may be due to
# pronounced systematic noise in the chromatogram or there are peaks which are not
# being detected. Call `Chromatogram.show()` and visually inspect where things
# have gone awry.""")
            
#             elif self._assessment == 3:
#                 bad_regions = score_df[(score_df['window_id'] > 0) &
#                                        (score_df['accepted']==False)]
#                 if len(bad_regions) > 1:
#                     num = f'{len(bad_regions)} out of {len(bad_regions) - 1} peak windows were'
#                 else:
#                     num  = 'The peak window was'
#                 termcolor.cprint(f"""
# !!!Failure: {num} not properly reconstructed!!!""", 
# "black", "on_red", attrs=["bold"])
#                 print("""
# The following regions are problematic. Try adjusting parameter bounds and peak 
# prominence filters or check if some locations need to be manually enforced (if 
# you have shouldered peaks, for example).
# """)
#                 for i in range(len(bad_regions)):
#                     info = bad_regions.iloc[i].to_dict()
#                     print(f"Window {info['window_id']} (t: {info['time_start']} - {info['time_end']}). R = {info['reconstruction_score']}")
#             print("""
# --------------------------------------------------------------------------------""")
#         self.scores = score_df
#         return score_df

 
    def show(self, time_range=[]):
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
        ax.set_ylabel(self.int_col)

        # Plot the raw chromatogram
        ax.plot(self.df[self.time_col], self.df[self.int_col], 'k-',
                label='raw chromatogram') 

        # Compute the skewnorm mix 
        if self.peaks is not None:
            time = self.df[self.time_col].values
            # Plot the mix
            convolved = np.sum(self.unmixed_chromatograms, axis=1)
            ax.plot(time, convolved, 'r--', label='inferred mixture') 
            for g, d in self.peaks.groupby('peak_id'):
                label = f'peak {int(g)}'
                if self._mapped_peaks is not None: 
                    if g in self._mapped_peaks.keys():
                        d = self.quantified_peaks[self.quantified_peaks['compound']==self._mapped_peaks[g]]
                        label = f"{self._mapped_peaks[g]}\n[{d.concentration.values[0]:0.3g}"
                        if 'unit' in d.keys():
                            label += f" {d['unit'].values[0]}]"
                        else:
                            label += ']'    
                    else:
                        label = f'peak {int(g)}'

                ax.fill_between(time, self.unmixed_chromatograms[:, int(g) - 1], label=label, 
                                alpha=0.5)
        if 'estimated_background' in self.df.keys():
            ax.plot(self.df[self.time_col], self.df['estimated_background'],  color='dodgerblue', label='estimated background', zorder=1)
        
        if self._added_peaks is not None:
            ymax = ax.get_ylim()[1]
            for loc in self._added_peaks:
                ax.vlines(loc, 0, ymax, linestyle='--', color='dodgerblue', label="suggested peak location")
        ax.legend(bbox_to_anchor=(1.5,1))
        fig.patch.set_facecolor((0, 0, 0, 0))
        if len(time_range) == 2:
          ax.set_xlim(time_range)
          # Determine the max min and max value of the chromatogram within range.
          _y = self.df[(self.df[self.time_col] >= time_range[0]) & (self.df[self.time_col] <= time_range[1])][self.int_col].values
          ax.set_ylim([ax.get_ylim()[0], 1.1 * _y.max()])

        return [fig, ax]