import numpy as np 
import pandas as pd

def scrape_chromatogram(fname, colnames=[], delimiter=',', dropna=False):
    R"""
    Parses a file containing a chromatogram and returns it as a Pandas DataFrame.

    Parameters 
    -----------
    fname: `str`
        The path to the file containing the chromatogram. This must be a text
        file (i.e. not `.xslx`!) 
    colnames : `list` or `dict`
        The desired columns present in the file. If provided as a dict, columns will
        be renamed as `key` -> `value`. If not provided, it will be assumed 
        that the chromatogram begins without needing to skip any lines. 
    delimiter : 'str' 
        The delimiter character separating columns in the chromatogram.
    dropna: `bool`
        If True, NaN's will be dropped from the chromatogram.  

    Returns
    -------
    df : `pandas.core.frame.DataFrame`
        The chromatogram loaded as a Pandas DataFrame with the desired columns.
    """
    if type(colnames) != dict:
        _colnames = list(colnames.keys())
    else:
        _colnames = colnames
    skip = 0
    if len(_colnames) != 0:
        with open(fname, 'r') as f:
            _lines = f.readlines()
            halted = False
            for line in _lines:   
                if np.array([nom.lower() in line.lower() for nom in _colnames]).all():
                    halted = True
                    break
                else:
                    skip += 1
            if not halted:
                raise ValueError("Provided column name(s) were not found in file. Check spelling?")

    # Given skips, load in dataframe and rename if necessary
    df = pd.read_csv(fname, skiprows=skip, delimiter=delimiter)
    if type(colnames) == dict:
        df.rename(columns=colnames, inplace=True)
        _colnames = list(colnames.values())
    if dropna:
        df.dropna(inplace=True)

    # Keep only the desired columns.
    df = df[_colnames]
    return df
    
