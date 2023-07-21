#%%
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 
import hplc.quant

chrom_df = pd.read_csv('../tests/test_chrom.csv')
peak_df = pd.read_csv('../tests/test_peak_table.csv')

props = ['retention_time', 'amplitude', 'area', 'scale', 'skew']
for g, d in chrom_df.groupby('iter'):
    if g == chrom_df.iter.max():
        chrom = hplc.quant.Chromatogram(d, cols={'time': 'x', 'intensity':'y'}, 
                                        bg_subtract=False)
        peaks = chrom.quantify(prominence=1E-4)
        assert len(peaks) == peak_df['peak_idx'].max()

        break

#%%
def test_peak_fitting():
    # Load test data
    chrom_df = pd.read_csv('../tests/test_chrom.csv')
    peak_df = pd.read_csv('../tests/test_peak_table.csv')

    # Define constants used in each test
    props = ['retention_time', 'amplitude', 'area', 'scale', 'skew']
    colnames = {'time':'x', 'intensity':'y'}
    n_peaks = peak_df['peak_idx'].max()
    tol = 1E-3
    for g, d in chrom_df.groupby('iter'):
        # Select true peak info
        truth = peak_df[peak_df['iter']==g]

        # Execute analysis
        chrom = hplc.quant.Chromatogram(d, bg_subtract=False, cols=colnames)
        peaks = chrom.quantify(prominence=1E-3, verbose=False)

        ## APPLY TESTING
        # Enforce that the correct number of peaks have been identified
        assert len(peaks) == n_peaks

        # Enforce that the inferred peak parameters are within a tolerance of 0.1%
        for p in props:
            print(g)
            assert np.isclose(peaks[p], truth[p], rtol=tol).all()
    



test_peak_fitting()