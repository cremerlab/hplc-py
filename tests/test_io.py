# %%
import numpy as np
import hplc.io


def test_load_chromatogram():
    """
    Tests that chromatogram loading acts as advertised and handles only one
    chromatogram at a time.
    """
    bad_cols = ['column1', '2nmuloc']
    columns_list = ['column1', 'column2']
    columns_dict = {'column1': 'col1', 'column2': 'col2'}

    # Ensure that loading aborts if not all provided column names are found
    try:
        df = hplc.io.load_chromatogram('./tests/test_data/valid_chrom.txt',
                                       bad_cols)
        assert False
    except ValueError:
        assert True

    # Ensure that invalid chromatograms are rejected
    try:
        df = hplc.io.load_chromatogram(
            './tests/test_data/invalid_chrom1.txt', columns_list)
        assert False
    except RuntimeError:
        assert True
    try:
        df = hplc.io.load_chromatogram(
            './tests/test_data/invalid_chrom2.txt', columns_list)
        assert False
    except ValueError:
        assert True

    # Load a valid chromatogram
    df = hplc.io.load_chromatogram(
        './tests/test_data/valid_chrom.txt', columns_list)
    # Ensure that only provided columns are provided
    assert np.array([c in df.keys() for c in columns_list]).all()
    assert 'column3' not in df.keys()

    # Test that empty line dropping is being performed as expected
    df = hplc.io.load_chromatogram(
        './tests/test_data/valid_chrom.txt', columns_list, dropna=False)
    assert df.isnull().any().any()

    # Ensure that columns are renamed if a dictionary is provided
    df = hplc.io.load_chromatogram(
        './tests/test_data/valid_chrom.txt', columns_dict)
    assert np.array([c in df.keys() for c in columns_dict.values()]).all()
