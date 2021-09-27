import numpy as np
import os
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import shutil
import tempfile
import urllib
import warnings
import zipfile



def check_numeric(data, dtype="float", copy=None, suppress_errors=False):
    """Check a matrix contains only numeric data.

    Parameters
    ----------
    data : array-like, shape=[n_samples, n_features]
        Input data
    dtype : str or `np.dtype`, optional (default: 'float')
        Data type to which to coerce the data
    copy : bool or None, optional (default: None)
        Copy the data before coercion. If None, default to
        False for all datatypes except pandas.SparseDataFrame
    suppress_errors : bool, optional (default: False)
        Suppress errors from non-numeric data

    Returns
    -------
    data : array-like, shape=[n_samples, n_features]
        Output data as numeric type

    Raises
    ------
    TypeError : if `data` cannot be coerced to `dtype`
    """
    if copy is None:
        copy = is_SparseDataFrame(data)
    try:
        return data.astype(dtype, copy=copy)
    except TypeError as e:
        if is_SparseDataFrame(data):
            if not copy:
                raise TypeError(
                    "pd.SparseDataFrame does not support "
                    "copy=False. Please use copy=True."
                )
            else:
                return data.astype(dtype)
        else:
            raise e
    except ValueError:
        if suppress_errors:
            warnings.warn(
                "Data is not numeric. Many scprep functions will not work.",
                RuntimeWarning,
            )
            return data
        else:
            raise

def is_SparseDataFrame(X):
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            "The SparseDataFrame class is removed from pandas. Accessing it from the "
            "top-level namespace will also be removed in the next version",
            FutureWarning,
        )
        try:
            return isinstance(X, pd.SparseDataFrame)
        except AttributeError:
            return False