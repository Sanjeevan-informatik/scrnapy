from scipy import sparse

import numbers
import numpy as np
import pandas as pd
import warnings


def matrix_transform(data, fun, *args, **kwargs):
    """Perform a numerical transformation to data.

    Parameters
    ----------
    data : array-like, shape=[n_samples, n_features]
        Input data
    fun : callable
        Numerical transformation function, `np.ufunc` or similar.
    args, kwargs : additional arguments, optional
        arguments for `fun`. `data` is always passed as the first argument

    Returns
    -------
    data : array-like, shape=[n_samples, n_features]
        Transformed output data
    """
    if is_sparse_dataframe(data) or is_SparseDataFrame(data):
        data = data.copy()
        for col in data.columns:
            data[col] = fun(data[col], *args, **kwargs)
    elif sparse.issparse(data):
        if isinstance(data, (sparse.lil_matrix, sparse.dok_matrix)):
            data = data.tocsr()
        else:
            # avoid modifying in place
            data = data.copy()
        data.data = fun(data.data, *args, **kwargs)
    else:
        data = fun(data, *args, **kwargs)
    return data

def sqrt(data):
    """Square root transform.

    Parameters
    ----------
    data : array-like, shape=[n_samples, n_features]
        Input data

    Returns
    -------
    data : array-like, shape=[n_samples, n_features]
        Square root transformed output data

    Raises
    ------
    ValueError : if data has negative values
    """

    return matrix_transform(data, np.sqrt)


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


def is_sparse_dataframe(x):
    if isinstance(x, pd.DataFrame) and not is_SparseDataFrame(x):
        try:
            x.sparse
            return True
        except AttributeError:
            pass
    return False



data = pd.read_csv(r"C:\Users\maky\Desktop\SRNAPY\data\test_10X\output\data.csv")  

index_data= data["0"]
del data['0']

data = data.set_index([pd.Index(index_data)])

data=sqrt(data)

print(data)

data.to_csv(r"C:\Users\maky\Desktop\SRNAPY\data\test_10X\output\data.csv")
