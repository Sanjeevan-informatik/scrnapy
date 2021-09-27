from scipy import sparse

import numbers
import numpy as np
import pandas as pd
import warnings



def toarray(x):
    """Convert an array-like to a np.ndarray.

    Parameters
    ----------
    x : array-like
        Array-like to be converted
    Returns
    -------
    x : np.ndarray
    """
    
    """
    type(x)
    print("x is the")
    print(x)
    print(type(x))
      """
    if isinstance(x, (pd.DataFrame, pd.Series, pd.Index)):
        x = x.to_numpy()
        #print(x)

    elif isinstance(x, list):
        x_out = []
        for xi in x:
            try:
                xi = toarray(xi)
            except TypeError:
                # recursed too far
                pass
            x_out.append(xi)
        # convert x_out from list to array
        x = np.array(x_out, dtype=_check_numpy_dtype(x_out))
    elif isinstance(x, (np.ndarray, numbers.Number)):
        pass
    else:
        raise TypeError("Expected array-like. Got {}".format(type(x)))
    return x

def gene_capture_count(data, cutoff=0):
    """Measure the number of cells in which each gene has non-negligible counts.

    Parameters
    ----------
    data : array-like, shape=[n_samples, n_features]
        Input data
    cutoff : float, optional (default: 0)
        Number of counts above which expression is deemed non-negligible

    Returns
    -------
    capture-count : list-like, shape=[m_features]
        Capturecount for each gene
    """
    #print(matrix_sum(data > cutoff, axis=0))
    gene_sums = np.array(matrix_sum(data > cutoff, axis=0)).reshape(-1)
    if isinstance(data, pd.DataFrame):
        gene_sums = pd.Series(gene_sums, index=data.columns, name="capture_count")
    return gene_sums


def _get_column_length(data):
    try:
        return data.shape[1]
    except (IndexError, AttributeError):
        return len(data)


def select_cols(
    data,
    idx=None,
    starts_with=None,
    ends_with=None,
    exact_word=None,
    regex=None,
):
    """Select columns from a data matrix.

    Parameters
    ----------
    data : array-like, shape=[n_samples, n_features]
        Input data
    extra_data : array-like, shape=[any, n_features], optional
        Optional additional data objects from which to select the same rows
    idx : list-like, shape=[m_features]
        Integer indices or string column names to be selected
    starts_with : str, list-like or None, optional (default: None)
        If not None, select columns that start with this prefix.
    ends_with : str, list-like or None, optional (default: None)
        If not None, select columns that end with this suffix.
    exact_word : str, list-like or None, optional (default: None)
        If not None, select columns that contain this exact word.
    regex : str, list-like or None, optional (default: None)
        If not None, select columns that match this regular expression.

    Returns
    -------
    data : array-like, shape=[n_samples, m_features]
        Subsetted output data.
    extra_data : array-like, shape=[any, m_features]
        Subsetted extra data, if passed.

    Examples
    --------
    data_subset = scprep.select.select_cols(
        data,
        idx=np.random.choice([True, False],
        data.shape[1])
    )
    data_subset, metadata_subset = scprep.select.select_cols(
        data,
        metadata,
        starts_with="MT"
    )

    Raises
    ------
    UserWarning : if no columns are selected
    """
    if (
        idx is None
        and starts_with is None
        and ends_with is None
        and exact_word is None
        and regex is None
    ):
        warnings.warn(
            "No selection conditions provided. Returning all columns.", UserWarning
        )
    if idx is None:
        if not isinstance(data, pd.DataFrame):
            raise ValueError(
                "Can only select based on column names with DataFrame input. "
                "Please set `idx` to select specific columns."
            )

    if not isinstance(idx, (numbers.Integral, str)):
        idx = toarray(idx)
        _check_idx_1d(idx)
        idx = idx.flatten()

    if is_SparseDataFrame(data):
        # evil deprecated dataframe; get rid of it
        data = SparseDataFrame(data)

    input_1d = _is_1d(data)
    if isinstance(data, pd.DataFrame):
        try:
            if isinstance(idx, (numbers.Integral, str)):
                data = data.loc[:, idx]
                print(working)
            else:
                if np.issubdtype(idx.dtype, np.dtype(bool).type):
                   
                    raise TypeError
                data = data.loc[:, idx]
        except (KeyError, TypeError):
            if isinstance(idx, str):
                raise
            if (
                isinstance(idx, numbers.Integral)
                or np.issubdtype(idx.dtype, np.dtype(int))
                or np.issubdtype(idx.dtype, np.dtype(bool))
            ):

                data = data.loc[:, np.array(data.columns)[idx]]
            else:
                raise
    elif isinstance(data, pd.Series):
        try:
            if np.issubdtype(idx.dtype, np.dtype(bool).type):
                # temporary workaround for pandas error
                raise TypeError
            data = data.loc[idx]
        except (KeyError, TypeError):
            if (
                isinstance(idx, numbers.Integral)
                or np.issubdtype(idx.dtype, np.dtype(int))
                or np.issubdtype(idx.dtype, np.dtype(bool))
            ):
                data = data.loc[np.array(data.index)[idx]]
            else:
                raise
    elif _is_1d(data):
        if isinstance(data, list):
            # can't numpy index a list
            data = np.array(data)
        data = data[idx]
    else:
        if isinstance(
            data,
            (
                sparse.coo_matrix,
                sparse.bsr_matrix,
                sparse.lil_matrix,
                sparse.dia_matrix,
            ),
        ):
            data = data.tocsr()
        if isinstance(idx, pd.Series):
            idx = toarray(idx)
        data = data[:, idx]
    if _get_column_length(data) == 0:
        warnings.warn("Selecting 0 columns.", UserWarning)
    elif isinstance(data, pd.DataFrame) and not input_1d:
        # convert to series if possible
        data = _convert_dataframe_1d(data, silent=True)
    return data


def filter_rare_genes(data, cutoff=0, min_cells=5):
    """Filter all genes with negligible counts in all but a few cells.

    Parameters
    ----------
    data : array-like, shape=[n_samples, n_features]
        Input data
    extra_data : array-like, shape=[any, n_features], optional
        Optional additional data objects from which to select the same rows
    cutoff : float, optional (default: 0)
        Number of counts above which expression is deemed non-negligible
    min_cells : int, optional (default: 5)
        Minimum number of cells above `cutoff` in order to retain a gene

    Returns
    -------
    data : array-like, shape=[n_samples, m_features]
        Filtered output data, where m_features <= n_features
    extra_data : array-like, shape=[any, m_features]
        Filtered extra data, if passed.
    """
    gene_sums = gene_capture_count(data, cutoff=cutoff)
    #print("gene_sums ")  
    #print(gene_sums)
    keep_genes_idx = gene_sums >= min_cells
    data = select_cols(data, idx=keep_genes_idx)
    return data


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

def matrix_sum(data, axis=None):
    """Get the column-wise, row-wise, or total sum of values in a matrix.

    Parameters
    ----------
    data : array-like, shape=[n_samples, n_features]
        Input data
    axis : int or None, optional (default: None)
        Axis across which to sum. axis=0 gives column sums,
        axis=1 gives row sums. None gives the total sum.

    Returns
    -------
    sums : array-like or float
        Sums along desired axis.
    """
    if axis not in [0, 1, None]:
        raise ValueError("Expected axis in [0, 1, None]. Got {}".format(axis))
    if isinstance(data, pd.DataFrame):
        if is_SparseDataFrame(data):
            if axis is None:
                sums = data.to_coo().sum()
            else:
                index = data.index if axis == 1 else data.columns
                sums = pd.Series(
                    np.array(data.to_coo().sum(axis)).flatten(), index=index
                )
        elif is_sparse_dataframe(data):
            if axis is None:
                sums = data.sparse.to_coo().sum()
            else:
                index = data.index if axis == 1 else data.columns
                sums = pd.Series(
                    np.array(data.sparse.to_coo().sum(axis)).flatten(), index=index
                )
        elif axis is None:
            sums = data.to_numpy().sum()
        else:
            sums = data.sum(axis)
    else:
        sums = np.sum(data, axis=axis)
        if isinstance(sums, np.matrix):
            sums = np.array(sums).flatten()
    return sums

def _convert_dataframe_1d(idx, silent=False):
    if _check_idx_1d(idx, silent=silent):
        idx = idx.iloc[:, 0] if idx.shape[1] == 1 else idx.iloc[0, :]
    return idx

def _check_idx_1d(idx, silent=False):
    if (not _is_1d(idx)) and np.prod(idx.shape) != np.max(idx.shape):
        if silent:
            return False
        else:
            raise ValueError("Expected idx to be 1D. Got shape {}".format(idx.shape))
    else:
        return True



def _is_1d(data):
    print("_is_1d")
    try:
        return len(data.shape) == 1
    except AttributeError:
        return True



data = pd.read_csv(r"C:\Users\maky\Desktop\SRNAPY\data\test_10X\output\data.csv")  

index_data= data["0"]
del data['0']

data = data.set_index([pd.Index(index_data)])

data=filter_rare_genes(data, cutoff=0, min_cells=5)

print(data)
data.to_csv(r"C:\Users\maky\Desktop\SRNAPY\data\test_10X\output\data.csv")
