from scipy import sparse

import numbers
import numpy as np
import pandas as pd
import warnings


def filter_library_size(
    data,
    cutoff=None,
    percentile=None,
    keep_cells=None,
  
):
    """Remove all cells with library size above or below a certain threshold.

    It is recommended to use :func:`~scprep.plot.plot_library_size` to
    choose a cutoff prior to filtering.

    Parameters
    ----------
    data : array-like, shape=[n_samples, n_features]
        Input data
    extra_data : array-like, shape=[n_samples, any], optional
        Optional additional data objects from which to select the same rows
    cutoff : float or tuple of floats, optional (default: None)
        Library size above or below which to retain a cell. Only one of `cutoff`
        and `percentile` should be specified.
    percentile : int or tuple of ints, optional (Default: None)
        Percentile above or below which to retain a cell.
        Must be an integer between 0 and 100. Only one of `cutoff`
        and `percentile` should be specified.
    keep_cells : {'above', 'below', 'between'} or None, optional (default: None)
        Keep cells above, below or between the cutoff.
        If None, defaults to 'above' when a single cutoff is given and
        'between' when two cutoffs are given.
    return_library_size : bool, optional (default: False)
        If True, also return the library sizes corresponding to the retained cells
    sample_labels : Deprecated
    filter_per_sample : Deprecated

    Returns
    -------
    data : array-like, shape=[m_samples, n_features]
        Filtered output data, where m_samples <= n_samples
    filtered_library_size : list-like, shape=[m_samples]
        Library sizes corresponding to retained samples,
        returned only if return_library_size is True
    extra_data : array-like, shape=[m_samples, any]
        Filtered extra data, if passed.
    """
    

    cell_sums = library_size(data)
    
    return filter_values(
        data,
        values=cell_sums,
        cutoff=cutoff,
        percentile=percentile,
        keep_cells=keep_cells,
        
    )


def library_size(data):
    """Measure the library size of each cell.

    Parameters
    ----------
    data : array-like, shape=[n_samples, n_features]
        Input data

    Returns
    -------
    library_size : list-like, shape=[n_samples]
        Sum over all genes for each cell
    """
    library_size = matrix_sum(data, axis=1)

    if isinstance(library_size, pd.Series):
        library_size.name = "library_size"
    return library_size


def filter_values(
    data,
    values=None,
    cutoff=None,
    percentile=None,
    keep_cells="above",

):
    """Remove all cells with `values` above or below a certain threshold.

   

    Parameters
    ----------
    data : array-like, shape=[n_samples, n_features]
        Input data
    extra_data : array-like, shape=[n_samples, any], optional
        Optional additional data objects from which to select the same rows
    values : list-like, shape=[n_samples]
        Value upon which to filter
    cutoff : float or tuple of floats, optional (default: None)
        Value above or below which to retain cells. Only one of `cutoff`
        and `percentile` should be specified.
    percentile : int or tuple of ints, optional (Default: None)
        Percentile above or below which to retain cells.
        Must be an integer between 0 and 100. Only one of `cutoff`
        and `percentile` should be specified.
    keep_cells : {'above', 'below', 'between'} or None, optional (default: None)
        Keep cells above, below or between the cutoff.
        If None, defaults to 'above' when a single cutoff is given and
        'between' when two cutoffs are given.
    return_values : bool, optional (default: False)
        If True, also return the values corresponding to the retained cells
    sample_labels : Deprecated
    filter_per_sample : Deprecated

    Returns
    -------
    data : array-like, shape=[m_samples, n_features]
        Filtered output data, where m_samples <= n_samples
    filtered_values : list-like, shape=[m_samples]
        Values corresponding to retained samples,
        returned only if return_values is True
    extra_data : array-like, shape=[m_samples, any]
        Filtered extra data, if passed.
    """
  
    assert values is not None
    keep_cells_idx = _get_filter_idx(values, cutoff, percentile, keep_cells)

    #print(data)
    #print(extra_data)
    #print(keep_cells_idx)
    data = select_rows(data,  idx=keep_cells_idx)
    return data


def _get_filter_idx(values, cutoff, percentile, keep_cells):
    """Return a boolean array to index cells based on a filter.

    Parameters
    ----------
    values : list-like, shape=[n_samples]
        Value upon which to filter
    cutoff : float or tuple of floats, optional (default: None)
        Value above or below which to retain cells. Only one of `cutoff`
        and `percentile` should be specified.
    percentile : int or tuple of ints, optional (Default: None)
        Percentile above or below which to retain cells.
        Must be an integer between 0 and 100. Only one of `cutoff`
        and `percentile` should be specified.
    keep_cells : {'above', 'below', 'between'} or None, optional (default: None)
        Keep cells above, below or between the cutoff.
        If None, defaults to 'above' when a single cutoff is given and
        'between' when two cutoffs are given.

    Returns
    -------
    keep_cells_idx : list-like
        Boolean retention array
    """
    cutoff = _get_percentile_cutoff(values, cutoff, percentile, required=True)
    if keep_cells is None:
        if isinstance(cutoff, numbers.Number):
            keep_cells = "above"
        else:
            keep_cells = "between"
    if keep_cells == "above":
        if not isinstance(cutoff, numbers.Number):
            raise ValueError(
                "Expected a single cutoff with keep_cells='above'."
                " Got {}".format(cutoff)
            )
        keep_cells_idx = values > cutoff

    elif keep_cells == "below":
        if not isinstance(cutoff, numbers.Number):
            raise ValueError(
                "Expected a single cutoff with keep_cells='below'."
                " Got {}".format(cutoff)
            )
        #print("values ",values)
        #print("cutoff ",cutoff)
        keep_cells_idx = values < cutoff
    elif keep_cells == "between":
        if isinstance(cutoff, numbers.Number) or len(cutoff) != 2:
            raise ValueError(
                "Expected cutoff of length 2 with keep_cells='between'."
                " Got {}".format(cutoff)
            )
        keep_cells_idx = np.logical_and(
            values > np.min(cutoff), values < np.max(cutoff)
        )
    else:
        raise ValueError(
            "Expected `keep_cells` in ['above', 'below', 'between']. "
            "Got {}".format(keep_cells)
        )
    #print(len(keep_cells_idx[keep_cells_idx==True]))
    return keep_cells_idx



def _get_percentile_cutoff(data, cutoff=None, percentile=None, required=False):
    """Get a cutoff for a dataset.

    Parameters
    ----------
    data : array-like
    cutoff : float or None, optional (default: None)
        Absolute cutoff value. Only one of cutoff and percentile may be given
    percentile : float or None, optional (default: None)
        Percentile cutoff value between 0 and 100.
        Only one of cutoff and percentile may be given
    required : bool, optional (default: False)
        If True, one of cutoff and percentile must be given.

    Returns
    -------
    cutoff : float or None
        Absolute cutoff value. Can only be None if required is False and
        cutoff and percentile are both None.
    """
    if percentile is not None:
        if cutoff is not None:
            raise ValueError(
                "Only one of `cutoff` and `percentile` should be given."
                "Got cutoff={}, percentile={}".format(cutoff, percentile)
            )
        if not isinstance(percentile, numbers.Number):
            return [_get_percentile_cutoff(data, percentile=p) for p in percentile]
        if percentile < 1:
            warnings.warn(
                "`percentile` expects values between 0 and 100."
                "Got {}. Did you mean {}?".format(percentile, percentile * 100),
                UserWarning,
            )
        cutoff = np.percentile(np.array(data).reshape(-1), percentile)
    elif cutoff is None and required:
        raise ValueError("One of either `cutoff` or `percentile` must be given.")
    return cutoff



def select_rows(
    data,
    idx=None,
    starts_with=None,
    ends_with=None,
    exact_word=None,
    regex=None,
):
    """Select rows from a data matrix.

    Parameters
    ----------
    data : array-like, shape=[n_samples, n_features]
        Input data
    extra_data : array-like, shape=[n_samples, any], optional
        Optional additional data objects from which to select the same rows
    idx : list-like, shape=[m_samples], optional (default: None)
        Integer indices or string index names to be selected
    starts_with : str, list-like or None, optional (default: None)
        If not None, select rows that start with this prefix.
    ends_with : str, list-like or None, optional (default: None)
        If not None, select rows that end with this suffix.
    exact_word : str, list-like or None, optional (default: None)
        If not None, select rows that contain this exact word.
    regex : str, list-like or None, optional (default: None)
        If not None, select rows that match this regular expression.

    Returns
    -------
    data : array-like, shape=[m_samples, n_features]
        Subsetted output data
    extra_data : array-like, shape=[m_samples, any]
        Subsetted extra data, if passed.

    Examples
    --------
    data_subset = scprep.select.select_rows(
        data,
        idx=np.random.choice([True, False],
        data.shape[0])
    )
    data_subset, labels_subset = scprep.select.select_rows(
        data,
        labels,
        end_with="batch1"
    )

    Raises
    ------
    UserWarning : if no rows are selected
    """
    
    if (
        idx is None
        and starts_with is None
        and ends_with is None
        and exact_word is None
        and regex is None
    ):
        warnings.warn(
            "No selection conditions provided. " "Returning all rows.", UserWarning
        )

    if idx is None:
        if not isinstance(data, pd.DataFrame):
            raise ValueError(
                "Can only select based on row names with DataFrame input. "
                "Please set `idx` to select specific rows."
            )
    #print("idx ",idx)
    if isinstance(idx, pd.DataFrame):
        idx = _convert_dataframe_1d(idx)
    elif not isinstance(idx, (numbers.Integral, str)):
        idx = toarray(idx)
        _check_idx_1d(idx)
        
        idx = idx.flatten()


    if isinstance(data, (pd.DataFrame, pd.Series)):
        
        if isinstance(idx, str):
            raise
        if (
            isinstance(idx, numbers.Integral)
            or np.issubdtype(idx.dtype, np.dtype(int))
            or np.issubdtype(idx.dtype, np.dtype(bool))
            ):
            data = data.loc[np.array(data.index)[idx]]
        else:
                raise
    return data


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
data = filter_library_size(data, cutoff=100, keep_cells='below')

print(data)

data.to_csv(r"C:\Users\maky\Desktop\SRNAPY\data\test_10X\output\data.csv")





