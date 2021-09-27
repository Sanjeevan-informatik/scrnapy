import utils

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


def _find_gz_file(*path):
    """Find a file that could be gzipped."""
    path = os.path.join(*path)
    if os.path.isfile(path):
        return path
    else:
        return path + ".gz"


def load_10X(data_dir, sparse=True):
    """Load data produced from the 10X Cellranger pipeline.

    A default run of the `cellranger count` command will generate gene-barcode
    matrices for secondary analysis. For both "raw" and "filtered" output,
    directories are created containing three files:
    'matrix.mtx', 'barcodes.tsv', 'genes.tsv'.
    Running `scprep.io.load_10X(data_dir)` will return a Pandas DataFrame with
    genes as columns and cells as rows.

    Parameters
    ----------
    data_dir: string
        path to input data directory
        expects 'matrix.mtx(.gz)', '[genes/features].tsv(.gz)', 'barcodes.tsv(.gz)'
        to be present and will raise an error otherwise
    sparse: boolean
        If True, a sparse Pandas DataFrame is returned.
    gene_labels: string, {'id', 'symbol', 'both'} optional, default: 'symbol'
        Whether the columns of the dataframe should contain gene ids or gene
        symbols. If 'both', returns symbols followed by ids in parentheses.
    allow_duplicates : bool, optional (default: None)
        Whether or not to allow duplicate gene names. If None, duplicates are
        allowed for dense input but not for sparse input.

    Returns
    -------
    data: array-like, shape=[n_samples, n_features]
        If sparse, data will be a pd.DataFrame[pd.SparseArray]. Otherwise, data will
        be a pd.DataFrame.
    """
    gene_labels="symbol"
    
    if gene_labels not in ["id", "symbol", "both"]:
        raise ValueError(
            "gene_labels='{}' not recognized. "
            "Choose from ['symbol', 'id', 'both']".format(gene_labels)
        )

    if not os.path.isdir(data_dir):
        raise FileNotFoundError("{} is not a directory".format(data_dir))

    try:
        m = sio.mmread(_find_gz_file(data_dir, "matrix.mtx"))
        try:
            genes = pd.read_csv(
                _find_gz_file(data_dir, "genes.tsv"), delimiter="\t", header=None
            )
        except FileNotFoundError:
            genes = pd.read_csv(
                _find_gz_file(data_dir, "features.tsv"), delimiter="\t", header=None
            )
        if genes.shape[1] == 2:
            # Cellranger < 3.0
            genes.columns = ["id", "symbol"]
      
        barcodes = pd.read_csv(
            _find_gz_file(data_dir, "barcodes.tsv"), delimiter="\t", header=None
        )

    except (FileNotFoundError, IOError):
        raise FileNotFoundError(
            "'matrix , genes/features and barcodes file is not found' {}".format(data_dir)
        )

    cell_names = barcodes[0]
    gene_names = _parse_10x_genes(
        genes["symbol"].values.astype(str),
        genes["id"].values.astype(str),
        gene_labels=gene_labels,
    )

    data = _matrix_to_data_frame(
        m.T, cell_names=cell_names, gene_names=gene_names, sparse=sparse
    )
    return data

def _parse_10x_genes(symbols, ids, gene_labels="symbol"):
    assert gene_labels in ["symbol", "id", "both"]
    if gene_labels == "symbol":
        columns = symbols
        if len(np.unique(columns)) < len(columns):
            warnings.warn(
                "Duplicate gene names detected!",
                RuntimeWarning,
            )
            gene_labels = "both"

    if gene_labels == "both":
        columns = _combine_gene_id(symbols, ids)
    elif gene_labels == "id":
        columns = ids
    return columns


def _matrix_to_data_frame(data, gene_names=None, cell_names=None, sparse=None):
    """Return the optimal data type given data, gene names and cell names.

    Parameters
    ----------
    data : array-like
    gene_names : `str`, array-like or `None` (default: None)
        Either a filename or an array containing a list of gene symbols or ids.
    cell_names : `str`, array-like or `None` (default: None)
        Either a filename or an array containing a list of cell barcodes.
    sparse : `bool` or `None` (default: None)
        If not `None`, overrides default sparsity of the data.
    """
     # print(data)
     # print(gene_names)  
     # print(cell_names)  
    if gene_names is None and cell_names is None and not isinstance(data, pd.DataFrame):
        # just a matrix
            warnings.warn(
                "the data is not correct",
                RuntimeWarning,
            )
    else:
        
        if sparse:
            # return pandas.DataFrame[SparseArray]
             #print(data)
             #print(pd.DataFrame)
            
            if sp.issparse(data):
                print("issparse")
                data = pd.DataFrame.sparse.from_spmatrix(
                    data, index=cell_names, columns=gene_names
                )

       
    #print(data)
    data =  utils.check_numeric(data, suppress_errors=True)
    #print(data)
    return data


data_uf = load_10X(r"C:\Users\maky\Desktop\SRNAPY\data\test_10X", sparse=True)

print(data_uf)

data_uf.to_csv(r"C:\Users\maky\Desktop\SRNAPY\data\test_10X\output\data.csv")