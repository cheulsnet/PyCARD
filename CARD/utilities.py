####################################################################################################
## Package : CARD with Python
## Version : 1.0.1
## Date    : 2024-5-12 11:30:00
## Modified:
## Title   : Spatially Informed Cell Type Deconvolution for Spatial Transcriptomics by CARD.
## Authors : Ying Ma, Faith Shim
## Contacts: ying_ma@brown.edu & faith_shim@brown.edu
##           Brown University, Department of Biostatistics & Computer Science Concentrator
####################################################################################################

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix


class CARD:
    """
    Each CARD object has a number of slots which store information. Key slots to access are listed below.
    """

    def __init__(self):
        self.sc_eset = None
        self.spatial_countMat = None
        self.spatial_location = pd.DataFrame()
        self.Proportion_CARD = np.array([])
        self.project = "deconvolution"
        self.info_parameters = {}
        self.algorithm_matrix = {}
        self.refined_prop = np.array([])
        self.refined_expression = np.array([])


def sc_QC(counts_in, metaData, ct_varname, ct_select, sc_gene, min_cells=0, min_genes=0):
    """
    Quality control of scRNA-seq count data.
    """
    coldf = metaData.copy()
    counts = counts_in.copy()
    sc_gene_index = np.asarray(sc_gene)

    # Filter cells by total counts
    if min_genes >= 0:
        nfeatures = counts.sum(axis=0).A1
        keep_cells = nfeatures > min_genes
        counts = counts[:, keep_cells]
        coldf = coldf.loc[keep_cells, :]

    # Filter genes by total counts
    if min_cells >= 0:
        row_sum = counts.sum(axis=1).A1
        keep_rows = row_sum > min_cells
        counts = counts[keep_rows, :]
        sc_gene_index = sc_gene_index[keep_rows]

    row_indices_csr, _ = counts.nonzero()
    unique_row_indices = np.unique(row_indices_csr)
    fdata = pd.DataFrame({"Gene": unique_row_indices}).astype(str)

    keep_cell = coldf[ct_varname].astype(str).isin(ct_select)
    counts = counts[:, keep_cell]
    coldf = coldf.loc[keep_cell, :]

    keep_gene = counts.sum(axis=1).A1 > 0
    fdata = fdata[fdata.index.isin(np.flatnonzero(keep_gene))].reset_index(drop=True)
    counts = counts[keep_gene, :].tocsc()
    sc_gene_index = sc_gene_index[keep_gene]

    if coldf.index.dtype != "object":
        coldf.index = coldf.index.astype(str)
    if fdata.index.dtype != "object":
        fdata.index = fdata.index.astype(str)

    adata = AnnData(X=counts.transpose(), obs=coldf, var=fdata)
    return adata, sc_gene_index


def create_CARDObject(
    sc_count,
    sc_meta,
    spatial_count,
    spatial_location,
    ct_varname,
    ct_select,
    sample_varname=None,
    minCountGene=100,
    minCountSpot=5,
):
    """
    Create the CARD object.
    """
    print("## QC on scRNASeq dataset! ...")

    sc_gene = np.asarray(sc_count.index)
    sc_countMat = _convert_to_sparse(sc_count)

    _check_scRNASeq_consistency(sc_countMat, sc_meta, sample_varname)
    ct_select = _check_cell_type_info(ct_varname, ct_select, sc_meta)

    sc_eset, sc_gene = sc_QC(sc_countMat, sc_meta, ct_varname, ct_select, sc_gene)

    spatial_gene = np.asarray(spatial_count.index)
    spatial_countMat = _convert_to_sparse(spatial_count)

    _check_spatial_dataset_consistency(spatial_countMat, spatial_location, spatial_gene, sc_gene)

    spatial_countMat, spatial_gene, spatial_location = _filter_spatial_counts(
        spatial_countMat,
        spatial_gene,
        spatial_location,
        minCountSpot,
        minCountGene,
    )

    spatial_countMat.index = spatial_gene
    spatial_countMat.index.name = "gene"

    obj = AnnData(X=sc_eset.X, obs=sc_eset.obs, var=sc_eset.var)
    _update_uns_attributes(
        obj,
        sc_eset,
        spatial_countMat,
        spatial_location,
        ct_varname,
        ct_select,
        sample_varname,
        sc_gene,
    )
    return obj


def _convert_to_sparse(data):
    if isinstance(data, pd.DataFrame):
        return csr_matrix(data.to_numpy(copy=False))
    elif isinstance(data, pd.Series):
        return csr_matrix(data.to_numpy(copy=False).reshape(1, -1))
    elif isinstance(data, csr_matrix):
        return data
    else:
        raise ValueError("Input has to be of following forms: DataFrame, Series, or sparseMatrix")


def _check_scRNASeq_consistency(sc_countMat, sc_meta, sample_varname):
    if sc_countMat is None:
        raise ValueError("Please provide scRNASeq count data")

    if sample_varname is None:
        sc_meta = pd.DataFrame(sc_meta)
        if "sampleID" not in sc_meta.columns:
            sc_meta["sampleID"] = "Sample"

    if sc_countMat.shape[1] != sc_meta.shape[0]:
        raise ValueError(
            "The number of cells in scRNA-seq counts and sc_meta should be consistent! "
            "(sc_count -- p x c; sc_meta -- c x 2)"
        )


def _check_cell_type_info(ct_varname, ct_select, sc_meta):
    if ct_varname is None:
        raise ValueError(
            "Please provide the column name indicating the cell type information in the meta data of scRNA-seq"
        )

    if ct_select is None:
        print("No cell types selected, we will use all the cell types in the scRNA-seq data")
        ct_select = sc_meta[ct_varname].unique()

    return [str(ct) for ct in ct_select if not pd.isna(ct)]


def _check_spatial_dataset_consistency(spatial_countMat, spatial_location, spatial_gene, sc_gene):
    if spatial_location is None:
        raise ValueError("Please provide the matched spatial location data frame")

    if spatial_countMat.shape[1] != spatial_location.shape[0]:
        raise ValueError(
            "The number of spatial locations in spatial_count and spatial_location should be consistent! "
            "(spatial_count -- p x n; spatial_location -- n x 2)"
        )

    common_genes = np.intersect1d(np.asarray(spatial_gene).astype(str), np.asarray(sc_gene).astype(str))
    if len(common_genes) == 0:
        raise ValueError("There are no common gene names in spatial count data and single cell RNAseq count data")

    if not np.array_equal(np.arange(spatial_countMat.shape[1]).astype(str), spatial_location.index.astype(str)):
        raise ValueError("The gene names in spatial_countMat do not match the index of spatial_location")


def _filter_spatial_counts(spatial_countMat, spatial_gene, spatial_location, minCountSpot, minCountGene):
    # Filter genes by row sum
    row_sum = spatial_countMat.sum(axis=1).A1
    valid_rows = row_sum > minCountSpot
    spatial_countMat = spatial_countMat[valid_rows, :]
    spatial_gene = np.asarray(spatial_gene)[valid_rows]

    # Filter spots by column sum
    col_sums = np.asarray(spatial_countMat.sum(axis=0)).ravel()
    valid_cols = (col_sums >= minCountGene) & (col_sums <= 1e6)
    spatial_countMat = spatial_countMat[:, valid_cols]
    spatial_location = spatial_location.loc[spatial_location.index.isin(np.flatnonzero(valid_cols)), :]

    # Build dense DataFrame once for downstream code
    spatial_count_dense = spatial_countMat.toarray()
    spatial_count_df = pd.DataFrame(
        spatial_count_dense,
        index=spatial_gene,
        columns=np.flatnonzero(valid_cols),
    )
    return spatial_count_df, spatial_gene, spatial_location


def _update_uns_attributes(
    obj,
    sc_eset,
    spatial_countMat,
    spatial_location,
    ct_varname,
    ct_select,
    sample_varname,
    sc_gene,
):
    obj.uns["Class"] = "CARD"
    obj.uns["sc_eset"] = sc_eset
    obj.uns["sc_gene"] = np.asarray(sc_gene)
    obj.uns["spatial_countMat"] = spatial_countMat
    obj.uns["spatial_location"] = spatial_location
    obj.uns["project"] = "Deconvolution"
    obj.uns["info_parameters"] = {
        "ct.varname": ct_varname,
        "ct.select": ct_select,
        "sample.varname": sample_varname,
    }


class CARDfree:
    """
    Each CARDfree object has a number of slots which store information. Key slots to access are listed below.
    """

    def __init__(self):
        self.spatial_countMat = None
        self.spatial_location = pd.DataFrame()
        self.Proportion_CARD = np.array([])
        self.estimated_refMatrix = np.array([])
        self.project = None
        self.markerList = []
        self.info_parameters = {}
        self.algorithm_matrix = {}
        self.refined_prop = np.array([])
        self.refined_expression = np.array([])


def _filter_spatial_counts_CARDfree(spatial_countMat, spatial_gene, spatial_location, minCountSpot, minCountGene):
    row_sums = spatial_countMat.getnnz(axis=1)
    selected_rows = row_sums > minCountSpot
    spatial_countMat = spatial_countMat[selected_rows]
    spatial_gene = np.asarray(spatial_gene)[selected_rows]

    col_sums = np.asarray(spatial_countMat.sum(axis=0)).ravel()
    selected_cols = (col_sums >= minCountGene) & (col_sums <= 1e6)
    spatial_countMat = spatial_countMat[:, selected_cols]
    spatial_location = spatial_location[selected_cols]

    dense_array = spatial_countMat.toarray()
    spatial_countMat = pd.DataFrame(dense_array, index=spatial_gene, columns=spatial_location.index)

    return spatial_countMat, spatial_gene, spatial_location


def create_CARDfreeObject(
    markerList,
    spatial_count,
    spatial_location,
    nmfSelect=1,
    minCountGene=100,
    minCountSpot=5,
    marker_ct=None,
):
    """
    Create the CARDfree object.
    """
    spatial_gene = np.asarray(spatial_count.index)

    if isinstance(spatial_count, np.ndarray):
        if spatial_count.ndim == 2:
            spatial_countMat = csr_matrix(spatial_count)
        elif spatial_count.ndim == 1:
            spatial_countMat = csr_matrix(spatial_count.reshape(1, -1))
        else:
            raise ValueError("Invalid spatial count format. Must be a vector or matrix.")
    elif isinstance(spatial_count, pd.DataFrame):
        spatial_countMat = csr_matrix(spatial_count.to_numpy(copy=False))
    elif isinstance(spatial_count, csr_matrix):
        spatial_countMat = spatial_count
    else:
        raise ValueError(
            "Spatial resolved transcriptomic counts must be of following forms: vector, matrix, or sparseMatrix"
        )

    if spatial_location is None:
        raise ValueError("Please provide the matched spatial location data frame")

    if spatial_countMat.shape[1] != spatial_location.shape[0]:
        raise ValueError(
            "The number of spatial locations in spatial_count and spatial_location should be consistent! "
            "(spatial_count -- p x n; spatial_location -- n x 2)"
        )

    spatial_countMat, spatial_gene, spatial_location = _filter_spatial_counts_CARDfree(
        spatial_countMat,
        spatial_gene,
        spatial_location,
        minCountSpot,
        minCountGene,
    )

    spatial_countMat.index = spatial_gene
    spatial_countMat.index.name = "gene"

    if marker_ct is None:
        n = len(markerList)
        marker_ct = [f"CT{i+1}" for i in range(n)]

    spatial_location.index = spatial_location.index.astype(str)
    card_free_object = AnnData(X=spatial_countMat.T, obs=spatial_location)

    card_free_object.uns["Class"] = "CARDfree"
    card_free_object.uns["markerList"] = markerList
    card_free_object.uns["marker_ct"] = marker_ct
    card_free_object.uns["spatial_gene"] = spatial_gene
    card_free_object.uns["nmfSelect"] = nmfSelect
    card_free_object.uns["project"] = "Deconvolution (reference-free)"
    card_free_object.uns["info_parameters"] = {}

    return card_free_object