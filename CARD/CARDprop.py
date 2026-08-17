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

import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import dirichlet

build_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../build"))
sys.path.insert(0, build_dir)

from CARDref_module import CARDref


DEFAULT_PHI_GRID = [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99]


def _to_dense_float64(x):
    if hasattr(x, "toarray"):
        arr = x.toarray()
    else:
        arr = np.asarray(x)
    return np.asarray(arr, dtype=np.float64)


def _normalize_phi_grid(phi_grid):
    if phi_grid is None:
        return list(DEFAULT_PHI_GRID)
    if isinstance(phi_grid, (list, tuple, np.ndarray)):
        vals = [float(x) for x in phi_grid]
    else:
        raise ValueError("phi_grid must be a list, tuple, or numpy array of floats")
    if len(vals) == 0:
        raise ValueError("phi_grid cannot be empty")
    return vals


def row_grp_means_by_chunk(merged_df, chunk_size):
    mean_chunks = []
    for _, group_data in merged_df.groupby("ct_sample_id"):
        chunks = [group_data.iloc[i : i + chunk_size] for i in range(0, len(group_data), chunk_size)]
        mean_chunks.append(pd.concat(chunks).groupby("ct_sample_id").mean())
    return pd.concat(mean_chunks)


def row_grp_means(df_data):
    return df_data.groupby("ct_sample_id").mean()


def createscRef(sc_eset, ct_varname, ct_select=None, sample_varname=None):
    X = _to_dense_float64(sc_eset.X)  # cells x genes

    if ct_select is None:
        ct_select = pd.unique(sc_eset.obs[ct_varname]).tolist()

    ct_select = [ct for ct in ct_select if pd.notna(ct)]
    ct_select = pd.Series(ct_select).unique()

    ct_arr = sc_eset.obs[ct_varname].astype(str).to_numpy()

    if sample_varname is None:
        sample_arr = np.full(ct_arr.shape[0], "Sample", dtype=object)
    else:
        sample_arr = sc_eset.obs[sample_varname].astype(str).to_numpy()

    ct_sample_id = np.char.add(np.char.add(ct_arr.astype(str), "$*$"), sample_arr.astype(str))
    unique_groups, inv = np.unique(ct_sample_id, return_inverse=True)

    n_groups = unique_groups.shape[0]
    n_genes = X.shape[1]

    group_sums = np.zeros((n_groups, n_genes), dtype=np.float64)
    np.add.at(group_sums, inv, X)

    group_counts = np.bincount(inv).astype(np.float64)
    group_means = group_sums / group_counts[:, None]

    group_ct = np.array([g.split("$*$", 1)[0] for g in unique_groups], dtype=object)

    basis = pd.DataFrame(group_means, index=group_ct)
    return {"basis": basis}


def selectInfo(Basis, sc_eset, sc_gene, ct_select, ct_varname):
    if Basis.empty:
        return []

    basis_index = Basis.index.to_numpy()
    basis_cols = np.asarray(Basis.columns.to_list(), dtype=object)
    B = Basis.to_numpy(dtype=np.float64, copy=False)

    unique_basis_cols = set(basis_cols.tolist())
    ct_select_clean = [ct for ct in ct_select if ct in unique_basis_cols]

    gene_keep_mask = np.zeros(B.shape[0], dtype=bool)

    for ict in ct_select_clean:
        target_mask = basis_cols == ict
        rest_mask = ~target_mask

        if not np.any(target_mask):
            continue

        target_vals = B[:, target_mask]
        target_mean = target_vals if target_vals.ndim == 1 else target_vals.mean(axis=1)

        if np.any(rest_mask):
            rest_mean = B[:, rest_mask].mean(axis=1)
        else:
            rest_mean = np.zeros(B.shape[0], dtype=np.float64)

        fc = np.log(target_mean + 1e-6) - np.log(rest_mean + 1e-6)
        gene_keep_mask |= (fc > 1.25) & (target_mean > 0)

    gene1 = basis_index[gene_keep_mask]
    if gene1.size == 0:
        return []

    counts = _to_dense_float64(sc_eset.X).T  # genes x cells
    sc_gene_arr = np.asarray(sc_gene)

    gene_to_idx = {g: i for i, g in enumerate(sc_gene_arr)}
    gene1_idx = np.array([gene_to_idx[g] for g in gene1 if g in gene_to_idx], dtype=int)

    if gene1_idx.size == 0:
        return []

    counts_sel = counts[gene1_idx, :]

    cell_types = sc_eset.obs[ct_varname].astype(str).to_numpy()
    cell_type_counts = pd.Series(cell_types).value_counts()
    ct_valid = cell_type_counts.index[cell_type_counts > 1].tolist()

    if len(ct_valid) == 0:
        return gene1.tolist()

    sd_cols = []
    for ict in ct_valid:
        mask = cell_types == ict
        if mask.sum() <= 1:
            continue

        temp = counts_sel[:, mask]
        mean = temp.mean(axis=1)
        var = temp.var(axis=1, ddof=1)

        with np.errstate(divide="ignore", invalid="ignore"):
            sd = np.divide(
                var,
                mean,
                out=np.full_like(var, np.nan, dtype=np.float64),
                where=mean != 0,
            )
        sd_cols.append(sd)

    if not sd_cols:
        return gene1.tolist()

    sd_within_mat = np.column_stack(sd_cols)
    means = np.nanmean(sd_within_mat, axis=1)
    threshold = np.nanpercentile(means, 99)

    condition = np.isnan(means) | (means < threshold)
    gene2 = gene1[condition]
    return gene2.tolist()


def CARD_deconvolution(
    CARD_object,
    phi_grid=None,
    warm_start_phi=False,
    max_iter=1000,
    epsilon=1e-04,
):
    profiling = {
        "createscRef_seconds": 0.0,
        "selectInfo_seconds": 0.0,
        "kernel_build_seconds": 0.0,
        "phi_seconds": {},
        "phi_objectives": {},
        "phi_total_seconds": 0.0,
        "postprocess_seconds": 0.0,
        "total_internal_seconds": 0.0,
        "warm_start_phi": bool(warm_start_phi),
        "max_iter": int(max_iter),
        "phi_grid": [],
        "winning_phi": None,
        "winning_obj": None,
    }
    t_total0 = time.perf_counter()

    phi_grid = _normalize_phi_grid(phi_grid)
    profiling["phi_grid"] = list(phi_grid)

    info = CARD_object.uns["info_parameters"]
    ct_select = info["ct.select"]
    ct_varname = info["ct.varname"]
    sample_varname = info["sample.varname"]

    sc_eset = CARD_object.uns["sc_eset"]
    sc_gene = np.asarray(CARD_object.uns["sc_gene"])

    print("## create reference matrix from scRNASeq...")

    t0 = time.perf_counter()
    Basis_ref = createscRef(sc_eset, ct_varname, ct_select, sample_varname)
    profiling["createscRef_seconds"] = time.perf_counter() - t0

    Basis = Basis_ref["basis"].T
    Basis.index = sc_gene
    Basis = Basis.loc[:, Basis.columns.isin(ct_select)]

    spatial_count = CARD_object.uns["spatial_countMat"]
    spatial_count = spatial_count.loc[spatial_count.sum(axis=1) > 0, :]
    spatial_count = spatial_count.loc[:, spatial_count.sum(axis=0) > 0]

    common_gene = spatial_count.index.intersection(Basis.index)
    Basis = Basis.loc[common_gene]

    print("## Select Informative Genes! ...")
    t0 = time.perf_counter()
    selected_common = selectInfo(Basis, sc_eset, sc_gene, ct_select, ct_varname)
    profiling["selectInfo_seconds"] = time.perf_counter() - t0

    if len(selected_common) == 0:
        raise ValueError("No informative genes were selected for deconvolution.")

    Basis = Basis.loc[selected_common]
    spatial_count = spatial_count.loc[selected_common]

    col_sum_vec = spatial_count.sum(axis=0)
    spatial_count = spatial_count.divide(col_sum_vec, axis=1)

    spatial_location = CARD_object.uns["spatial_location"]
    spatial_location = spatial_location.loc[spatial_location.index.isin(spatial_count.columns), :]
    spatial_location = spatial_location.reindex(spatial_count.columns)

    norm_cords = spatial_location[["x", "y"]].copy().astype(np.float64)
    norm_cords["x"] = norm_cords["x"] - norm_cords["x"].min()
    norm_cords["y"] = norm_cords["y"] - norm_cords["y"].min()
    scale_factor = max(norm_cords["x"].max(), norm_cords["y"].max())
    if scale_factor != 0:
        norm_cords["x"] = norm_cords["x"] / scale_factor
        norm_cords["y"] = norm_cords["y"] / scale_factor

    print("## Deconvolution Starts! ...")

    np.random.seed(20200107)
    Vint1 = dirichlet.rvs(alpha=[10] * Basis.shape[1], size=spatial_count.shape[1])
    spot_names = (spatial_location["x"].astype(str) + "x" + spatial_location["y"].astype(str)).to_numpy()

    isigma = 0.1

    t0 = time.perf_counter()
    coords = norm_cords.to_numpy(dtype=np.float64, copy=False)
    ED = cdist(coords, coords, metric="euclidean")
    kernel_mat = np.exp(-(ED ** 2) / (2 * isigma ** 2))
    np.fill_diagonal(kernel_mat, 0.0)
    profiling["kernel_build_seconds"] = time.perf_counter() - t0

    X = spatial_count.to_numpy(dtype=np.float64, copy=True)
    U = Basis.to_numpy(dtype=np.float64, copy=True)

    mean_X = X.mean()
    mean_B = U.mean()

    if mean_X != 0:
        X = X * (1e-01 / mean_X)
    if mean_B != 0:
        U = U * (1e-01 / mean_B)

    initV = np.asarray(Vint1, dtype=np.float64)
    b0 = np.zeros(U.shape[1], dtype=np.float64)
    lam0 = np.repeat(10.0, len(ct_select)).astype(np.float64)

    best_obj = -np.inf
    best_phi = None
    best_res = None

    phi_total_t0 = time.perf_counter()
    for phi_value in phi_grid:
        phi_t0 = time.perf_counter()
        res = CARDref(
            X,
            U,
            kernel_mat,
            float(phi_value),
            int(max_iter),
            float(epsilon),
            initV,
            b0,
            0.1,
            lam0,
        )
        phi_elapsed = time.perf_counter() - phi_t0

        phi_key = str(phi_value)
        profiling["phi_seconds"][phi_key] = phi_elapsed
        profiling["phi_objectives"][phi_key] = float(res["Obj"])

        current_obj = float(res["Obj"])
        if current_obj >= best_obj:
            best_obj = current_obj
            best_phi = float(phi_value)
            best_res = res

        if warm_start_phi:
            initV = np.asarray(res["V"], dtype=np.float64)

    profiling["phi_total_seconds"] = time.perf_counter() - phi_total_t0

    if best_res is None:
        raise RuntimeError("CARD deconvolution failed to produce a valid result.")

    profiling["winning_phi"] = best_phi
    profiling["winning_obj"] = best_obj

    t0 = time.perf_counter()
    best_V = pd.DataFrame(best_res["V"], index=spot_names, columns=Basis.columns)

    print("## Deconvolution Finish! ...")

    CARD_object.uns["info_parameters"]["phi"] = best_phi
    CARD_object.uns["info_parameters"]["phi_grid"] = list(phi_grid)
    CARD_object.uns["info_parameters"]["warm_start_phi"] = bool(warm_start_phi)
    CARD_object.uns["info_parameters"]["max_iter"] = int(max_iter)
    CARD_object.uns["info_parameters"]["epsilon"] = float(epsilon)
    CARD_object.uns["info_parameters"]["Proportion_CARD"] = best_V.div(best_V.sum(axis=1), axis=0)

    optimal_B = U * (mean_B / 1e-01) if mean_B != 0 else U
    optimal_B = pd.DataFrame(optimal_B, index=Basis.index, columns=Basis.columns).sort_index()

    xinput_norm = X * (mean_X / 1e-01) if mean_X != 0 else X
    xinput_norm = pd.DataFrame(xinput_norm, index=spatial_count.index, columns=spot_names)

    best_res_out = dict(best_res)
    best_res_out["V"] = best_V

    CARD_object.uns["info_parameters"]["algorithm_matrix"] = {
        "B": optimal_B,
        "Xinput_norm": xinput_norm,
        "Res": best_res_out,
    }

    CARD_object.uns["spatial_countMat"] = xinput_norm
    spatial_location_out = spatial_location.copy()
    spatial_location_out.index = spot_names
    CARD_object.uns["spatial_location"] = spatial_location_out

    profiling["postprocess_seconds"] = time.perf_counter() - t0
    profiling["total_internal_seconds"] = time.perf_counter() - t_total0

    if "profiling" not in CARD_object.uns:
        CARD_object.uns["profiling"] = {}
    CARD_object.uns["profiling"]["card_deconvolution_internal"] = profiling

    return CARD_object