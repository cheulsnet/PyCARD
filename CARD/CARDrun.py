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
import matplotlib.pyplot as plt
import pandas as pd

from CARD.utilities import create_CARDObject
from CARD.utilities import create_CARDfreeObject
from CARD.CARDprop import CARD_deconvolution
from CARD.CARDimputation import CARD_imputation
from CARD.CARDrefFree import CARD_refFree
from CARD.CARDSCMapping import CARD_SCMapping

from CARD.visualization import CARD_visualize_pie
from CARD.visualization import CARD_visualize_prop
from CARD.visualization import CARD_visualize_prop_2CT
from CARD.visualization import CARD_visualize_Cor
from CARD.visualization import plot_location_imputation
from CARD.visualization import CARD_visualize_gene
from CARD.visualization import CARD_visualize_SCMapping


def _run_CARD_deconvolution(
    sc_count,
    sc_meta,
    spatial_count,
    spatial_location,
    phi_grid=None,
    warm_start_phi=False,
    max_iter=1000,
    epsilon=1e-04,
    verbose=True,
):
    CARD_obj = create_CARDObject(
        sc_count=sc_count,
        sc_meta=sc_meta,
        spatial_count=spatial_count,
        spatial_location=spatial_location,
        ct_varname="cellType",
        ct_select=sc_meta["cellType"].unique(),
        sample_varname="sampleInfo",
        minCountGene=100,
        minCountSpot=5,
    )

    CARD_obj = CARD_deconvolution(
        CARD_object=CARD_obj,
        phi_grid=phi_grid,
        warm_start_phi=warm_start_phi,
        max_iter=max_iter,
        epsilon=epsilon,
    )

    if verbose:
        print(CARD_obj.uns["info_parameters"]["Proportion_CARD"].iloc[0:2])

        profiling = CARD_obj.uns.get("profiling", {}).get("card_deconvolution_internal", {})
        if profiling:
            print("## Internal profiling summary")
            print({
                "createscRef_seconds": profiling.get("createscRef_seconds"),
                "selectInfo_seconds": profiling.get("selectInfo_seconds"),
                "kernel_build_seconds": profiling.get("kernel_build_seconds"),
                "phi_total_seconds": profiling.get("phi_total_seconds"),
                "winning_phi": profiling.get("winning_phi"),
                "warm_start_phi": profiling.get("warm_start_phi"),
                "max_iter": profiling.get("max_iter"),
            })

    return CARD_obj


def _run_CARD_imputation(CARD_obj, num_grids=2000, in_neighbor=10, exclude=None, showGrid=False):
    CARD_obj = CARD_imputation(
        CARD_object=CARD_obj,
        num_grids=num_grids,
        in_neighbor=in_neighbor,
        exclude=exclude
    )

    if showGrid is True:
        _ = plot_location_imputation(
            CARD_obj.uns["spatial_location"],
            CARD_obj.uns["info_parameters"]["location_imputation"]
        )
        plt.show()

    return CARD_obj


def _run_CARDfree_deconvolution(
    spatial_count,
    spatial_location,
    nmfSelect,
    marker_list_path,
):
    markerList_df = pd.read_csv(marker_list_path)
    markerList = [markerList_df[column].tolist() for column in markerList_df.columns]

    CARDfree_obj = create_CARDfreeObject(
        markerList=markerList,
        spatial_count=spatial_count,
        spatial_location=spatial_location,
        nmfSelect=nmfSelect,
        minCountGene=100,
        minCountSpot=5,
        marker_ct=None
    )

    CARDfree_obj = CARD_refFree(CARDfree_obj)
    print(CARDfree_obj.uns["info_parameters"]["Proportion_CARD"].iloc[0:2])

    return CARDfree_obj


def _run_CARD_scMapping(CARD_obj, shapeSpot, numCell, ncore):
    scMapping = CARD_SCMapping(
        CARD_object=CARD_obj,
        shapeSpot=shapeSpot,
        numCell=numCell,
        ncore=ncore
    )

    print(scMapping)
    return scMapping


def _run_CARD_visualization(vis_functions):
    visual_functions = {
        "visual_pie": CARD_visualize_pie,
        "visual_prop": CARD_visualize_prop,
        "visual_prop_2CT": CARD_visualize_prop_2CT,
        "visual_Cor": CARD_visualize_Cor,
        "visual_imput_loc": plot_location_imputation,
        "visual_gene": CARD_visualize_gene,
        "visual_scmapping": CARD_visualize_SCMapping
    }

    for function_name, function_args in vis_functions:
        if function_name in visual_functions:
            _ = visual_functions[function_name](**function_args)
            plt.show()
        else:
            print(f"Warning: '{function_name}' is not a valid visualization function.")