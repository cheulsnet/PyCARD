####################################################################################################
## Package : CARD with Python
## Version : 1.0.1
## Optimized main.py
##
## Practical Recommendation
## (Standard fast running)
##      python main.py
##
## (Standard + Extras)
##      python main.py --run-scmapping --run-imputation
##
## (CARDfree only)
##      python main.py --run-cardfree --set-r-libs-user ~/R/library  (if NFM is not installed)
##      python main.py --run-cardfree (if NFM is already installed)
##
####################################################################################################

import argparse
import os
import random

from CARD.dataload import _load_data
from CARD.CARDrun import _run_CARD_deconvolution
from CARD.CARDrun import _run_CARDfree_deconvolution
from CARD.CARDrun import _run_CARD_imputation
from CARD.CARDrun import _run_CARD_visualization
from CARD.CARDrun import _run_CARD_scMapping

### global random seed ###
random.seed(20200107)


def parse_phi_grid(phi_grid_str):
    if phi_grid_str is None or phi_grid_str.strip() == "":
        return None
    return [float(x.strip()) for x in phi_grid_str.split(",") if x.strip()]


def parse_args():
    parser = argparse.ArgumentParser(description="Run CARD pipeline.")

    parser.add_argument("--data-dir", default="./data", help="Directory containing input data files")
    parser.add_argument("--set-r-libs-user", default=None,
                        help="Optional R_LIBS_USER path for stable R package lookup, e.g. ~/R/library")

    parser.add_argument("--run-deconvolution", action="store_true", default=True,
                        help="Run standard CARD deconvolution")
    parser.add_argument("--run-scmapping", action="store_true",
                        help="Run CARD scMapping")
    parser.add_argument("--run-imputation", action="store_true",
                        help="Run CARD imputation")
    parser.add_argument("--run-cardfree", action="store_true",
                        help="Run reference-free CARDfree")

    parser.add_argument("--phi-grid", type=str, default="0.9",
                        help='Comma-separated phi grid for standard CARD, e.g. "0.9" or "0.1,0.5,0.9"')
    parser.add_argument("--warm-start-phi", action="store_true",
                        help="Warm-start across phi values in standard CARD")
    parser.add_argument("--max-iter", type=int, default=200,
                        help="Maximum iterations for standard CARD deconvolution")
    parser.add_argument("--epsilon", type=float, default=1e-4,
                        help="Convergence epsilon for standard CARD deconvolution")

    parser.add_argument("--scmapping-num-cell", type=int, default=20)
    parser.add_argument("--scmapping-ncore", type=int, default=10)

    parser.add_argument("--imputation-num-grids", type=int, default=2000)
    parser.add_argument("--imputation-in-neighbor", type=int, default=10)
    parser.add_argument("--imputation-show-grid", action="store_true")

    parser.add_argument("--cardfree-nmf-select", type=int, default=4,
                        help="NMF backend selection for CARDfree")
    parser.add_argument("--verbose", action="store_true",
                        help="Print extra profiling and result details")

    return parser.parse_args()


def main():
    args = parse_args()

    if args.set_r_libs_user:
        os.environ["R_LIBS_USER"] = os.path.expanduser(args.set_r_libs_user)

    phi_grid = parse_phi_grid(args.phi_grid)

    #### load data
    spatial_count = _load_data(f"{args.data_dir}/spatial_count.csv", "csv")
    spatial_location = _load_data(f"{args.data_dir}/spatial_location.csv", "csv")
    sc_count = _load_data(f"{args.data_dir}/sc_count.csv", "csv")
    sc_meta = _load_data(f"{args.data_dir}/sc_meta.csv", "csv")

    ## visualization settings
    colors = [
        "#FFD92F", "#4DAF4A", "#FCCDE5", "#D9D9D9", "#377EB8", "#7FC97F", "#BEAED4",
        "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666", "#1B9E77",
        "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D"
    ]

    visualize_ct = [
        "Acinar_cells", "Cancer_clone_A", "Cancer_clone_B", "Ductal_terminal_ductal_like",
        "Ductal_CRISP3_high-centroacinar_like", "Ductal_MHC_Class_II",
        "Ductal_APOL1_high-hypoxic", "Fibroblasts"
    ]

    gene_visualize = ["Tm4sf1", "S100a4", "Tff3", "Apol1", "Crisp3", "CD248"]

    CARD_obj = None

    #### standard CARD deconvolution
    if args.run_deconvolution:
        CARD_obj = _run_CARD_deconvolution(
            sc_count=sc_count,
            sc_meta=sc_meta,
            spatial_count=spatial_count,
            spatial_location=spatial_location,
            phi_grid=phi_grid,
            warm_start_phi=args.warm_start_phi,
            max_iter=args.max_iter,
            epsilon=args.epsilon,
            verbose=args.verbose,
        )

    #### scMapping
    if args.run_scmapping:
        if CARD_obj is None:
            raise RuntimeError("scMapping requires standard CARD deconvolution to run first.")

        scMapping = _run_CARD_scMapping(
            CARD_obj,
            shapeSpot="Square",
            numCell=args.scmapping_num_cell,
            ncore=args.scmapping_ncore,
        )

        if args.verbose:
            print(scMapping)

        # Uncomment to visualize
        # MapCellCords = scMapping.obs
        # color_map = [
        #     "#8DD3C7", "#CFECBB", "#F4F4B9", "#CFCCCF", "#D1A7B9", "#E9D3DE", "#F4867C", "#C0979F",
        #     "#D5CFD6", "#86B1CD", "#CEB28B", "#EDBC63", "#C59CC5", "#C09CBF", "#C2D567", "#C9DAC3",
        #     "#E1EBA0", "#FFED6F", "#CDD796", "#F8CDDE"
        # ]
        # CARDscMapping_visual_functions = [
        #     ("visual_scmapping", {"MapCellCords": MapCellCords, "colors": color_map})
        # ]
        # _run_CARD_visualization(CARDscMapping_visual_functions)

    #### imputation
    if args.run_imputation:
        if CARD_obj is None:
            raise RuntimeError("Imputation requires standard CARD deconvolution to run first.")

        CARDimput_obj = _run_CARD_imputation(
            CARD_obj,
            num_grids=args.imputation_num_grids,
            in_neighbor=args.imputation_in_neighbor,
            exclude=None,
            showGrid=args.imputation_show_grid,
        )

        if args.verbose:
            refined_prop = CARDimput_obj.uns["info_parameters"]["refined_prop"]
            print(refined_prop.iloc[:2, :2])

        # Uncomment to visualize
        # refined_prop = CARDimput_obj.uns["info_parameters"]["refined_prop"]
        # refined_spatial = CARDimput_obj.uns["info_parameters"]["location_imputation"]
        # refined_expression = CARDimput_obj.uns["info_parameters"]["refined_expression"]
        # origin_expression = CARDimput_obj.uns["spatial_countMat"]
        # origin_spatial = CARDimput_obj.uns["spatial_location"]
        #
        # CARDimput_visual_functions = [
        #     ("visual_prop", {
        #         "proportion": refined_prop,
        #         "spatial_location": refined_spatial,
        #         "ct_visualize": visualize_ct,
        #         "colors": ["lightblue", "lightyellow", "red"],
        #         "num_cols": 4,
        #         "point_size": 20.0
        #     }),
        #     ("visual_gene", {
        #         "spatial_expression": refined_expression,
        #         "spatial_location": refined_spatial,
        #         "gene_visualize": gene_visualize,
        #         "colors": None,
        #         "num_cols": 6,
        #         "point_size": 50.0
        #     }),
        #     ("visual_gene", {
        #         "spatial_expression": origin_expression,
        #         "spatial_location": origin_spatial,
        #         "gene_visualize": gene_visualize,
        #         "colors": None,
        #         "num_cols": None,
        #         "point_size": 18.0
        #     })
        # ]
        # _run_CARD_visualization(CARDimput_visual_functions)

    #### CARDfree
    if args.run_cardfree:
        CARDfree_obj = _run_CARDfree_deconvolution(
            spatial_count=spatial_count,
            spatial_location=spatial_location,
            nmfSelect=args.cardfree_nmf_select,
        )

        if args.verbose:
            free_prop = CARDfree_obj.uns["info_parameters"]["Proportion_CARD"]
            print(free_prop.iloc[0:2])

        # Uncomment to visualize
        # free_prop = CARDfree_obj.uns["info_parameters"]["Proportion_CARD"]
        # free_prop = free_prop.iloc[:, [7, 9, 13, 1, 0, 5, 11, 17, 6, 12, 19, 18, 15, 16, 10, 14, 3, 8, 2, 4]]
        # free_prop.columns = [f"CT{i}" for i in range(1, 21)]
        # free_spatial = CARDfree_obj.uns["spatial_location"]
        # free_visualize_ct = ["CT1", "CT2", "CT3", "CT4", "CT5", "CT6", "CT7", "CT8"]
        #
        # CARDfree_visual_functions = [
        #     ("visual_pie", {
        #         "proportion": free_prop,
        #         "spatial_location": free_spatial,
        #         "colors": colors,
        #         "radius": 1.4,
        #         "sort": False,
        #         "fontsize": 12,
        #         "legend_col": 5
        #     }),
        #     ("visual_prop", {
        #         "proportion": free_prop,
        #         "spatial_location": free_spatial,
        #         "ct_visualize": free_visualize_ct,
        #         "colors": ["lightblue", "lightyellow", "red"],
        #         "num_cols": 4,
        #         "point_size": 25.0
        #     }),
        #     ("visual_prop_2CT", {
        #         "proportion": free_prop,
        #         "spatial_location": free_spatial,
        #         "ct2_visualize": ["CT1", "CT2"],
        #         "colors": [["lightblue", "lightyellow", "red"], ["lightblue", "lightyellow", "black"]]
        #     }),
        #     ("visual_Cor", {"proportion": free_prop, "colors": None}),
        # ]
        # _run_CARD_visualization(CARDfree_visual_functions)


if __name__ == "__main__":
    main()