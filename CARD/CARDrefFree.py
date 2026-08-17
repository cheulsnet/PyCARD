
####################################################################################################
## Package : CARD
## Version : 1.0.1
## Date    : 2021-1-7 09:10:08
## Modified: 2021-12-13 16:18:07
## Title   : Spatially Informed Cell Type Deconvolution for Spatial Transcriptomics by CARD.
## Authors : Ying Ma
## Contacts: yingma@umich.edu
##           University of Michigan, Department of Biostatistics
####################################################################################################

import numpy as np
import pandas as pd
from nimfa import Nmf
from sklearn.decomposition import NMF
import mlpack as ml
from scipy.spatial.distance import cdist
from scipy.sparse import csr_matrix

import os
import sys
build_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../build'))
sys.path.insert(0, build_dir)

from CARDfree_module import CARDfree
from CARD.pyNMFpackage import pyNMFselection

def CARD_refFree(CARDfree_object):
    """
    Extension of CARD into a reference-free version of deconvolution: CARDfree.

    Parameters
    ----------
    CARDfree_object : CARDfree
        CARDfree object created by the createCARDfreeObject function.

    Returns
    -------
    CARD
        A CARD object with estimated cell type proportion stored in object@Proportion_CARD. 
        Because this is a reference-free version, the columns of estimated proportion are not cell types but cell type clusters.

    Examples
    --------
    card_ref_free(my_CARDfree_object)

    Notes
    -----
    Additional notes or explanations, if needed.

    See Also
    --------
    createCARDfreeObject : Function to create a CARDfree object.

    References
    ----------
    Any relevant references or citations.

    """    
    # Retrieve relevant data from the CARDfree_object, such as spatial_countMat and spatial_location.
    spatial_countMat = pd.DataFrame(CARDfree_object.X).T 
    spatial_gene = CARDfree_object.uns['spatial_gene']
    spatial_countMat.index = spatial_gene
    
    # Get the selected NMF package id
    nmfSelect = CARDfree_object.uns['nmfSelect']
    
    # Calculate the number of cell types (numK) and create a unique list of markers.
    markerList = CARDfree_object.uns['markerList']
    numK = len(markerList)
    
    # First, flatten the list of lists into a single list and convert it into upper-case, 
    # then obtain the list of unique genes from the flattened list
    marker = list(set(val.upper() for sublist in markerList for val in sublist))
        
    print(f"## Number of unique marker genes: {len(marker)} for {numK} cell types ...")

    # Convert spatial_countMat index into upper case
    spatial_countMat.index = spatial_countMat.index.str.upper()
    
    # Get the intersect upper-cased gene names
    commonGene = np.intersect1d(spatial_countMat.index, marker)
   
    # Remove mitochondrial and ribosomal genes
    commonGene = [gene for gene in commonGene if not gene.startswith('MT-')]
    
    # Check if the number of common genes is sufficient; raise an error if not.
    if len(commonGene) < numK * 10:
        raise ValueError("STOP! The average number of unique marker genes for each cell type is less than 20.")

    # Sort rows and columns
    Xinput = spatial_countMat.sort_index()
    Xinput = Xinput.sort_index(axis=1)
    
    # Filter rows by commonGene (uppercase)
    Xinput = Xinput[Xinput.index.str.upper().isin(commonGene)]

    # Filter rows by row sums > 0
    Xinput = Xinput[(Xinput.sum(axis=1) > 0)]
    # Filter columns by column sums > 0
    Xinput = Xinput.loc[:, (Xinput.sum(axis=0) > 0)]

    # Normalize by column sums
    Xinput_norm = Xinput.div(Xinput.sum(axis=0), axis=1)
    
    ###### Initialization ###############
    ##### Use NMF to factorize the normalized count matrix into matrices B and Vint1.
    ##### Select NMF Package #######
    # 1: NMF from sklearn.decomposition import NMF
    # 2: NMF from nimfa import Nmf (not recommended due to lack of reproductivity)
    # 3: NMF import mlpack as ml
    # 4: NMF customized c++ module using mlpack
    # 5: NMF from the original theory
    # 6: NMF from R
    # 7: NMF from RcppML
    #################################

    # Build NMF package class
    pyNMFselect = pyNMFselection(Xinput_norm, numK=20, seed=20200107, max_iter=1000, iter=50, is_print=False)   
    if Xinput_norm.shape[1] < 5000: 
        # Select NMF package to use
        Basis, Vint1 = pyNMFselect.pyNMF(nmfSelect)    
    else:
        # Use custom mlpack NMF (section id: 4)
        Basis, Vint1 = pyNMFselect.pyNMF(4)

    Basis = pd.DataFrame(Basis, index=commonGene)
    Vint1 = pd.DataFrame(Vint1, index=Xinput_norm.columns)
        
    ##### Spatial Location Normalization
    spatial_location = CARDfree_object.obs
    
    # convert the spatial_location indices into an integer type
    spatial_location.index = spatial_location.index.astype(int)
    spatial_location = spatial_location[spatial_location.index.isin(Xinput_norm.columns)]
    spatial_location = spatial_location.reindex(Xinput_norm.columns)
    
    ##### Normalize the spatial coordinates without changing the shape and relative position
    # Extract 'x' and 'y' columns from spatial_location
    norm_cords = spatial_location[['x', 'y']]
    # Subtract the minimum values from 'x' and 'y'
    norm_cords['x'] = norm_cords['x'] - norm_cords['x'].min()
    norm_cords['y'] = norm_cords['y'] - norm_cords['y'].min()
    # Calculate scaleFactor
    scaleFactor = max(norm_cords['x'].max(), norm_cords['y'].max())
    # Normalize 'x' and 'y' by scaleFactor
    norm_cords['x'] = norm_cords['x'] / scaleFactor
    norm_cords['y'] = norm_cords['y'] / scaleFactor
    
    # Compute Euclidean distances between spatial locations.
    # Euclidean distance matrix
    ED = cdist(norm_cords.values, norm_cords.values, metric='euclidean')

    ###### Set parameters for Gaussian kernel calculation.
    # Construct Gaussian kernel with the default scale/length parameter to be 0.1
    isigma = 0.1
    # Convergence epsilon
    epsilon = 1e-04
    # Gridded values for phi
    phi = [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99]
    # Kernel matrix
    kernel_mat = np.exp(-ED**2 / (2 * isigma**2))
    np.fill_diagonal(kernel_mat, 0) # Set diagonal elements to zero
    
    ###### scale the Xinput_norm and B to speed up the convergence.
    # Calculate the mean of Xinput_norm and B along the rows
    mean_X = np.mean(Xinput_norm)
    mean_B = np.mean(Basis)
    
    # Scale Xinput_norm and B
    Xinput_norm = Xinput_norm * 1e-01 / mean_X
    Basis = Basis * 1e-01 / mean_B
    
    # Build colunms of Xinput_norm and Basis
    Basis.columns = CARDfree_object.uns['marker_ct']
    # using 10.0x10.0 format
    Xinput_norm.columns = spatial_location['x'].astype(str) + 'x' + spatial_location['y'].astype(str)
    # using 10x10 format (integer format)
    #Xinput_norm.columns = spatial_location['x'].astype(int).astype(str) + 'x' + spatial_location['y'].astype(int).astype(str)
    
    ####### Deconvolution Iteration iterating over different values of phi and using card_free function
    ResList = []
    Obj = []

    for phi_value in phi:
        res = CARDfree(
            Xinput_norm.values,
            Basis.values,
            kernel_mat,
            phi_value,
            1000,
            epsilon,
            Vint1.values,
            np.zeros(Basis.shape[1]),
            0.1,
            np.repeat(10, Basis.shape[1])
        )
        
        # Set row and column names for the result's proportion matrix
        res['V'] = pd.DataFrame(res['V'])  # origin
        
        # Store the result in the ResList
        ResList.append(res)

        # Store the objective function value in the Obj list
        Obj.append(res['Obj'])

    ##### Optimal Phi Selection
    # Step 1: Find indices where Obj is equal to its maximum
    optimal_indices = np.where(np.array(Obj) == max(Obj))[0]
    # Step 2: If there are multiple indices, choose the last one
    optimal_index = optimal_indices[-1]
    # Step 3: Get the optimal phi value
    optimal_phi = phi[optimal_index]
    # Step 4: Get the corresponding result
    optimal_res = ResList[optimal_index]
    # Step 5: Set index and columns to optimal_res['V']
    optimal_res['V'].index = Xinput_norm.columns
    optimal_res['V'].columns=Basis.columns
    
    
    print("## Deconvolution Finish! ...")
    
    ######---> 7. Output and Cleanup
    # Updates the CARD_object with the optimal phi
    CARDfree_object.uns['info_parameters']['phi'] = optimal_phi
    
    # Updating proportion matrix: The proportions of cell types (Proportion_CARD) are calculated 
    # by normalizing the results (V) of the optimal deconvolution (OptimalRes) by row sums
    CARDfree_object.uns['info_parameters']['Proportion_CARD'] = optimal_res['V'].div(optimal_res['V'].sum(axis=1), axis=0)
    
    # Updating Algorithm Matrix
    # B: The reference matrix (basis) scaled back to the original scale
    # Xinput_norm: The normalized spatial count matrix (x_input_norm) scaled back to the original scale.
    # Res: The results for the optimal phi (optimal_res).
    optimal_B = np.array(optimal_res['B']) * mean_B / 1e-01
    optimal_B = pd.DataFrame(optimal_B, index=Basis.index, columns=Basis.columns)
    CARDfree_object.uns['info_parameters']['algorithm_matrix'] = {
        'B': optimal_B,
        'Xinput_norm': Xinput_norm * mean_X / 1e-01,
        'Res': optimal_res
    }
    
    # The spatial location is assigned to the spatial_location property of CARD_object
    CARDfree_object.uns['spatial_location'] = spatial_location
    
    # Set the estimated reference matrix attribute of the CARDfree object
    CARDfree_object.uns['info_parameters']['estimated_refMatrix'] = optimal_B.copy()
    
    return CARDfree_object
