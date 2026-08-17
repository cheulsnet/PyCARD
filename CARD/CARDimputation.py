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
import pandas as pd
import numpy as np
import geopandas as gpd

from alphashape import alphashape
from shapely.geometry import Polygon, Point, MultiPolygon
from shapely.ops import unary_union
from scipy.sparse import csr_matrix, diags
from sklearn.neighbors import NearestNeighbors

def sample_grid_within(location, num_sample, concavity=1.0):
    """
    Make new spatial locations on unmeasured tissue through grids.

    Parameters
    ----------
    location : DataFrame
        Spatial location data frame of the original spatial resolved transcriptomics dataset, stored in the CARD_object@spatial_location.
    numSample : Numeric
        Approximate number of cells in grid within the shape of the spatial location data frame.
    concavity : Numeric, optional
        A relative measure of concavity. The default is 2.0, which can produce detailed enough shapes. Infinity results in a convex hull while 1 results in a more detailed shape.

    Returns
    -------
    list of DataFrame
        A list of DataFrames with newly gridded points. Each DataFrame contains randomly sampled points within the concave polygon.

    Examples
    --------
    # Replace location with an actual spatial location DataFrame.
    result = sample_grid_within(location, num_sample=100)

    """  
    # Create a GeoDataFrame from the location DataFrame
    gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(location['x'], location['y']))
    
    # Step 1: Construct the concave hull using alphashape
    # higher concavity at R concaveman, lower alpha at python alphashape
    alpha_shape = alphashape(gdf.geometry, concavity)
    
    if alpha_shape.is_empty:
        raise ValueError("Cannot compute concave hull. Change concavity value.")
    
    # Step 2: Extract the exterior boundary of the concave hull
    boundary_polygons = []
    if isinstance(alpha_shape, MultiPolygon):
        for polygon in alpha_shape.geoms:
            boundary_polygons.append(Polygon(polygon.exterior))
    else:
        boundary_polygons.append(Polygon(alpha_shape.exterior))
    
    # Connect the constituent polygons into a single Polygon
    unified_polygon = unary_union(boundary_polygons)
      
    # Step 3: Fill the new num_sample points inside the exterior boundary
    random_points = []
    minx, miny, maxx, maxy = unified_polygon.bounds
    
    # Adjust the grid spacing to ensure the desired number of points
    num_points = int(np.sqrt(num_sample))
    x_points = np.linspace(minx, maxx, num_points)
    y_points = np.linspace(miny, maxy, num_points)
    
    for x in x_points:
        for y in y_points:
            point = Point(x, y)
            if point.within(unified_polygon):
                random_points.append(point)
    
    # Create the output GeoDataFrame
    data = gpd.GeoDataFrame(geometry=random_points)
    
    return data

    
def norm_coords_train_test(location_combine, train_ind, test_ind):
    """
    Normalize the new spatial locations without changing the shape and relative positions.

    Parameters:
    - locationOrig (pd.DataFrame): combined spatial location data frame of the original spatial locations and imputation locations.
    - trainInd (pd.Series or list): Index of the original spatial locations.
    - testInd (pd.Series or list): Index of the newly gridded spatial locations.

    Returns:
    pd.DataFrame: The normalized spatial location data frame containing the normalized spatial coordinates
                  for both training and test data. Normalization is done to a 0-1 scale.

    Note:
    The normalization process typically involves scaling the coordinates to fit within a range of 0 to 1,
    preserving the shape and relative positions of the spatial locations.

    Example usage:
    normalized_data = norm_coords_train_test(location_orig, train_ind, test_ind)

    """
    # Normalize to 0-1 scale
    
    # (1) Normalize training coordinates
    # extract norm_coords_train containing the spatial coordinates (x and y) for the training indices.
    norm_coords_train = location_combine.iloc[train_ind - 1, [1, 2]].copy()
    # calculate the minimum values for both x and y coordinates from the training data.
    location_factor_x = norm_coords_train['x'].min()
    location_factor_y = norm_coords_train['y'].min() 
    # subtract the minimum values from the x and y coordinates in the training data
    norm_coords_train['x'] = norm_coords_train['x'] - location_factor_x
    norm_coords_train['y'] = norm_coords_train['y'] - location_factor_y
    # calculate the scale factor as the maximum value between the maximum x and y coordinates in the training data.
    scale_factor = norm_coords_train[['x', 'y']].max().max()
    # scale the x and y coordinates in the training data by dividing them by the computed scale factor.
    norm_coords_train['x'] = norm_coords_train['x'] / scale_factor
    norm_coords_train['y'] = norm_coords_train['y'] / scale_factor
    
    # (2) Normalize test coordinates
    norm_coords_test = location_combine.iloc[test_ind - 1, [1, 2]].copy()
    # subtract the minimum values from the x and y coordinates in the test data
    norm_coords_test['x'] = norm_coords_test['x'] - location_factor_x
    norm_coords_test['y'] = norm_coords_test['y'] - location_factor_y
    # scale the x and y coordinates in the training data by dividing them by the computed scale factor.
    norm_coords_test['x'] = norm_coords_test['x'] / scale_factor
    norm_coords_test['y'] = norm_coords_test['y'] / scale_factor
    
    return pd.concat([norm_coords_test, norm_coords_train])


def sigma(location_combine, train_ind, test_ind, optimal_phi, in_neighbor):
    """
    Calculate the variance covariance matrix used in the imputation of the newly gridded locations.

    Parameters:
    - location_combine (pd.DataFrame): combined spatial location data frame of the original spatial locations and imputation locations..
    - trainInd (pd.Series or list): Index of the original spatial locations.
    - testInd (pd.Series or list): Index of the newly gridded spatial locations.
    - optimalPhi (float): The optimal phi value stored in CARD_object.
    - inNeighbor (int): Number of neighbors used in the imputation on newly gridded spatial locations. Default is 10.

    Returns:
    dict: A dictionary containing the computed spatial covariance matrices and the weighted adjacency matrix.
          The keys in the dictionary correspond to different components of the covariance matrices.

    Notes:
    - The function relies on the 'norm_coords_train_test' function to normalize the spatial coordinates.

    Example usage:
    cov_matrices = sigma(location_orig, train_ind, test_ind, optimal_phi, in_neighbor)

    """
    # call the norm_coords_train_test function to normalize the spatial coordinates.
    norm_coords_temp = norm_coords_train_test(location_combine, train_ind, test_ind)
    
    # (1) find neighbors
    # use NearestNeighbors from scikit-learn to find the nearest neighbors for each point in the normalized coordinates.
    near_data = NearestNeighbors(n_neighbors=in_neighbor + 1).fit(norm_coords_temp.iloc[:, :2])
    distances, neighbors = near_data.kneighbors(norm_coords_temp.iloc[:, :2])
    
    # Remove the location itself as the neighbor (as the first column corresponds to the location itself.)
    neighbors = neighbors[:, 1:]
    
    # (2) Construct WTemp
    # Flatten the neighbors array and create row indices
    row_indices = np.repeat(np.arange(len(neighbors)), in_neighbor)
    # Flatten the neighbors array and use them as column indices
    column_indices = neighbors.flatten()
    n_rows = len(neighbors)
    n_cols = len(neighbors)
    
    isigma = 0.1
    kernel_neibors = np.exp(-distances[:, 1:] ** 2 / (2 * isigma ** 2))
    wtemp = csr_matrix((kernel_neibors.flatten(), (row_indices, column_indices)), shape=(n_rows, n_cols))
    
    # (3) Construct SigmaTemp
    dtemp = diags(wtemp.sum(axis=1).A.ravel())
    w22 = wtemp[len(test_ind):, len(test_ind):]
    sigma_temp = dtemp - optimal_phi * wtemp

    # (4) Build the output sigmas
    sigma11 = sigma_temp[:len(test_ind), :len(test_ind)]
    sigma12 = sigma_temp[:len(test_ind), len(test_ind):]
    sigma21 = sigma_temp[len(test_ind):, :len(test_ind)]
    sigma22 = sigma_temp[len(test_ind):, len(test_ind):]

    # return a dictionary containing the computed spatial covariance matrices and the weighted adjacency matrix. 
    # The keys in the dictionary correspond to different components of the covariance matrices.
    return {
        'SigmaTemp': sigma_temp,
        'Sigma11': sigma11,
        'Sigma12': sigma12,
        'Sigma21': sigma21,
        'Sigma22': sigma22,
        'W22': w22,
    }


def mvn_cv(Vtrain, location_combine, train_ind, test_ind, B, optimal_b, optimal_phi, lambd, in_neighbor):
    """
    Imputation and Construction of High-Resolution Spatial Maps for Cell Type Composition and Gene Expression
    by the spatial correlation structure between original spatial locations and new gridded spatial locations.

    Parameters:
    - Vtrain (np.DataFrame): Matrix, estimated V matrix from CARD.
    - location_combine (pd.GeoDataFrame): combined spatial location data frame of the original spatial locations and imputation locations.
    - train_ind (np.ndarray): Index of the original spatial locations.
    - test_ind (np.ndarray): Index of the newly gridded spatial locations.
    - B (np.DataFrame): Matrix, used in the deconvolution as the reference basis matrix.
    - XinputNorm (np.ndarray): Matrix, used in the deconvolution as the normalized spatial count data.
    - optimal_b (np.ndarray): Vector, vector of the intercept for each cell type estimated based on the original spatial resolution.
    - optimal_phi (float): The optimal phi value stored in CARD_object.
    - lambd (list): Vector, vector of cell type-specific scalar in the CAR model.
    - in_neighbor (int): Number of neighbors used in the imputation on newly gridded spatial locations. Default is 10.

    Returns:
    dict: A dictionary with the imputed Cell type composition Vtest matrix on the newly gridded spatial locations
          and predicted normalized gene expression.

    Notes:
    - The function relies on the 'sigma' function to calculate spatial covariance matrices.

    Example usage:
    result_dict = mvn_cv(Vtrain, location_orig, train_ind, test_ind, B, xinput_norm, optimal_b, optimal_phi, lambd, in_neighbor)

    """
    #  calculate spatial covariance matrices (sigma11, sigma21, sigma12, sigma22, w22) using the sigma function
    SigmaList = sigma(location_combine, train_ind, test_ind, optimal_phi, in_neighbor)
    Sigma11 = SigmaList["Sigma11"]
    Sigma21 = SigmaList["Sigma21"]
    Sigma12 = SigmaList["Sigma12"]
    Sigma22 = SigmaList["Sigma22"]
    W22 = SigmaList["W22"]
    
    Mu_cond = np.zeros((len(test_ind), len(lambd)))
    
    for ict in range(len(lambd)):
        temp = Sigma12.dot((Vtrain.iloc[:, ict].values - optimal_b[:Vtrain.shape[0], ict]).reshape(-1, 1))
        Sigma11_dense = Sigma11.toarray()
        solver = np.linalg.solve(Sigma11_dense, temp)
        Mu = optimal_b[:len(test_ind), ict] - solver.flatten()
        Mu_cond[:, ict] = Mu
    
    # calculate Vtest
    Vtest = Mu_cond
    Vtest = np.divide(Vtest, np.sum(Vtest, axis=1, keepdims=True))
    colnames_Vtrain = list(Vtrain.columns)
    Vtest = pd.DataFrame(Vtest, index=test_ind, columns=colnames_Vtrain)
    # reset the index with default starting from 0
    Vtest.reset_index(drop=True, inplace=True)
    
    # calculate XtestHat
    XtestHat = B.dot(Vtest.T)
    
    return {"Vtest": Vtest, "XtestHat": XtestHat}


def CARD_imputation(CARD_object, num_grids, in_neighbor=10, exclude=None):
    """
    Construct an enhanced spatial expression map on the unmeasured tissue locations.

    Parameters:
    - card_object (CARD Object): CARD Object with estimated cell type compositions on the original spatial resolved transcriptomics data.
    - num_grids (int): Initial number of newly gridded spatial locations. 
                      The final number of newly gridded spatial locations will be lower than this value 
                      since the newly gridded locations outside the shape of the tissue will be filtered.
    - in_neighbor (int): Number of neighbors used in the imputation on newly gridded spatial locations. Default is 10.
    - exclude (list or None): Vector, the row names of spatial location data on the original resolution that you want to exclude. 
                             This is to avoid the weird detection of the shape.

    Returns:
    CARD Object: The updated CARD object with the refined cell type compositions estimated for newly gridded spots 
                 and the refined predicted gene expression (normalized).

    Notes:
    - The function relies on the 'mvn_cv' function for imputation.

    Example usage:
    result_card = card_imputation(CARD_object, num_grids, in_neighbor, exclude)

    """
    # Extract results from CARD object
    # The B, xinput_norm, and Vtrain are matrices, 
    # and location is a DataFrame containing spatial information.
    B = CARD_object.uns['info_parameters']['algorithm_matrix']['B']
    Vtrain = CARD_object.uns['info_parameters']['algorithm_matrix']['Res']['V']
    location = CARD_object.uns['spatial_location']
    
    # Check  if the row indices of the spatial locations and Vtrain match. 
    # If they do, it prints a message indicating that the row names are matched.
    if sum(location.index == Vtrain.index) == len(Vtrain):
        print("## The row names of locations are matched ...")

    # Make new grids
    # This part prepares the spatial locations for sampling. 
    # If exclude is provided, it removes the specified indices from the location DataFrame.
    print("## Make grids on new spatial locations ...")
    if exclude is not None:
        location_use_to_sample = location[~location.index.isin(exclude)]
    else:
        location_use_to_sample = location

    # Uses the sample_grid_within function to generate new grids on the spatial locations.
    # higher concavity at R concaveman, lower alpha at python alphashape
    # recommended concavity(alpha) which less than or equal to 1.0 
    data = sample_grid_within(location_use_to_sample, num_grids, concavity=1.0)
    
    # Extract x, y coordinates from the 'geometry' column in 'data'
    data['x'] = data['geometry'].x
    data['y'] = data['geometry'].y
    
    # Align 'location' format to match 'data'
    location['geometry'] = gpd.points_from_xy(location['x'], location['y'])

    # Concatenate 'location' and 'data' based on the x, y coordinates
    location_combine = pd.concat([data, location]).reset_index(drop=True)
    
    # Drop duplicates based on x, y coordinates
    location_combine = location_combine.drop_duplicates(subset=['x', 'y']).reset_index(drop=True)
    
    # Get indices from 1 to the number of rows in 'data'
    test_ind = np.arange(1, len(data) + 1)
    # Get indices from the number of rows in 'data' + 1 to the number of rows in 'location_combine'
    train_ind = np.arange(len(data) + 1, len(location_combine) + 1)
    
    # Results from the CARD
    # Extract optimal_b from CARD_object
    optimal_b = np.array(CARD_object.uns['info_parameters']['algorithm_matrix']['Res']['b'])
    
    # Create a vector of ones with the same length as the number of rows in 'location_combine'
    vecOne = np.ones((location_combine.shape[0], 1))
    # Reshape optimal_b to ensure proper matrix multiplication
    optimal_b = np.matmul(vecOne, optimal_b.T)
    # Extract OptimalPhi and lambda from CARD_object
    optimal_phi = CARD_object.uns['info_parameters']['phi']
    lambd = CARD_object.uns['info_parameters']['algorithm_matrix']['Res']['lambda']
    
    # calls the mvn_cv function to perform imputation based on the provided parameters and data.
    imputation = mvn_cv(Vtrain, location_combine, train_ind, test_ind, B, optimal_b, optimal_phi, lambd, in_neighbor)
    
    # Build index of imputation['Vtest'] and columns of imputation['XtestHat']
    Vtest_col = data.iloc[:, [1, 2]].copy()
    
    # using 10.0x10.0 format
    refined_prop_columns = Vtest_col['x'].astype(str) + 'x' + Vtest_col['y'].astype(str)
    # using 10x10 format (integer format)
    #refined_prop_columns = Vtest_col['x'].astype(int).astype(str) + 'x' + Vtest_col['y'].astype(int).astype(str)
    
    imputation['Vtest'].index = refined_prop_columns
    imputation['XtestHat'].columns = refined_prop_columns
    
    # updates the card_object with the imputed values for refined proportions (refined_prop) 
    # and refined expression values (refined_expression).
    CARD_object.uns['info_parameters']['refined_prop'] = imputation['Vtest']
    CARD_object.uns['info_parameters']['refined_expression'] = imputation['XtestHat']
    
    # updates the refined spatial locations including the new spatial points
    CARD_object.uns['info_parameters']['location_imputation'] = Vtest_col
    
    return CARD_object



