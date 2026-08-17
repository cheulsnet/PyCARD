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
from scipy.optimize import nnls
from scipy.spatial.distance import cdist
from anndata import AnnData
from concurrent.futures import ThreadPoolExecutor

def getWeightForCell(sc_eset, sc_gene, B):
    """
    Estimate the cell type composition signature for each single cell in the scRNaseq reference data.

    Parameters
    ----------
    sc_eset : AnnData
        The scRNaseq reference data stored in the CARD object.
    sc_gene : AnnData
        the index of sc_count (gene names).
    B : pd.DataFrame
        Reference basis matrix (stored in the CARD object).

    Returns
    -------
    pd.DataFrame
        Matrix of the cell type composition signature for each single cell in the scRNaseq reference.

    Examples
    --------
    getWeightForCell(sc_eset, sc_gene, ct_varname, ct_select, sample_varname, B)
    """
    # Extract count matrix from AnnData object and transpose it (genes x cellIDs)
    count = pd.DataFrame(sc_eset.X.T.toarray(), index=sc_gene)
    
    # Filter columns with non-zero sums
    col_sums = count.sum(axis=0)  # Calculate column sums
    count = count.loc[:, col_sums > 0]
    
    # Normalize each column by dividing by its sum
    count = count.div(col_sums, axis=1)
    
    # Keep only rows present in B
    count = count[count.index.isin(B.index)]
    # Reorder rows to match B
    count = count.reindex(B.index)
    
    # Initialize an empty array to store the results
    Mean_Cell = np.empty((B.shape[1], count.shape[1]))
    
    # Perform nnls (non-negative least squares regression) for each column of count
    for i in range(count.shape[1]):
        mod1 = nnls(B, count.iloc[:, i])[0]
        Mean_Cell[:, i] = mod1
 
    # Transpose Mean_Cell
    Mean_Cell = Mean_Cell.T

    # Set row names to match column names of count
    Mean_Cell = pd.DataFrame(Mean_Cell, index=np.arange(count.shape[1]))

    # Set column names to match column names of B
    Mean_Cell.columns = B.columns
    
    return Mean_Cell 
    
# This helper function generates points within a square centered at (Cords[i, 0], Cords[i, 1]) with side length 2 * min_distance. 
# It returns a DataFrame with columns 'x' and 'y'.
def getPointsWithinSquare(Cords, i, numCell, min_distance):
    minX = min_distance
    minY = min_distance
    rectangular_x = np.random.uniform(low=Cords.iloc[i, 0] - minX, high=Cords.iloc[i, 0] + minX, size=numCell)
    rectangular_y = np.random.uniform(low=Cords.iloc[i, 1] - minY, high=Cords.iloc[i, 1] + minY, size=numCell)
    df = pd.DataFrame({'x': rectangular_x, 'y': rectangular_y})
    return df

# This helper function generates random points within a circle. 
def runifdisc(num_points, radius, centre=(0, 0)):
    # Generate random angles
    angles = np.random.uniform(0, 2 * np.pi, num_points)
    # Generate random radii uniformly within the given radius
    radii = np.sqrt(np.random.uniform(0, radius**2, num_points))
    # Convert polar coordinates to Cartesian coordinates
    x = radii * np.cos(angles) + centre[0]
    y = radii * np.sin(angles) + centre[1]
    return np.column_stack((x, y))

# This helper function generates points within a circle centered at (Cords[i, 0], Cords[i, 1]) with radius min_distance. 
# It returns a DataFrame with columns 'x' and 'y'.
def getPointsWithinCircle(Cords, i, numCell, min_distance):
    circle = runifdisc(numCell, radius=min_distance, centre=(Cords.iloc[i, 0], Cords.iloc[i, 1]))
    df = pd.DataFrame({'x': circle[:, 0], 'y': circle[:, 1]})
    return df

def getHighresCords(Cords, numCell, shape='Square'):
    """
    Sample the spatial location information for each single cell.

    Parameters
    ----------
    Cords : pd.DataFrame
        Spatial location information in the measured spatial locations, with the first and second columns representing the 2-D x-y coordinate system.
    numCell : int
        Numeric value indicating the number of single cells in each measured location.
    shape : str, optional
        Character indicating whether the sampled spatial coordinates for single cells are located in a Square-like region or a Circle-like region. 
        The center of this region is the measured spatial location in the non-single cell resolution spatial transcriptomics data.
        Default is 'Square'; the other shape is 'Circle'.

    Returns
    -------
    pd.DataFrame
        DataFrame with the sampled spatial location information for each single cell.

    Examples
    --------
    getHighresCords(spatial_location_data, 20, 'Square')
    """
    # Use scipy.spatial.distance.cdist to compute pairwise Euclidean distances (ED) between the coordinates in Cords.
    # Calculate the minimum distance to the nearest cell for each cell and store it in the array dis.
    ED = cdist(Cords, Cords, metric='euclidean')

    n = Cords.shape[0]
    dis = np.array([min(ED[i, np.arange(n) != i]) for i in range(n)])

    # Set min_distance to half of the median of the minimum distances.
    min_distance = np.median(dis) / 2

    # Initialize an empty Pandas DataFrame Cords_new to store the high-resolution coordinates
    Cords_new = pd.DataFrame()

    # Loop through each cell and generate high-resolution points based on the specified shape.
    for i in range(n):
        if shape == 'Square':
            df = getPointsWithinSquare(Cords, i, numCell, min_distance)
        elif shape == 'Circle':
            df = getPointsWithinCircle(Cords, i, numCell, min_distance)

        # Add additional columns to the DataFrame, including the centerSPOT identifier and center coordinates.
        df['centerSPOT'] = Cords.iloc[i, 0].astype(str) + "x" + Cords.iloc[i, 1].astype(str)
        df['centerx'] = Cords.iloc[i, 0]
        df['centery'] = Cords.iloc[i, 1]

        # Concatenate the generated DataFrame df to the main DataFrame Cords_new.
        Cords_new = pd.concat([Cords_new, df], ignore_index=True)

    # Remove duplicate rows based on the 'x' and 'y' columns.
    Cords_new = Cords_new[~Cords_new.duplicated(subset=['x', 'y'])]
    # Set the index of Cords_new to a combination of 'x' and 'y'.
    Cords_new['row_names'] = Cords_new['x'].astype(str) + "x" + Cords_new['y'].astype(str)
    
    return Cords_new.set_index('row_names')


def AssignSCcords(MappingSpotCellCor, Cords_new, numCell, sc_eset, ct_varname):
    """
    This function assigns spatial location information for each single cell.

    Parameters
    ----------
    MappingSpotCellCor : pd.DataFrame
        A mapped correlation matrix indicating the relationship between each measured spatial location and the single cell in the scRNAseq reference.
    Cords_new : pd.DataFrame
        Output from the function getHighresCords.
    numCell : int
        A numeric value indicating the number of single cells in each measured location.
    sc_eset : SingleCellExperiment
        A single-cell experiment object stored in the CARD object.
    ct_varname : str
        The name of the column in metaData that specifies the cell type annotation information, stored in the CARD object.

    Returns
    -------
    pd.DataFrame
        DataFrame containing assigned spatial location information.

    Examples
    --------
    AssignSCcords(MappingSpotCellCor, Cords_new, 20, sc_eset, 'cell_type')
    """
    # initializes an empty Pandas DataFrame MapCellCords to store the results.
    MapCellCords = pd.DataFrame()
    
    # Loop over each spot (ispot) in MappingSpotCellCor.
    for ispot in range(1, MappingSpotCellCor.shape[0] + 1):
        mapCell = MappingSpotCellCor.iloc[ispot - 1, :]
        
        # Sort the values in the row corresponding to the spot in descending order.
        mapCell = mapCell.sort_values(ascending=False)
        
        # Extract the coordinates of the spot (centerspot) from the index.
        centerspot = MappingSpotCellCor.index[ispot - 1]
        ispot_cords = pd.DataFrame({
            'x': [float(centerspot.split('x')[0])],
            'y': [float(centerspot.split('x')[1])]
        })

        # Extract the subset of Cords_new where centerSPOT matches the current spot.
        subCordsCell = Cords_new[Cords_new['centerSPOT'] == centerspot]
        
        # Calculate Euclidean distances (EDwithCenter) between each cell in the subset and the spot's coordinates.
        EDwithCenter = cdist(subCordsCell[['x', 'y']], ispot_cords)
        
        # Add the EDwithCenter column to subCordsCell
        subCordsCell = subCordsCell.copy()
        subCordsCell.loc[:, 'EDwithCenter'] = EDwithCenter.flatten()
        
        # sort subCordsCell by the calculated distances (EDwithCenter).
        subCordsCell = subCordsCell.sort_values(by='EDwithCenter')
        
        # Create a temporary DataFrame (mapCellCordsTemp) with the top cells associated with the spot.
        mapCellCordsTemp = pd.DataFrame({'CorwithSpot': mapCell[:numCell]})
        
        # Add columns for correlation with the spot (CorwithSpot)
        # Concatenating along the column axis and ignoring the indexes (not alligned)
        mapCellCordsTemp = pd.concat([mapCellCordsTemp.reset_index(drop=True), subCordsCell.reset_index(drop=True)], axis=1)
        
        # Add columns for correlation with cell type (CT)
        mapCellCordsTemp['CT'] = sc_eset.obs['cellType'].iloc[mapCell[:numCell].index].values
        
        # Add columns for sequence of cell
        mapCellCordsTemp['Cellseq'] = mapCell[:numCell].index
        # Add columns for correlation with cell identifier (cellID)
        mapCellCordsTemp['Cell'] = sc_eset.obs['cellID'].iloc[mapCell[:numCell].index].values
        
        # Concatenate the temporary DataFrame to the main result DataFrame (MapCellCords).
        MapCellCords = pd.concat([MapCellCords, mapCellCordsTemp])

    return MapCellCords


def CARD_SCMapping(CARD_object, numCell, shapeSpot="Square", ncore=10):
    """
    Extension of CARD for performing single-cell mapping from non-single cell spatial transcriptomics dataset.

    Parameters
    ----------
    CARD_object : CARD
        CARD object created by the createCARDObject function.
    numCell : int
        Numeric value indicating the number of single cells in each measured location.
    shapeSpot : str, optional
        Character indicating whether the sampled spatial coordinates for single cells are located in a Square-like region or a Circle-like region.
        The center of this region is the measured spatial location in the non-single cell resolution spatial transcriptomics data.
        Default is "Square"; the other shape is "Circle".
    ncore : int, optional
        Numeric value indicating the number of cores used to accelerate the procedure. Default is 10.

    Returns
    -------
    AnnData
        AnnData object SCE representing the mapped scRNAseq data (the mapped expression at single-cell resolution and the spatial location information of each single cell).

    Examples
    --------
    CARD_SCMapping(my_CARD_object, 20, 'Square', 10)
    """
    # Load sc count data stored in CARD_object
    sc_eset = CARD_object.uns['sc_eset']
    sc_gene = CARD_object.uns['sc_gene']
    
    # Load genes-to-cell type matrix, B created at CARD deconvolution
    B = CARD_object.uns['info_parameters']['algorithm_matrix']['B']
    # sort by column name (unique cell types) in a case-insensitive way
    B = B.sort_index(axis=1, key=lambda x: x.str.casefold())
    
    # Load cell type variance name
    ct_varname = CARD_object.uns['info_parameters']['ct.varname']
    
    # Load proportion_CARD matrix created at CARD deconvolution
    res_CARD = CARD_object.uns['info_parameters']['Proportion_CARD'].copy()
    # Sort by the same column order as B
    res_CARD = res_CARD[B.columns]
    
    # (Step 1) Get weighted cell means sorted by the same column order as B
    Mean_Cell = getWeightForCell(sc_eset, sc_gene, B)
    
    # (Step 2) Get spatial coordinates
    # Extracting spatial_location from CARD_object
    spatial_location = CARD_object.uns['spatial_location']
    Cords = spatial_location.copy()
    Cords_new = getHighresCords(Cords, numCell=numCell, shape=shapeSpot)
    
    # Transpose Mean_Cell dataframe
    Mean_Cell_transposed = Mean_Cell.T

    # Calculate correlations between each row of res_CARD and each row of Mean_Cell_transposed
    MappingSpotCellCor = res_CARD.apply(lambda row1: Mean_Cell_transposed.apply(lambda row2: row1.corr(row2)), axis=1)
     
    # drop the 'NaN'
    MappingSpotCellCor.dropna(axis=1, inplace=True)
           
    # Create row names based on 'x' and 'y' columns of Cords DataFrame
    MappingSpotCellCor.index = Cords['x'].astype(str) + "x" + Cords['y'].astype(str)
    
    # Assign high-resolution coordinates
    MapCellCords = AssignSCcords(MappingSpotCellCor, Cords_new, numCell, sc_eset, ct_varname)
    
    ### The scRNAseq count data is mapped to high-resolution coordinates.
    counts = sc_eset.X.T
    # Convert sparse matrix to dense array
    counts = counts.toarray()
    # Create DataFrame  
    count_sc = pd.DataFrame(counts)
    
    # Filter out columns with zero sums
    count_sc = count_sc.loc[:, count_sc.sum(axis=0) > 0]
    # count_sc = sc_eset.X[:, sc_eset.X.sum(axis=0).A1 > 0] # simpler version without dataframe conversion
    
    # intialize the count_CT
    count_CT = []
    
    # Define the function to process each row of res_CARD in parallel
    def process_row(ispot, res_CARD, MapCellCords, count_sc):
        spot = res_CARD.index[ispot]
        MapCellCords_spot = MapCellCords[MapCellCords['centerSPOT'] == spot]
        # Pick the cells specified in 'Cellseq'
        df = count_sc.iloc[:, MapCellCords_spot['Cellseq']]
        colnames = [f"{MapCellCords_spot['Cell'].iloc[i]}:{spot}:{MapCellCords_spot['x'].iloc[i]}x{MapCellCords_spot['y'].iloc[i]}" for i in range(len(MapCellCords_spot))]
        df.columns = colnames
        return df

    # Perform parallel processing for each row of res_CARD
    with ThreadPoolExecutor(max_workers=ncore) as executor:
        count_CT_parallel = list(executor.map(lambda i: process_row(i, res_CARD, MapCellCords, count_sc), range(len(res_CARD))))

    # The results are concatenated.
    count_CT = pd.concat(count_CT_parallel , axis=1)
    
    # Create AnnData object representing the mapped scRNAseq data
    MapCellCords.index = count_CT.columns
    sce = AnnData(
        X=count_CT.values.T,  # counts (transpose to match dimensions)
        obs=MapCellCords
    )
    
    return sce


