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
import pyreadr
from pyreadr.custom_errors import LibrdataError

def _load_data(file_path, file_format):
    """
    Load data from file into a DataFrame.

    Parameters:
        file_path (str): Path to the data file.
        file_format (str): Format of the data file (e.g., 'csv', 'feather', 'xlsx', 'RData').

    Returns:
        DataFrame: Loaded data as a pandas DataFrame.
    """
    if file_format.lower() == 'csv':
        return pd.read_csv(file_path, index_col=0) # Assuming the first column should be the index
    elif file_format.lower() == 'feather':
        df = pd.read_feather(file_path)
        df.set_index(df.columns[0], inplace=True)  # Assuming the first column should be the index
        return df
    elif file_format.lower() == 'xlsx':
        return pd.read_excel(file_path, index_col=0) # Assuming the first column should be the index
    elif file_format.lower() == 'rdata':
        try:
            result = pyreadr.read_r(file_path)
            if result:
                return next(iter(result.values()))
            else:
                raise ValueError("The file '{}' does not contain any data.".format(file_path))
        except LibrdataError:
            raise ValueError("The file '{}' does not contain recognizable data.".format(file_path))
          
    else:
        raise ValueError(f"Unsupported file format: {file_format}")

