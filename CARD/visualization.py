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
import random
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.cm as cm
from matplotlib.colors import Normalize

def create_custom_colormap(colors):
    cmap_segments = [(i / (len(colors) - 1), colors[i]) for i in range(len(colors))]
    return LinearSegmentedColormap.from_list("custom_cmap", cmap_segments)

def draw_pie(dist, xpos, ypos, size, colors, ax=None):
    """
    Drawing the pie chart

    Parameters
    ----------
    dist: 
        List of values representing the distribution of different pie slices.
    xpos: 
        x-coordinate of the center of the pie chart.
    ypos: 
        y-coordinate of the center of the pie chart.
    size: 
        Size of each pie slice.
    colors: 
        List of colors for each pie slice.
    
    Returns
    -------
    pie chart

    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10,8))

    # for incremental pie slices
    cumsum = np.cumsum(dist)
    cumsum = cumsum/ cumsum[-1]
    pie = [0] + cumsum.tolist()

    for i in range(len(pie) - 1):
        r1, r2 = pie[i], pie[i + 1]
        angles = np.linspace(2 * np.pi * r1, 2 * np.pi * r2)
        x = [0] + np.sin(angles).tolist()
        y = [0] + np.cos(angles).tolist()
        
        xy = np.column_stack([x, y])
        ax.scatter([xpos], [ypos], marker=xy, s=size, color=colors[i])

    return ax

def CARD_visualize_prop(proportion, spatial_location, ct_visualize, colors=["lightblue", "lightyellow", "red"], point_size=25.0, num_cols=None):
    """
    Visualize the spatial distribution of cell type proportion.

    Parameters
    ----------
    proportion : pd.DataFrame
        Data frame, cell type proportion estimated by CARD in either original resolution or enhanced resolution.
    spatial_location : pd.DataFrame
        Data frame, spatial location information.
    ct_visualize : list
        Vector of selected cell type names that are interested to visualize.
    colors : list, optional
        Vector of color names that you want to use. If None, the default color scale will be used (["lightblue", "lightyellow", "red"]).
    num_cols : int, optional
        Number of columns in the figure panel. It depends on the number of cell types you want to visualize.
    point_size : float, optional
        Size of the points in the plot.

    Returns
    -------
    plt
        Returns a plt object.

    Examples
    --------
    CARD_visualize_prop(proportion_df, spatial_location_df, ['TypeA', 'TypeB'], colors=["blue", "yellow"], num_cols=2, point_size=4.0)
    """
    # Check if the colors parameter is provided. If not, use the default colors.
    if colors is None:
        colors = ["lightblue", "lightyellow", "red"]
    else:
        colors = colors

    # Create dataframes (res_CARD and location) from the input proportion and spatial_location.
    res_CARD = pd.DataFrame(proportion)
    res_CARD = res_CARD[sorted(res_CARD.columns, key=str.casefold)]
    
    location = spatial_location.copy()
    
    # using 10.0x10.0 format
    loc_index = location['x'].astype(str) + 'x' + location['y'].astype(str)
    # using 10x10 format (integer format)
    #loc_index = location['x'].astype(int).astype(str) + 'x' + location['y'].astype(int).astype(str)
    
    # Check if the row indices of res_CARD and location match. If not, raise a ValueError.
    if not (res_CARD.index == loc_index).all():
        raise ValueError("The rownames of proportion data do not match with the rownames of spatial location data")

    # Scale the data
    res_CARD_scale = (res_CARD - res_CARD.min()) / (res_CARD.max() - res_CARD.min())
    res_CARD_scale['x'] = location['x'].values[:len(res_CARD_scale)]
    res_CARD_scale['y'] = location['y'].values[:len(res_CARD_scale)]
    
    # Calculate the number of rows and columns for the grid
    num_plots = len(ct_visualize)
    num_rows = min(4, (num_plots + 3) // 4)  # Ensure a maximum of 4 rows
    num_cols = min(4, num_plots)  # Maximum of 4 columns
    
    # Create a plot for each cell type
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(4*num_cols, 3*num_rows), sharex='col', sharey='row', squeeze=False, gridspec_kw={'hspace': 0.16, 'wspace': 0.1})
        
    custom_cmap = create_custom_colormap(colors)
    
    for i, cell_type in enumerate(ct_visualize):
        # Calculate the row and column index for the current plot
        row_index = i // num_cols
        col_index = i % num_cols
        
        # Select data for the current cell type
        cell_data = res_CARD_scale[['x', 'y', cell_type]]
        
        # Plotting
        ax = axes[row_index, col_index]
        sc = ax.scatter(cell_data['x'], cell_data['y'], c=cell_data[cell_type], cmap=custom_cmap, s=point_size)
        
        # Set plot properties
        ax.set(xticks=[], yticks=[], xlabel=None, ylabel=None, title=cell_type)
    
    # Create a colorbar legend for all plots
    cbar_ax = fig.add_axes([0.875, 0.4, 0.01, 0.2])  # [left, bottom, width, height]
    cbar = fig.colorbar(sc, cax=cbar_ax, orientation='vertical', label='value')
    
    # Remove grid lines, edge lines, and tick lines from colorbar legend
    cbar.ax.yaxis.grid(False)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(axis='y', which='both', length=0)
    # Set font size of tick labels at colorbar legend
    cbar.ax.tick_params(axis='y', labelsize=10)
            
    # Center the plot and legend in the canvas
    fig.subplots_adjust(left=0.1, right=0.85, top=0.9, bottom=0.1)
    
    return fig

def CARD_visualize_prop_2CT(proportion, spatial_location, ct2_visualize, colors=None):
    """
    Visualize the spatial distribution of two cell type proportions on the same plot.

    Parameters
    ----------
    proportion : pd.DataFrame
        Data frame, cell type proportion estimated by CARD in either original resolution or enhanced resolution.
    spatial_location : pd.DataFrame
        Data frame, spatial location information.
    ct2_visualize : list
        Vector of selected two cell type names that are interested in visualizing. Here we only focus on two cell types.
    colors : list, optional
        List of color names that you want to use for each cell type. 
        If None, the default color scale list will be used: 
        list(["lightblue", "lightyellow", "red"], ["lightblue", "lightyellow", "black"]).

    Returns
    -------
    plt
        Returns a plt object.

    Examples
    --------
    CARD_visualize_prop_2CT(proportion_df, spatial_location_df, ['TypeA', 'TypeB'], colors=["blue", "yellow"])
    """
    # Check if the colors parameter is provided. If not, use default color schemes.
    if colors is None:
        colors = [
            ["lightblue", "lightyellow", "red"],
            ["lightblue", "lightyellow", "black"]
        ]
    else:
        colors = colors
    
    # Create dataframes (res_CARD and location) from the input proportion and spatial_location.
    res_CARD = pd.DataFrame(proportion)
    res_CARD = res_CARD[sorted(res_CARD.columns)]
    
    location = spatial_location.copy()
    
    # using 10.0x10.0 format
    loc_index = location['x'].astype(str) + 'x' + location['y'].astype(str)
    # using 10x10 format (integer format)
    #loc_index = location['x'].astype(int).astype(str) + 'x' + location['y'].astype(int).astype(str)
    
    # Check if the row indices of res_CARD and location match. If not, raise a ValueError.
    if not (res_CARD.index == loc_index).all():
        raise ValueError("The rownames of proportion data do not match with the rownames of spatial location data")

    # Scale the data
    res_CARD_scale = (res_CARD - res_CARD.min()) / (res_CARD.max() - res_CARD.min())
    res_CARD_scale['x'] = location['x'].values[:len(res_CARD_scale)]
    res_CARD_scale['y'] = location['y'].values[:len(res_CARD_scale)]
    
    # Melt the dataframe
    mData = pd.melt(res_CARD_scale, id_vars=["x", "y"], var_name="Cell_Type")
    mData.rename(columns={"value": "Value"}, inplace=True)
    
    # Create custom colormaps
    cmap1 = create_custom_colormap(colors[0])
    cmap2 = create_custom_colormap(colors[1])

    # Create the plot
    fig, ax = plt.subplots(figsize=(8, 6))
 
    for i, ct in enumerate(ct2_visualize):
        subset_data = mData[mData["Cell_Type"] == ct]
        if i == 0:
            ax.scatter(subset_data["x"], subset_data["y"], c='none', edgecolors=cmap1(subset_data["Value"]), linewidth=0.9, s=115, label=ct)
            cbar_ax = fig.add_axes([0.8, 0.25, 0.02, 0.17])  # [left, bottom, width, height]
            cbar = fig.colorbar(cm.ScalarMappable(norm=None, cmap=cmap1), cax=cbar_ax, orientation='vertical')
            cbar.set_label(ct)
            # Remove grid lines, edge lines, and tick lines from colorbar legend
            cbar.ax.yaxis.grid(False)
            cbar.outline.set_visible(False)
            cbar.ax.tick_params(axis='y', which='both', length=0)
            # Set font size of tick labels at colorbar legend
            cbar.ax.tick_params(axis='y', labelsize=8)
        else:
            ax.scatter(subset_data["x"], subset_data["y"], c=subset_data["Value"], cmap=cmap2, marker='s', s=12, label=ct)
            cbar_ax = fig.add_axes([0.8, 0.55, 0.02, 0.17])  # [left, bottom, width, height]
            cbar = fig.colorbar(cm.ScalarMappable(norm=None, cmap=cmap2), cax=cbar_ax, orientation='vertical')
            cbar.set_label(ct)
            # Remove grid lines, edge lines, and tick lines from colorbar legend
            cbar.ax.yaxis.grid(False)
            cbar.outline.set_visible(False)
            cbar.ax.tick_params(axis='y', which='both', length=0)
            # Set font size of tick labels at colorbar legend
            cbar.ax.tick_params(axis='y', labelsize=8)
            
    # Set plot properties
    ax.set_aspect('equal')
    ax.set(xticks=[], yticks=[], xlabel=None, ylabel=None)
    
    # Center the plot and legend in the canvas
    fig.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.1)
    
    return fig


def CARD_visualize_pie(proportion, spatial_location, colors=None, radius=None, seed=None, sort=None, fontsize=None, legend_col=None):
    """
    Visualize the spatial distribution of cell type proportion in a geom scatterpie plot.

    Parameters
    ----------
    proportion : pd.DataFrame
        Data frame, cell type proportion estimated by CARD in either original resolution or enhanced resolution.
    spatial_location : pd.DataFrame
        Data frame, spatial location information.
    colors : list, optional
        Vector of color names that you want to use. If None, the color palette "Spectral" from the RColorBrewer package will be used.
    radius : float, optional
        Numeric value about the radius of each pie chart. If None, it will be calculated inside the function.
    seed : int, optional
        Seed number for generating colors if users do not provide the colors. If None, it will be generated inside the function.
    sort : boolean, optional
        sort by columns if it's not False 
    legend_col : int, optional
        number of legend columns
        
    Returns
    -------
    plt
        Returns a plt object.

    Examples
    --------
    CARD_visualize_pie(proportion, spatial_location, colors=["blue", "yellow", "red"], radius=0.02, seed=42, legend_col=5)
    """
    # read the input proportion and spatial_location. 
    # sort columns of res_CARD alphabetically if sort flag isn't set to False.
    res_CARD = proportion
    if sort is not False:
        res_CARD = res_CARD[sorted(res_CARD.columns, key=str.casefold)]
    
    location = spatial_location.copy()

    # using 10.0x10.0 format
    loc_index = spatial_location['x'].astype(str) + 'x' + spatial_location['y'].astype(str)
    # using 10x10 format (integer format)
    #loc_index = spatial_location['x'].astype(int).astype(str) + 'x' + spatial_location['y'].astype(int).astype(str)
    
    # Check if the row indices of res_CARD and location match. If not, raise a ValueError.
    if not (res_CARD.index == loc_index).all():
        raise ValueError("The rownames of proportion data do not match with the rownames of spatial location data")

    # Define color options: If colors is not provided, generate colors randomly or using a seed if specified. 
    # If provided, use the given colors.
    color_candidate = [
        "#1e77b4", "#ff7d0b", "#ceaaa3", "#2c9f2c", "#babc22", "#d52828", "#9267bc",
        "#8b544c", "#e277c1", "#d42728", "#adc6e8", "#97df89", "#fe9795", "#4381bd",
        "#f2941f", "#5aa43a", "#cc4d2e", "#9f83c8", "#91675a", "#da8ec8", "#929292",
        "#c3c237", "#b4e0ea", "#bacceb", "#f7c685", "#dcf0d0", "#f4a99f", "#c8bad8",
        "#F56867", "#FEB915", "#C798EE", "#59BE86", "#7495D3", "#D1D1D1", "#6D1A9C",
        "#15821E", "#3A84E6", "#997273", "#787878", "#DB4C6C", "#9E7A7A", "#554236",
        "#AF5F3C", "#93796C", "#F9BD3F", "#DAB370", "#877F6C", "#268785", "#f4f1de",
        "#e07a5f", "#3d405b", "#81b29a", "#f2cc8f", "#a8dadc", "#f1faee", "#f08080"
    ]
    
    if colors is None:
        if res_CARD.shape[1] > len(color_candidate):
            colors = color_candidate
        else:
            iseed = seed if seed is not None else 12345
            random.seed(iseed)
            colors = random.sample(color_candidate, res_CARD.shape[1])
    else:
        colors = colors

    # Set scaling factor for radius 
    scaling_factor = 145

    if radius is None:
        radius = ((location['x'].max() - location['x'].min()) * (location['y'].max() - location['y'].min())) / len(location)
        radius = np.sqrt(radius) * scaling_factor
    else:
        radius = radius * scaling_factor

    #### Design Canvas for main plot and legend
    ### (Step 1) Draw the main plot
    fig, axes = plt.subplots(figsize=(8,8))

    # Iterate through each index in 'location' and create pie chart using draw_pie
    for index in location.index:
        row = location.loc[index]
        x, y = row['x'], row['y']
        dist = res_CARD.loc[index].values
        draw_pie(dist, x, y, size=radius, colors=colors, ax=axes)
        
    # Set plot properties
    axes.axis('off')

    # Theme elements
    axes.set_facecolor('none')  # Panel Background
    axes.xaxis.label.set_visible(False)  # Axis Title
    axes.yaxis.label.set_visible(False)
    axes.tick_params(axis='both', which='both', length=0)  # Axis Ticks

    ### (Step 2) Draw legend below the main plot
    # Calculate the number of rows and columns for the legend
    legend_labels = list(res_CARD.columns)
    legend_handles = [Line2D([0], [0], marker='s', linestyle='None', markersize=10, markerfacecolor=colors[i]) for i in range(len(colors))]
    
    if legend_col is None:
        legend_col = 5  # Maximum 5 columns per a row
    num_columns = min(len(legend_labels), legend_col)
    
    if fontsize is None:
        fontsize = 8.2  # Default legend font size
    # Adjust the height ratio of the legend relative to the main plot
    legend_height_ratio = 0.05  # smaller value reduce the space between main plot and legend
    bbox_to_anchor_height = -legend_height_ratio / (1 + legend_height_ratio)

    # Add legend outside the plot in a horizontal way
    axes.legend(legend_handles, legend_labels, bbox_to_anchor=(0.5, bbox_to_anchor_height), 
                title="Cell Type", loc='upper center', frameon=False, 
                title_fontsize=12, ncol=num_columns, columnspacing=0.5, 
                handletextpad=0.5, handlelength=1, fontsize=fontsize)
    
    # Adjust layout and show the plot
    plt.tight_layout()

    return fig, axes


def CARD_visualize_Cor(proportion, colors=["#91a28c", "white", "#8f2c37"]):
    """
    Visualize the cell type proportion correlation.

    Parameters
    ----------
    proportion : pd.DataFrame
        Data frame, cell type proportion estimated by CARD in either original resolution or enhanced resolution.
    colors : list, optional
        Vector of color names that you want to use. If None, the default color scale will be used (["#91a28c", "white", "#8f2c37"]).

    Returns
    -------
    ggplot
        Returns a ggplot object.

    Examples
    --------
    CARD_visualize_Cor(proportion_df, colors=["blue", "white", "red"])
    """
    # Sort columns without case sensitivity
    proportion = proportion.reindex(sorted(proportion.columns, key=str.casefold), axis=1)
    
    # Compute correlation matrix
    cor_CARD = proportion.corr()
    
    # Default colors if not provided
    if colors is None:
        colors = ["#91a28c", "white", "#8f2c37"]
    
     # Normalize correlation values to the range [-1, 1]
    norm = Normalize(vmin=-1.0, vmax=1.0)
    
    # Create custom colormap
    cmap_segments = [(i / (len(colors) - 1), colors[i]) for i in range(len(colors))]
    custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", cmap_segments)
    
    # Create the correlation plot
    fig, ax = plt.subplots(figsize=(16, 14))
    custom_cmap = create_custom_colormap(colors)
    cax = ax.matshow(cor_CARD, cmap=custom_cmap, aspect='auto', norm=norm)
    
    # Set up the colorbar
    cbar_ax = fig.add_axes([0.875, 0.36, 0.015, 0.17])  # [left, bottom, width, height]
    cbar = fig.colorbar(cax, cax=cbar_ax, orientation='vertical')
    cbar.set_label('Corr', fontsize=12)  # Change the fontsize parameter as needed
    
    # Set colorbar range
    cbar.set_ticks([-1.0, 0.0, 1.0])
    
    # Remove grid lines, edge lines, and tick lines from colorbar legend
    cbar.ax.yaxis.grid(False)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(axis='y', which='both', length=0)
    # Set font size of tick labels at colorbar legend
    cbar.ax.tick_params(axis='y', labelsize=12)
    
    # Set plot properties
    ax.set_title("Correlation")
    ax.set_xticks(range(len(cor_CARD.columns)))
    ax.set_yticks(range(len(cor_CARD.columns)))
    ax.set_xticklabels(cor_CARD.columns, rotation=90)
    ax.set_yticklabels(cor_CARD.columns)
    
    # Remove grid lines, edge lines, and tick lines from colorbar legend
    ax.tick_params(axis='both', which='both', length=0) # Remove tick lines
    ax.grid(False)                                      # Remove grid lines
    # Remove all edge lines from the plot
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
     # Set font size of tick labels
    ax.tick_params(axis='both', labelsize=12)

    # Center the plot and legend in the canvas
    fig.subplots_adjust(left=0.2, right=0.85, top=0.80, bottom=0.05)
    
    return fig


def plot_location_imputation(origin_location, location_imputation):
    """
    Visualize the newly grided spatial locations to see if the shape is correctly detected. 
    
    Parameters
    ----------
    origin_location : pd.DataFrame
        Data frame, original spatial gene location.
    
    location_imputation : pd.DataFrame
        Data frame, spatial gene location newly added for an enhanced resolution.
    
    Returns
    -------
    ggplot
        Returns a ggplot object.
        
    location_imputation
        Return a pd.DataFrame storing the x, y coordinates
        
    Examples
    --------
    fig, location_imputation = plot_location_imputation(origin_location, location_imputation)
    """
    # Creating the scatter plot
    fig = plt.figure(figsize=(8, 6))
    plt.scatter(origin_location['x'], origin_location['y'], marker='o', s=6, color="#7dc7f5", label='Original Points')
    plt.scatter(location_imputation['x'], location_imputation['y'], marker='o', s=3, color='red', label='Imputation Points')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Location Imputation')
    plt.grid(False)
    plt.legend(fontsize='medium')

    return fig


def CARD_visualize_gene(spatial_expression, spatial_location, gene_visualize, point_size=None, colors=None, num_cols=None):
    """
    Visualize the marker gene expression.

    Parameters
    ----------
    spatial_expression : pd.DataFrame
        Data frame, spatial gene expression in either original resolution or enhanced resolution.
    spatial_location : pd.DataFrame
        Data frame, spatial location information.
    gene_visualize : list
        Vector of selected gene names that are interested in visualizing.
    colors : list, optional
        Vector of color names that you want to use. If None, the default color scale in the virdis palette will be used.
    NumCols : int, optional
        Numeric, number of columns in the figure panel. It depends on the number of cell types you want to visualize.

    Returns
    -------
    ggplot
        Returns a ggplot object.

    Examples
    --------
    CARD_visualize_gene(expression_df, location_df, ['GeneA', 'GeneB'], colors=["blue", "yellow", "red"], NumCols=2)
    """
    expression = spatial_expression.div(spatial_expression.sum(axis=1), axis=0)  # Normalize expression data
    location = spatial_location.copy()
    
    # using 10.0x10.0 format
    loc_index = location['x'].astype(str) + 'x' + location['y'].astype(str)
    # using 10x10 format (integer format)
    #loc_index = location['x'].astype(int).astype(str) + 'x' + location['y'].astype(int).astype(str)  
    
    if not set(expression.columns).issubset(loc_index):
        raise ValueError("The column names of the expression data do not match the row names of the spatial location data")
    
    gene_select = gene_visualize
    if sum(gene.upper() in map(str.upper, expression.index) for gene in gene_select) != len(gene_select):
        raise ValueError("There exist selected genes that are not in the expression data!")

    # Create a new dataframe (data) for plotting. Normalize the expression values for each gene and concatenate the data.
    data = pd.DataFrame()  # Initialize data
    for gene in gene_select:
        ind = np.where(expression.index.str.upper() == gene.upper())[0][0]
        df = pd.Series(expression.iloc[ind], name='value')
        df.index = expression.columns
        df = (df - df.min()) / (df.max() - df.min())
        
        d = pd.DataFrame({'value': df})
        d['x'] = location['x'].values[:len(d)]
        d['y'] = location['y'].values[:len(d)]
    
        d['gene'] = gene
        data = pd.concat([data, d], axis=0)  # Concatenate to data
        
    data['gene'] = pd.Categorical(data['gene'], categories=gene_select)
    
    # default num_cols (if it's None) is same as the number of gene_select
    num_cols = len(gene_select) if num_cols is None else num_cols
    
    # Create the FacetGrid
    fig = sns.FacetGrid(data, col="gene", col_wrap=num_cols)
    
    # Plot your data using fig.map_dataframe
    palette = "viridis" if colors is None else None         # default colors = viridis palette
    point_size = 50.0 if point_size is None else point_size # default point_size = 50
    fig.map_dataframe(sns.scatterplot, x="x", y="y", hue="value", palette=palette, marker="s", s=point_size, edgecolor="none")

    # Set x-axis and y-axis attributes
    for ax in fig.axes.flat:
        ax.set(xticks=[], yticks=[], xlabel=None, ylabel=None, title=gene)
        ax.spines['bottom'].set_visible(False)  # Turn off bottom spine (x-axis)
        ax.spines['left'].set_visible(False)    # Turn off left spine (y-axis)

    # Set the figure title and font size
    fig.set_titles("{col_name}", size=15)
    
    # control spacing between subplots
    plt.subplots_adjust(wspace=0.1, hspace=0.5)

    # Create a colorbar legend for all plots
    cbar_ax = plt.gcf().add_axes([0.5, 0.1, 0.1, 0.05])  # Adjust the position and size of the colorbar
    cbar = plt.colorbar(fig.facet_axis(0, 0).collections[0], cax=cbar_ax, orientation='horizontal')
    # Position the label at the left side of the colorbar
    cbar.ax.text(-4.5/num_cols, 0.9/num_cols, 'Expression', ha='left', va='center', fontsize=15)
    # Set the font size of the tick labels
    cbar.ax.tick_params(labelsize=10)
    
    return fig

def CARD_visualize_SCMapping(MapCellCords, colors, fontsize=None, legend_col=None):
    """
    Visualize the single cell resolution gene expression. 
    
    Parameters
    ----------
    MapCellCords : pd.DataFrame
        Data frame, mapped single cell spatial locations 
    
    colors : list
        Vector of color names that you want to use.
    
    fontsize : int
        integer, legend font size 
    
    legend_col : int
        integer, maximum number of columns per a legend row
        
    Returns
    -------
    ggplot
        Returns a ggplot object.
          
    Examples
    --------
    CARD_visualize_SCMapping(MapCellCords, colors, legend_col=5)
    
    """
    # extract the number of unique cell types
    unique_cells = sorted(MapCellCords['CT'].unique(), key=lambda x: x.lower())
    num_unique_cell = len(unique_cells)
    
    if num_unique_cell > len(colors):
        raise ValueError("Insufficient colors. Please add more colors.")
        
    #### (Step 1) Design Canvas for main plot and legend
    fig, axes = plt.subplots(figsize=(8,8))
    plt.scatter(MapCellCords['x'], MapCellCords['y'], c=[colors[unique_cells.index(ct)] for ct in MapCellCords['CT']], s=10)

    # Set plot properties
    axes.axis('off')

    # Theme elements
    axes.set_facecolor('none')  # Panel Background
    axes.xaxis.label.set_visible(False)  # Axis Title
    axes.yaxis.label.set_visible(False)
    axes.tick_params(axis='both', which='both', length=0)  # Axis Ticks
	
	#### (Step 2) Draw legend below the main plot
	# Calculate the number of rows and columns for the legend
    legend_labels = unique_cells
    legend_handles = [Line2D([0], [0], marker='o', linestyle='None', 
                        markersize=5, markerfacecolor=color, markeredgecolor='none') for color in colors]

    if legend_col is None:
        legend_col = 5  # Maximum 5 columns per a row
    num_columns = min(len(legend_labels), legend_col)
    
    if fontsize is None:
        fontsize = 8.2  # Default legend font size
    
    # Adjust the height ratio of the legend relative to the main plot
    legend_height_ratio = 0.05  # smaller value reduce the space between main plot and legend
    bbox_to_anchor_height = -legend_height_ratio / (1 + legend_height_ratio)
	
	# Add legend outside the plot in a horizontal way
    axes.legend(legend_handles, legend_labels, bbox_to_anchor=(0.5, bbox_to_anchor_height), 
                title="Cell Type", loc='upper center', frameon=False, 
                title_fontsize=12, ncol=num_columns, columnspacing=0.5, 
                handletextpad=0.5, handlelength=1, fontsize=fontsize)
    
    # Adjust layout and show the plot
    plt.tight_layout()

    return fig, axes
    

