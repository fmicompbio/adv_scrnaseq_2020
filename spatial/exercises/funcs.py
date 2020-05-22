#!/usr/bin/env python3



import anndata as ad
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from typing import Tuple,Union, List,Dict



def visualize_score_landscape():
    from mpl_toolkits.mplot3d import Axes3D 
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter


    fig = plt.figure(figsize=(10,10))
    ax = fig.gca(projection='3d')

    # Make data.
    X = np.arange(0, 1, 0.01)
    Y = np.arange(0, 1, 0.01)
    X, Y = np.meshgrid(X, Y)
    Z = X*Y
    # Plot the surface.
    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                        linewidth=0, antialiased=False)

    # Customize the z axis.
    ax.set_zlim(-1.01, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    ax.set_zlabel("TLS_score",fontsize = 15)
    ax.set_xlabel("B-cell proportion",fontsize = 15)
    ax.set_ylabel("T-cell proportion",fontsize = 15)


    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()

def join_data(data : Dict[str,ad.AnnData],
              )-> ad.AnnData:

    diter = iter(data.keys())

    uns = {n:{} for n in range(len(data))}

    for n,d in enumerate(data.values()):
        for k,v in d.uns.items():
            uns[n].update({k:v}) 

    k0 = next(iter(data.keys()))
    d0 = data.pop(k0)


    out = d0.concatenate(*list(data.values()),
                         join = "outer",
                         batch_key = "sample",
                         index_unique = "_")
    out.uns = uns

    return out



from scipy.stats import mannwhitneyu as utest
from statsmodels.stats.multitest import multipletests


def dge_test(data : ad.AnnData,
             labels : np.ndarray,
             contrast : Union[str,int],
             alpha : float = 0.01,
             positive : bool = True,
            )-> pd.DataFrame:
    
    """ conduct dge test
    
    Parameters:
    -----------
    
    data : ad.AnnData
        anndata object containing spatial data
    labels : np.ndarray
        labels indicating which group each observation
        belongs to. E.g., cluster indices.
    contrast : Union[str,int]
        label to be compared to the rest. Must
        be a member of "labels"
    alpha : float
        family-wise error rate (FWER) used for
        correction of MHT
    positive : bool
        only return genes with a positive log
        fold change. 
    
    Returns:
    --------
    Pandas data frame with 4 columns:
    pval, contains unadjusted p-values from
    the U test; pval_adj contains the adjusted
    p-values; l2fc, the log2FoldChange values,
    when comparing the mean of the selected
    group toward the rest; sig, inidcator
    of whether the gene has a significant
    p-values (1) or not (0). Genes are
    sorted w.r.t significance.
    
    
    """
    
    # object to hold result values
    res = np.zeros((data.X.shape[1],4))
    res[:,0] = np.nan
    res = pd.DataFrame(res,
                       index = data\
                               .var["name"]\
                               .values,
                       columns = ['pval',
                                  "pval_adj",
                                  "l2fc",
                                  "sig"],
                      )
    
    
    # library size normalization 
    libsize = data.X.sum(axis =1,
                         keepdims = True)
    
    norm_data = data.X / libsize
    
    # get indices for group members
    pos_1 = labels == contrast
    pos_2 = labels != contrast
     
    # test every gene
    for g in range(data.X.shape[0]):
        # skip if gene is not observed in any spot
        if norm_data[:,g].sum() <= 0:
            res.iloc[g,0] = np.nan
        else:
            # perform U test
            _,pval = utest(norm_data[pos_1,g],
                           norm_data[pos_2,g])
            # store p-value
            res.iloc[g,0] = pval 

            # compute average expression for
            # selected cluster and rest
            mu_1 = np.mean(norm_data[pos_1,g])
            mu_2 = np.mean(norm_data[pos_2,g])

            # compute log2FoldChange
            if mu_1 == 0: 
                l2fc = -np.inf
            elif mu_2 == 0:
                l2fc = np.inf
            else:
                l2fc = np.log2(mu_1) - np.log2(mu_2)
            
            # sotre l2FC
            res.iloc[g,2] = l2fc
    
    # remove bad genes
    res = res.dropna()
    # correct fot MHT, adjust p-values
    mht = multipletests(pvals = res["pval"].values,
                                method = "fdr_bh",
                                alpha = alpha)
    
    # store adjusted p-values
    res.loc[:,"pval_adj"] = mht[1]
    # add indicator of significane
    res.loc[:,"sig"] = (res["pval_adj"]\
                        .values <= alpha)\
                        .astype(int)
    # sort values
    res = res.sort_values(by = "pval_adj",
                          axis = 0)
    
    # remove negative l2FC values
    # if specified
    if positive:
        res = res.\
              iloc[res["l2fc"]\
              .values > 0,:]
    
    return res


def plot_data(ax : plt.Axes,
              data : ad.AnnData,
              feature : Tuple[str,str] = None,
              title : str = None,
              plt_args : dict = None,
              show_image : bool = True,
              clean : bool = True,
              index : str = None,
             )->None:
    
    """ Plot Spatial Data
    
    Parameters:
    ----------
    ax : plt.Axes
        matplotlib axes object to plot in
    data : ad.AnnData
       anndata object containing spatial data
    feature : Tuple[str,str]
        first element the type of feature (e.g., 'name') to be plotted
        second element is variable name (e.g., 'ERBB2')
    title : str
        string containing title of plot
    plt_args : dict
        matplotlib keyword arguments for scatter plot
        to customize the visualization
    show_image : bool
        overlay plot on image
    clean :
        remove spines and tick from plot
    
    """
    
    # get scale factor

    if index is not None:
        uns = data.uns[index]
    else:
        uns = data.uns

    sf = uns['tissue_hires_scalef']
    
    # plot image as background
    # if specified
    if show_image:
        ax.imshow(uns['image_hires']\
                  .transpose(1,0,2))
    else:
        ax.invert_yaxis()
    
   # default plot settings
    _plt_args = dict(alpha = 0.8,
                     edgecolor = None,
                     marker = 'o',
                     s = 20,
                     cmap = plt.cm.magma,
                    )
    # update plot settings if
    # provided as input
    if plt_args is not None:
        for k,v in plt_args.items():
            _plt_args[k] = v
            
    # set feature values as color
    # values to be plotted
    if feature is not None:
        # locate position of 
        # feature in the data matrix
        pos = np.where(data.var[feature[0]]\
                       .values == feature[1])[0]
        
        # grab feature values
        vals = data.X[:,pos].flatten()
        vals /= vals.max()
        
        # set colorspace
        if plt_args["cmap"] == "feature":
            # if alpha-values should
            # be used to to display
            # signal strength of feature
            rgba = np.zeros((vals.shape[0],4))
            rgba[:,0] = 1
            rgba[:,3] = vals 
            _plt_args["c"] = rgba
            _plt_args["cmap"] = None
            _plt_args["alpha"] = None
        else:
            # if colormap should be
            # used to display signal
            # values
            _plt_args["c"] = vals
            
    
    # generate scatter plot
    # from spot coordinates
    # adjusted to fit image using
    # the scaling factor
    ax.scatter(data.obs['x'] * sf,
               data.obs['y'] * sf,
               **_plt_args,
              )
    
    ax.set_aspect("equal")
    
    # set title if one is provided
    if title is not None:
        ax.set_title(title)
        
    # clean axes if specified
    if clean:
        clean_ax(ax)
            
    return None

def clean_ax(ax : plt.Axes,
             remove_spines : Union[str,List[str]] = "all",
             keep_ticks : bool = False,
            )->None:
    """clean axes
    
    removes ticks and 
    spines from an axes
    object
    
    Parameters:
    ----------
    ax : plt.Axes
        axes object to clean
    
    remove_spines : Union[str,List[str]]
       list of spines to remove. If "all" is
       given, all spines will be removed.
    
    keep_ticks : bool
        set to true if tick should be kept
    
    """
    
    if remove_spines is None:
        remove_spines = []
    elif isinstance(remove_spines,str):
        remove_spines = [remove_spines]
        
    if remove_spines[0] == "all":
        remove_spines = ["left",
                         "top",
                         "right",
                         "bottom"]
    
    if not keep_ticks:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    for sp in remove_spines:
        ax.spines[sp].set_visible(False)
    
def naive_normalize(X : np.ndarray,
                    )->np.ndarray:
    
    Y = 2*np.sqrt(X + 3./8.)
    Y /= Y.sum(axis=1,
               keepdims = True)
    return Y
