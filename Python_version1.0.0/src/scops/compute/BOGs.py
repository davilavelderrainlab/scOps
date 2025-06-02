# package specific function
import scanpy as sc
from scops.helper.FilterBOGs import FilterBOGs
from scops.helper.Dictionary2Memership import Dictionary2Memership
from scops.helper._panda_in_var_ import _panda_in_var_

# BOGs
def BOGs(adata,
    group,
    layer='logcounts',
    method='wilcoxon',
    correction='benjamini-hochberg',
    foldchange_th=0.05,
    adjpvalue_th=0.01,
    expression_th=0.25,
    sort_field='pvals_adj',
    return_dea_or_bog='both',
    field_name='BOGs',
    DE_out_field='DE_out',
    BOGDict_field='BOGDict',
    verbose=True,
    warningsOff=False,
    **kwargs):
    """Performs differentially expressed analysis (DEA), 1vsAll, for each group
    
    It returns the DEA statistics of differentially expressed genes (BOGs) in 
    .uns of anndata under ``'DE_out_field'``, also return in .uns the dictionary 
    of BOGs names under ``'BOGDict_field'``. 

    BOGs are stored in adata.varm under the field_name

    It uses sc.tl.rank_genes_groups to perform DEA with different defautls. 
    
    Parameters
    ----------
    adata       
        AnnData object: stores your data 
    group 
        str: column of .obs or array of .obsm, specifies the groups used for DEA 
    layer
        str='logcounts': layer of adata.layer. specifies the counts to use for 
            DEA
    method
        str='wilcoxon': method to be used to perform DEA. more at 
            sc.tl.rank_genes_groups
    correction
        str='benjamini-hochberg': correction of p-value for multiple testing. 
            more at sc.tl.rank_genes_groups
    foldchange_th
        float=0.05: FoldChange threshold for BOGs
    adjpvalue_th
        float=0.01: Adjusted p-value threshold for BOGs
    expression_th=0.25
        float: percentage of cells expressing the gene in the group for BOGs
    sort_field
        str='pvals_adj': field used to order BOGs
    return_dea_or_bog
        str='both': whether to return DEA statistics or BOGs as a result of DEA 
    field_name
        str=None: name used to store BOGs matrix
    DE_out_field
        str='DE_out': .uns field used to store DEA statistics of BOGs 
    BOGDict_field
        str='BOGDict': .uns field used to store dictionary of BOGs 
    verbose
        bool=True: wheter to have info from the function while it is running 
    warningsOff
        bool=False: wheter to soppress warnings
    **kwargs
        other named arguments to be passed to sc.tl.rank_genes_groups

    Returns
    -------
    adata  
        with new .uns fields and BOGs membership vectors in .var

    Examples
    --------
    adata = anndata.read_h5ad(path_to_anndata)
    adata = computeBOGs(adata, group = 'leiden')
    
    Notes
    -----
    Performance warnings will arise if BOGs are computed on many groups, due to 
    scanpy functions. To ignore them, set warnignsOff = True.
    """
    
    decomposition = list(adata.obsm.keys())
    partition = list(adata.obs.columns)

    if warningsOff : 
        import warnings
        warnings.filterwarnings("ignore")

    if group not in (partition + decomposition):
        raise(ValueError(f'{group} is not in adata.obs, please add it to adata.obs'))

    if layer not in adata.layers :
        raise(ValueError(f'{layer} not in adata, layers, please specifiy a valid layer'))

    # to handle matrices 
    if group in decomposition :
        if verbose : 
            print(f'computing membership vector from {group}')
        # check that all the groups have at least one cell, to avoid breaking DEA 
        array = adata.obsm[group]
        membership_vector = array.argmax(axis = 1)
        group = group + '_membership'
        adata.obs[group] = membership_vector

        # cheking for rare groups which will break DEA 
        table = adata.obs[group].value_counts()
        # prototypes with single-cell assigned 
        rare_group = table[table < 2].index.to_list()
        position = [obs in rare_group for obs in adata.obs[group]]
        if len(rare_group) != 0 : 
            rare_group = ', '.join(map(str, rare_group))
            print(f'single obseration groups from groups: {rare_group} have been aggregated to perform DEA.') 
        adata.obs.loc[position, group] = 'rare'

        # this is stupid but it does not break: converting to pandas 
        adata.obs[group] = adata.obs[group].astype(str).astype('category')

    # Differential Expression Anlysis 
    if verbose : 
        print(f'performing differential expression analysis')
    sc.tl.rank_genes_groups(adata, groupby = group, layer = layer, method = method, corr_method = correction, pts = True, use_raw = False, **kwargs)
    # do not specify group using a column of obs, since it works to isolate a single group! if you want all of them use None instead
    bog = sc.get.rank_genes_groups_df(adata, group = None) 

    # Filter DEA results 
    if verbose : 
        print(f'filtering out not differentially expressed genes')
    bog_dict, dea_dict = FilterBOGs(bog,
        foldchange_th = foldchange_th, 
        adjpvalue_th = adjpvalue_th,
        expression_th = expression_th,
        sort_field = sort_field, 
        return_dea_or_bog = return_dea_or_bog)

    # Membership Matrix 
    if verbose : 
        print(f'computing BOGs memberships')
    unique_groups = list(adata.obs[group].unique())
    features = adata.var_names.to_list()
    groups = unique_groups
    membership = Dictionary2Memership(bog_dict, features = features, groups = groups) 
        
    # Storing 
    
    # deprecated
    # description: For each group a BOGs membership vector is present in .var under ``'BOG_group_name'``. 
    #adata = _panda_in_var_(adata, membership, field_name = 'BOGs_', group_names = group_names)
    
    adata.varm[field_name] = membership
    adata.uns[DE_out_field] = dea_dict
    adata.uns[BOGDict_field] = bog_dict 

    return(adata)