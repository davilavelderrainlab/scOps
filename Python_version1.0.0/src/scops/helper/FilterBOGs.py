def FilterBOGs(bog,
    foldchange_th=0.05,
    adjpvalue_th=0.01,
    expression_th=0.25,
    sort_field='pvals_adj',
    return_dea_or_bog='both') :
    """Filters the results of DEA and gives only differentially expressed genes 
    
    Already implemented in computeBOGs.

    Parameters
    ----------
    bog       
        pandas Data Frame: obtained from sc.get.rank_genes_groups_df
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
    return_dea_or_bog
        str='both': what to return; can be one of 
            * ``'both'`` -, returns both BOGs and DEA dictionaries 
            * ``'bog'`` -, returns only BOGs dictionary
            * ``'dea'`` -, returns only DEA dictionary

    Returns
    -------
    Depens on return_dea_or_boh; By default it returns both BOGs and DEA
    dictionaries, where each key is a group and values are the names of the BOGs
    (BOGs dictionary) or the full statistics of DEA (DEA dictionary),
    respectively. 

    Notes
    -----
    none
    """

    # return_dea_or_bog can be dea; bog or both 
    fc_bad = abs(bog['logfoldchanges']) < foldchange_th
    adjp_bad = bog['pvals_adj'] > adjpvalue_th
    exp_bad = bog['pct_nz_group'] < expression_th
    remove = fc_bad.add(adjp_bad)
    keep = ~remove.add(exp_bad)
    bog_good = bog.loc[keep]

    # storing the results of DEA seprately for each group 
    bog_dict = {}
    dea_dict = {}
    for group in bog_good['group'].unique() :
        # copy to avoid sorting a view instead of a hard copy 
        split = bog_good.loc[bog_good['group'] == group, :].copy()
        split.sort_values(sort_field, ascending = False, inplace = True)
        # what info do you want to store?
        if return_dea_or_bog == 'bog' :
            bog_dict[group] = split['names']
        elif return_dea_or_bog == 'dea': 
            bog_dict[group] = split
        elif return_dea_or_bog == 'both' : 
            bog_dict[group] = split['names']
            dea_dict[group] = split
        else : 
            raise(ValueError('please provide a return_dea_or_bog which is: dea, bog or both'))

    return(bog_dict, dea_dict)
