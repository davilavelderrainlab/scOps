def _panda_in_var_(adata,
    data,
    field_name='BOG_',
    group_names=None) : 
    """Handy function to store pandas data frames in .var 
    
    For each column of the data frame, a new .var colum will be added 
    with a specified name given by field_name + group_name

    Parameters
    ----------
    adata       
        AnnData object: stores your data 
    data
        pandas DataFrame: which will be stored in .var
    field_name
        str='BOG_': name to identify the new columns of .var
    group_names
        list=None: name used to identify BOGs results for each group in .var 
            * ``'None'`` -, uses the deafult names of the groups found in groups
            * ``'List'`` -, needs to have one str entry for each group

    Returns
    -------
    AnnData object with new .var fields 

    Notes
    -----
    none
    """

    unique_groups = data.columns.to_list()

    # storing in the anndata 
    if group_names is None : 
        group_names = unique_groups

    for group, group_name in zip(unique_groups, group_names): 
        field = field_name + group_name
        adata.var[field] = data.loc[:, group]
        
    return(adata)