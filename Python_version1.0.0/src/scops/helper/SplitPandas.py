def SplitPandas(df,
    groups, 
    axis=1): 
    """Splits a DataFrame based on a column or o row, default is for column
    
    It works on columns for axis = 1 or on rows for axis = 0. 

    Parameters
    ----------
    df
        pandas DataFrame
    groups
        column or row used to split the data frame
    axis
        int and on of:
        * ``'0'`` -, split based on rows
        * ``'1'`` -, split based on columns
    
    Returns
    -------
    A dictionary where keys are the groups and values are the actual splits. 

    Notes
    -----
    none
    """

    if axis not in [0, 1] : 
        raise(ValueError('please pick an axis between 0 : rows and 1 : columns'))

    # col split 
    if axis == 1 : 
        if groups not in df.columns.values : 
            raise(ValueError(f'{group} not in data frame column names'))

    # row split 
    else : 
        # to make it work with column splits 
        df = df.T 

    # split 
    unique_groups = df[groups].unique()
    df_split = {} # named list == dictionary 
    for group in unique_groups : 
        msk = df[groups] == group
        split_i = df.loc[msk]

        # adjusting for row split 
        if axis == 0 : 
            split_i = split_i.T

        df_split[group] = split_i

    return(df_split)