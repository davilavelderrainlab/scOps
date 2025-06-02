def _gettoporbottom_(df,
    column,
    n=5,
    bottom=False): 
    """subsets a pandas data frame and get the top (bottom) n elements
    
    works on a specified column 

    Parameters
    ----------
    df
        pandas DataFrame
    column
        column used to order the entries of the data frame
    n
       int=5: number of elements to select
    bottom
        bool=False: wheter to get the bottom or top entries

    Returns
    -------
    Top(Bottom) n entries of the data frame

    Notes
    -----
    helper function; please use GetTopOrBottom instead 
    """

    if column not in df.columns.values : 
        raise(ValueError(f'{column} not in column names of the data frame'))

    # picking the top or bottom n columns: 
    # NOTE: the order will be highest to lowest for top and inverted for bottom. 
    res = df.sort_values(column, ascending = bottom).iloc[0 : n] # it is already -1 on n
    # NOTE it is already properly ordered for both top and bottom! 
    return(res)