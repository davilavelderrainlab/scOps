# package specific function 
import pandas as pd 
from scops.helper._gettoporbottom_ import _gettoporbottom_

def GetTopOrBottom(df,
    column,
    n=5,
    bottom=False) : 
    """Gets the top(bottom) elements of a DataFrame/Dictionary of PD/List of PD
    PD: Pandas Dataframes

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
    Returns the top(bottom) n entries of a pandas DataFrame.
    For dictionaries and lists returns a list. 

    Notes
    -----
    none
    """
    
    # single pandas data frame 
    if isinstance(df, pd.core.frame.DataFrame) : 
        res = _gettoporbottom_(df, column = column, n = n, bottom = bottom)

    # multiple pandas data frames 
    #       as a dictionary
    elif isinstance(df, dict) : 
        df = [values for values in df.values()]

    #  as a list
    if isinstance(df, list) : 
        res = []
        for split in df :
            top_bottom = _gettoporbottom_(split, column = column, n = n, bottom = bottom)
            res.append(top_bottom)
        res = pd.concat(res)

    return(res)
