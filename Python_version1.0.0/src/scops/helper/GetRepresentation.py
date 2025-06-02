
def GetRepresentation(adata, 
    field) :
    """Retrieves the specified representation fro .var in a pd DataFrame

    Parameters
    ----------
    adata
        AnnData object where you store your data 
    field 
        str: expression to be matched in .var. 
        All the columns matchin that pattern will be included in the final DataFrame
    
    Return
    ------
        data
            pandas data frame (features x representation)
    """
    data = adata.var[list([col for col in adata.var.columns if field in col])]

    return(data)