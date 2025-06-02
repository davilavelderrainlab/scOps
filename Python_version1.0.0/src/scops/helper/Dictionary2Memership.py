import pandas as pd

def Dictionary2Memership(dictionary,
    features,
    groups): 
    """Defines a membership matrix from a dictionary 
    
    It needs the features of the membership matrix and the groups(rows and cols)

    Parameters
    ----------
    dictionary       
        dictionary: each group is a key and features are the values
    features
        list: list of features (rows of the memebership matrix)
    groups
        list: list of groups (columns of the memebership matrix)

    Returns
    -------
    A binary pandas DataFrame (features x groups) specifying membership. 

    Notes
    -----
    none
    """
    
    # initialisation of the martix
    membership = pd.DataFrame(
        data = 0,
        index = features,
        columns = groups,
        dtype = 'int'
    )

    # ele_dict: element dictionary 
    for group, ele_dict in dictionary.items():
        position = [feature for feature in ele_dict]
        membership.loc[position, group] = 1

    return(membership)