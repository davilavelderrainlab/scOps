from scipy.sparse import csr_matrix 
import pandas as pd 
import numpy as np
from anndata import AnnData

# package specific function


def Profiles(adata, 
    group,
    layer='logcounts',
    field_name='Profiles',
    return_panda=False): 
    """Computes Profiles from AnnData object by specifying a representation

    Profiles are feature-wise averages of the counts for each representation. 

    Representations can be either membership vectors in .obs or arrays 
    in .obsm. In the latter case, a weighted mean will be used to compute 
    profiles, to account for all the representations.

    Results are stored in adata.varm

    Parameters
    ----------
    adata       
        AnnData object or Matrix: stores your data or is the matrix to use. 
    group 
        str: .obs column or .obsms array, specifies the groups/representation
            used to group or weigth observations
        list: that specifies the a membership of each sample
        matrix: that specifies the contribution of each prototype. 
    layer
        str='logcounts': layer of adata.layer. specifies the counts used
            to compute profiles
    field_name
        str='Profile': name used to store BOGs matrix
    return_panda
        bool=False: wheter to return the pandas DataFrame of profiles instead of
            the AnnData object

    Returns
    -------
    AnnData  
        with new Profiles in .var
        OR a pandas DataFrame features x profiles if return_pandas=True.

    Examples
    --------
    adata = anndata.read_h5ad(path_to_anndata)
    adata = computeProfiles(adata, group = 'leiden')
    
    Notes
    -----
    none
    """
    
    # check for matrix or anndata
    is_anndata = isinstance(adata, AnnData)

    # weighted profiles (Decomposition Based Groups)
    prototypes = None
    # on the anndata
    if is_anndata and (group in adata.obsm) : 
        array = adata.layers[layer]
        prototypes = adata.obsm[group] 
    # on the matrix
    elif (not isinstance(group, str)) and group.ndim > 1 :
        array = adata 
        prototypes = group  
    
    # average profiles (Partition Based prototypes)
    else :
        if is_anndata :
            if group in adata.obs.columns :
                # working on the layer 
                array = adata.layers[layer]
                membership_vector = adata.obs[group]
        # work on the matrix directly
        else : 
            array = adata 
            membership_vector = group

    # average profiles (Partition Based prototypes)
    if prototypes is None :
        # horizontal split of the array
        counts_splits = {group_x : array[np.where(membership_vector == group_x)] for group_x in membership_vector.unique()}

        # computing feature-wise averages 
        profiles = {}
        for group_i, counts in counts_splits.items() :
            # dense single layer conversion. 
            profiles[group_i] = np.array(counts.mean(axis = 0)).reshape(-1)

        profiles = pd.DataFrame.from_dict(profiles)
    
        # add gene names from the anndaya
        if is_anndata :
            profiles = profiles.set_index(adata.var_names) 
    
    # Weighted Profiles based on continous prototypes
    if prototypes is not None:
        # scaling
        # scaling between 0 and 1 the contribution of each prototype to any given cell
        # matrix operation goes speed! 
        scaled_prototypes = prototypes.T @ np.diag(1 / prototypes.sum(axis = 1))

        # weight 
        weighted = csr_matrix(scaled_prototypes) @ array # both sparse 

        # mean 
        # ugly rearranging due to 
        denominator = csr_matrix(np.diag( 1 / scaled_prototypes.sum(axis = 1)))
        weighted_profiles = denominator @ weighted

        # proper formatting for Pandas DataFrame
        profiles = pd.DataFrame.sparse.from_spmatrix(weighted_profiles.T)
        profiles.columns = [str(profile) for profile in profiles.columns]
        
        if is_anndata :
            profiles.index = adata.var_names

    # directly return the pandas Data Frame 
    if return_panda or (not is_anndata): 
        return(profiles)

    else : 
        adata.varm[field_name] = profiles
        return(adata)