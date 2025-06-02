from scipy.sparse import csr_matrix 
import pandas as pd 
import numpy as np

# package specific function
from scops.helper._panda_in_var_ import _panda_in_var_

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
        AnnData object: stores your data 
    group 
        str: .obs column or .obsms array, specifies the groups/representation
            used to group or weigth observations
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
    
    # profiles (partition based prototypes)
    if group in adata.obs.columns :

        # working on the layer 
        array = adata.layers[layer]
        membership_vector = adata.obs[group]

        # horizontal split of the array
        counts_splits = {group_x : array[np.where(membership_vector == group_x)] for group_x in membership_vector.unique()}

        # computing feature-wise averages 
        profiles = {}
        for group_i, counts in counts_splits.items() :
            # dense single layer conversion. 
            profiles[group_i] = np.array(counts.mean(axis = 0)).reshape(-1)
            # profiles[group] = np.mean(counts)

        profiles = pd.DataFrame.from_dict(profiles)
        profiles = profiles.set_index(adata.var_names) 

    # weighted profiles (decomposition based prototypes)
    elif group in adata.obsm : 
        # scaling
        # scaling between 0 and 1 the contribution of each prototype to any given cell
        prototypes = adata.obsm[group] 
        # matrix operation goes speed! 
        scaled_prototypes = prototypes.T @ np.diag(1 / prototypes.sum(axis = 1))

        # weight 
        expression = adata.layers[layer]
        weighted = csr_matrix(scaled_prototypes) @ expression # both sparse 

        # mean 
        # ugly rearranging due to 
        denominator = csr_matrix(np.diag( 1 / scaled_prototypes.sum(axis = 1)))
        weighted_profiles = denominator @ weighted

        # proper formatting for Pandas DataFrame
        profiles = pd.DataFrame.sparse.from_spmatrix(weighted_profiles.T)
        profiles.columns = [str(profile) for profile in profiles.columns]
        profiles.index = adata.var_names

    else : 
        raise(ValueError('please specify a field in .obs or .obsm'))

    # directly return the pandas Data Frame 
    if return_panda : 
        return(profiles)

    # storing in the anndata  

    # deprecated
    # description: Results are stored in .var; one column for each representation/group.  
    #adata = _panda_in_var_(adata, profiles, field_name = field_name, group_names = group_names)
            
    adata.varm[field_name] = profiles

    return(adata)