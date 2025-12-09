from scops.helper._panda_in_var_ import _panda_in_var_
from scops.helper.GetRepresentation import GetRepresentation


# signatures 
# from profiles 
def Signatures(adata,
    profile_field="Profiles",
    field_name='Signatures',
    return_panda=False) :
     """Computes Signatures from Profiles

    A signature is the summed pair-wise differences of a profile vs all the 
    others. 

    Signatures are computed from profiles (which can be computed 
    via computeProfiles) and stored in adata.varm 

    Parameters
    ----------
    adata       
        AnnData object: stores your data 
    profile_field 
        str"Profiles": field of adata.varm used to store profiles
    field_name
        str='Signatures': name used to store the signature matrix in adata.varm
    return_panda
        bool=False: wheter to return the pandas DataFrame of signatures instead 
            of the AnnData object

    Returns
    -------
    AnnData  
        with new Signatures in .var
        OR a pandas DataFrame features x signatures if return_pandas=True.

    Examples
    --------
    adata = anndata.read_h5ad(path_to_anndata)
    adata = computeProfiles(adata, group = 'leiden')
    adata = computeSignatures(adata)
    
    Notes
    -----
    none
    """
    
    # selecting profile
    profiles = adata.varm[profile_field].copy()

    if profiles.empty : 
        raise(ValueError('please provide a valid profile_field or compute profiles first, using computeProfiles'))

    # faster way to compute signatures 
    # n = len(profiles.columns)
    n_mean = profiles.mean(axis = 1)
    signatures = profiles.sub(n_mean, axis = 0)

    adata.varm[field_name] = signatures
    adata.varm[field_name].columns = adata.varm[field_name].columns.astype('str') 

    if return_panda : 
        return(signatures)

    return(adata)