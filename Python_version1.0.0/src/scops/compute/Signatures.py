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
        str"Profile_": field of adata.varm used to store profiles
    field_name
        str='Signature_': name used to store the signature matrix in adata.varm
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
    
    # selecting profiles

    # deprecated
    #profiles = GetRepresentation(adata, profile_field)

    profiles = adata.varm[profile_field]

    if profiles.empty : 
        raise(ValueError('please provide a valid profile_field or compute profiles first, using computeProfiles'))

    # faster way to compute signatures 
    #n = len(profiles.columns)
    n_mean = profiles.mean(axis = 1) #* n #TODO test 
    signatures = profiles.sub(n_mean, axis = 0)
    # splitting by groups 
    signature_groups = [col_name.split(profile_field)[-1] for col_name in signatures.columns]
    signatures.columns = signature_groups

    # deprecated 
    #adata = _panda_in_var_(adata, signatures, field_name = field_name, group_names = group_names)

    adata.varm[field_name] = signatures

    if return_panda : 
        return(signatures)

    return(adata)