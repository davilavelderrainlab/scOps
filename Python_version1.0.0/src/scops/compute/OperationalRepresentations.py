# package specific function
from scops.compute.BOGs import BOGs
from scops.compute.Profiles import Profiles
from scops.compute.Signatures import Signatures

def OperationalRepresentations(adata, 
    group,
    layer='logcounts',
    verbose=True,
    profile_field='Profiles',
    signature_field='Signatures',
    **kwargs) :
    """Computes operational representations for the chosen group representation

    The function combines computeProfiles, computeSignatures and computeBOGs 
    using a membership vector from .obs or an array from .obsm. 

    It returns Profiles, Signatures and BOGs in the anndata. 
    
    Parameters
    ----------
    adata       
        AnnData object: stores your data 
    group 
        str: .obs column or .obsms array, specifies the groups/representation
            used to obtain operational representations.
    layer
        str='logcounts': layer of adata.layer. specifies the counts used
            to compute operational representations. 
    verbose
        bool=True: wheter to have info from the function while it is running 
    profile_field
        str='Profile_': string used to store profiles in adata.varm
    signature_field
        str='Signature_': string used to store signatures in adata.varm
    **kwargs
        other named arguments to be passed to computeBOGs

    Returns
    -------
    AnnData  
        with new .uns fields for BOGs and Profiles, Signatures and BOGs in .var

    Examples
    --------
    adata = anndata.read_h5ad(path_to_anndata)
    adata = computeOperationalRepresentations(adata, group = 'leiden')
    
    Notes
    -----
    Performance warnings will arise if BOGs are computed on many groups, due to 
    scanpy functions. To ignore them, set warnignsOff = True.
    """

    # Profiles 
    if verbose : 
        print(f'--- computing profiles for {group} using {layer} ---')
    adata = Profiles(adata, group = group, layer = layer, field_name = profile_field)

    # Signatures 
    if verbose : 
        print('--- computing signatures from profiles ---')
    adata = Signatures(adata, profile_field = profile_field, field_name = signature_field)

    # Bags Of Genes (BOGs)
    if verbose : 
        print(f'--- computing BOGs for {group} using {layer} ---')
    adata = BOGs(adata, layer = layer, group = group, **kwargs)

    return(adata)
