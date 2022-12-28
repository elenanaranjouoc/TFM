import rpy2.robjects as ro
from rpy2.robjects import  pandas2ri, r

def create_expressionset_filter(expres, metadata):
    """
    This function is subdivided into two R functions, which we will access through the rpy2 library:
    - create_ExpSet: create an expressionSet
    - filter_ExpSet: reduces the size of expressionset using the nsfilter function

    Parameters:
    expres: Expression matrix
    metadata: Phenotypic variables

    Returns:
    df_expres_filt: Dataframe containing expression matrix after being filtered by the nsfilter function
    expresionSet: R object of class expressionSet 
    """

    pandas2ri.activate()
    expres_r = ro.conversion.py2rpy(expres)
    metadata_r = ro.conversion.py2rpy(metadata)

    ro.globalenv["expres_r"] = expres_r
    ro.globalenv["metadata_r"] = metadata_r
    
    r('library(dplyr)')
    r('library(Biobase)')
    r('library("genefilter")')
    r('library(hgu133plus2.db)')


    r("""
    create_ExpSet <- function(expres, metadata) {
        exprs_matrix <- as.matrix(expres[,-1])
        rownames(exprs_matrix) = expres[,1]
        metadatos = data.frame(labelDescription = c( "title","geo_accession","status", "submission_date","last_update_date", "type","channel_count","source_name_ch1","organism_ch1","taxid_ch1","ER","HER2","B-R grade","node status","LVI","tumor type","tumor size","treatment_protocol_ch1","molecule_ch1","extract_protocol_ch1","label_ch1","label_protocol_ch1","hyb_protocol","scan_protocol","description","data_processing","platform_id","contact_name","contact_email","contact_laboratory","contact_department","contact_institute","contact_address","contact_city","contact_state","contact_zip/postal_code","contact_country","supplementary_file","relation","series_id","data_row_count"),row.names=colnames(metadata))
        fenotipo <- AnnotatedDataFrame(data=metadata, varMetadata=metadatos)
        datosexperimento = new("MIAME",name="James Dirk Iglehart",lab="Iglehart lab",contact ="JIGLEHART@PARTNERS.ORG",title = "Predicting Features of Breast Cancer with Gene Expression Patterns")
        expr_Set <- ExpressionSet(assayData=exprs_matrix, phenoData=fenotipo, experimentData = datosexperimento, annotation = "hgu133plus2")
        return(expr_Set)
    }

    filter_ExpSet <- function(expr_Set) {

        filt1 = nsFilter(expr_Set,var.func=IQR,var.cutoff=0.5,require.GOBP=TRUE)
        filt2 = nsFilter(expr_Set,var.func=sd,var.cutoff=0.5,require.GOBP=TRUE)
        sel = intersect(featureNames(filt1), featureNames(filt2))
        expr_Set_filt <- expr_Set[sel,]

    return(expr_Set_filt)
    }

    expSet <- create_ExpSet(expres_r, metadata_r)
    expSet_filt <- filter_ExpSet(expSet)
    
    mtx_exprs = t(exprs(expSet_filt))
    df_expres_r <- data.frame(mtx_exprs)

    """)
    df_expres_r = ro.globalenv['df_expres_r']

    expresionSet = ro.globalenv['expSet_filt']

    df_expres_filt = ro.conversion.rpy2py(df_expres_r)

    return df_expres_filt, expresionSet
