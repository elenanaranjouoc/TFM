import pandas as pd
import numpy as np

from get_and_save_data import get_data, extract_data
from expresion import create_expressionset_filter
from reduccion_dimension import pca_f, tsne, plots
from grupos_genes import grupo_genes, gsa

"""
    OBTENCIÓN DE DATOS
"""

# Extraemos datos
extract_data('GSE5460', './data/')

# Guardamos datos
expres, metadata = get_data('./data/')

# Categorizamos tumor size
metadata['tumor size'] = metadata['tumor size'].astype('float64')
datos_unicos = metadata['tumor size'].unique()

conditionlist = [
    (metadata['tumor size']<1),
    (metadata['tumor size']>=1) & (metadata['tumor size']<2),
    (metadata['tumor size']>=2) & (metadata['tumor size']<3),
    (metadata['tumor size']>=3) & (metadata['tumor size']<4),
    (metadata['tumor size']>=4) & (metadata['tumor size']<5),
    (metadata['tumor size']>=5) & (metadata['tumor size']<6),
    (metadata['tumor size']>=6) & (metadata['tumor size']<7),
    (metadata['tumor size']>=7) & (metadata['tumor size']<8),
    (metadata['tumor size']>=8)
]

choicelist = ['0-1', '1-2', '2-3', '3-4', '4-5', '5-6', '6-7', '7-8', '8-9']
metadata['tumor_size_cat'] = np.select(conditionlist, choicelist, default='Not Specified')

# Filtramos datos
df_expres_filt, expresionSet = create_expressionset_filter(expres, metadata)

# Procedemos con PCA
pca, df_pca_trans = pca_f(df_expres_filt, 129)


# Volvemos a ejecutar con 2 componentes principales
pca_2, df_pca_trans_2 = pca_f(df_expres_filt, 2)

#Ploteando
plots(df_pca_trans_2, metadata, 'PCA', 'ER')
plots(df_pca_trans_2, metadata, 'PCA', 'HER2')
plots(df_pca_trans_2, metadata, 'PCA', 'B-R grade')
plots(df_pca_trans_2, metadata, 'PCA', 'node status')
plots(df_pca_trans_2, metadata, 'PCA', 'tumor_size_cat')
plots(df_pca_trans_2, metadata, 'PCA', 'tumor type')
plots(df_pca_trans_2, metadata, 'PCA', 'LVI')

# TSNE
tsne, df_tsne = tsne(df_expres_filt, 2)

#Ploteando
plots(df_tsne, metadata, 'TSNE', 'ER')
plots(df_tsne, metadata, 'TSNE', 'HER2')
plots(df_tsne, metadata, 'TSNE', 'B-R grade')
plots(df_tsne, metadata, 'TSNE', 'node status')
plots(df_tsne, metadata, 'TSNE', 'tumor type')
plots(df_tsne, metadata, 'TSNE', 'tumor_size_cat')
plots(df_tsne, metadata, 'TSNE', 'LVI')

# Grupo de genes
grupo_genes(expresionSet)

# Grupo de genes: GSA
gsa(expresionSet)
