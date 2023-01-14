from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def pca_f(df, n_compon):
    """
    Applies the pca algorithm to a dataframe that will be previously standardized. Saves the variance explained by the algorithm in a csv file.

    Parameters:
    df: the dataframe on which the pca is to be applied
    n_compon: number of pca components

    Returns:
    pca: trained pca model
    df_PCA: dataframe transformed using the pca
    """
    #normalizar datos
    scaler = StandardScaler()
    df_scaled = scaler.fit_transform(df)

    #comienzo pca
    pca = PCA(n_components = n_compon)
    X_trans = pca.fit_transform(df_scaled)
    df_PCA = pd.DataFrame(data = X_trans, columns = ["PC"+str(i+1) for i in range(n_compon)])
    df_PCA = df_PCA.set_index(df.index)
    #Guardamos los datos de la varianza explicada
    cumVar = pd.DataFrame(np.cumsum(pca.explained_variance_ratio_)*100, 
                      columns=["cumVarPerc"])
    expVar = pd.DataFrame(pca.explained_variance_ratio_*100, columns=["VarPerc"])
    df_var = pd.concat([expVar, cumVar], axis=1).rename(index={i:f"PC{i+1}" for i in range(pca.n_components_)})
    df_var.to_csv('./results/pca/Varianza explicada de las componentes tras pca con '+str(n_compon)+'.csv')
    return pca, df_PCA


def plots(df, metadata, method, colu):
    """
    Saves in the results folder the images of the dataframe transformed by the method used for variable reduction and compared with the selected phenotypic data column.

    Parameters:
    df: Dataframe transformed using the selected method
    metadata: Phenotypic variables
    metodo: Method used for variable reduction
    colu: Selected phenotypic variable

    """
    Y = pd.DataFrame(metadata[colu], columns=[colu])

    df.index= Y.index
    df_a = pd.concat([df, Y], axis=1)

    plt.figure()
    plt.figure(figsize=(12, 6))
    sns.scatterplot(data=df_a, x="PC1", y="PC2", hue=colu)
    plt.savefig("./results/"+method+"/"+colu+".jpg")




