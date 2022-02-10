import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn import preprocessing
# from sqlalchemy import column



def make_PCA(chr):
    print()
    #https://www.youtube.com/watch?v=Lsue2gEM9D0
    # center and scale the data
    # T = samples to be rows instead of columns
    scaled_data = StandardScaler().fit_transform(chr.iloc[1:,1:].T) #preprocessing.scale(chr.iloc[1:,1:].T)
    # create PCA object
    pca = PCA()
    # all the PCA math (i.e. calculate loading scores and the variation each PC accounts for)
    pca.fit(scaled_data)
    # generate coordinates for a PCA graph based on the loading scores and the scaled data
    pca_data = pca.transform(scaled_data)

    #scree plot to see how many PCs should go into the final plot
    # 1) calculate the percentage of variation that each PC accounts for
    per_var = np.round(pca.explained_variance_ratio_*100, decimals=1)
    # 2) create labels for the scree plot. (these are PC1, PC2 etc)
    labels = ['PC' + str(x) for x in range(1, len(per_var)+1)] #['PC' + str(x) for x in range(1, len(per_var)+1)]
    # Create bar plot
    plt.bar(x=range(1,len(per_var)+1), height=per_var) #, tick_label=labels)
    plt.ylabel('Percentage of Explained Variance')
    plt.xlabel('Principal Component (PC)')
    plt.title('Scree Plot')
    plt.show()
    # Create line plot
    plt.plot(per_var) #, tick_label=labels)
    plt.ylabel('Percentage of Explained Variance')
    plt.xlabel('Principal Component (PC)')
    plt.title('Scree Plot')
    plt.show()

    # # Draw a PCA plot
    # # 1) put the new coordinates, created bu pca.transform(scaled.data), into a nice
    # # matrix where the rows have sample labels and the columns have PC labels
    # # TODO vanaf 9.22
    # pca_df = pd.DataFrame(pca_data, index=chr.iloc[1,:].T, columns=labels)
    # # drwa a scatter plot with a title and nice axis labels
    # plt.scatter(pca_df.PC1, pca_df.PC2)
    # plt.title('My PCA graph')
    # plt.xlabel('PC1 - {0}%'.format(per_var[0]))
    # plt.ylabel('PC2 - {0}%'.format(per_var[1]))

    # # Loop adds sample names to the graph
    # for sample in pca_df.index:
    #     plt.annotate(sample, (pca_df.PC1.loc[sample], pca_df.PC2.loc[sample]))

    # plt.show()






    




    #https://towardsdatascience.com/principal-component-analysis-pca-from-scratch-in-python-7f3e2a540c51
    # X = chr.iloc[:,1:].to_numpy()
    # # print(X)
    # X_scaled = StandardScaler().fit_transform(X)
    # print(X_scaled[:5])
    # # features = X_scaled.T
    # cov_matrix = np.cov(X_scaled)
    # print(cov_matrix[:5])
    # values, vectors = np.linalg.eig(cov_matrix)
    # print(values[:5])
    # print(vectors[:5])
    # # explained_variances = []
    # # for i in range(len(values)):
    # #     explained_variances.append(values[i] / np.sum(values))
    
    # # print(np.sum(explained_variances), '\n', explained_variances)
    # #DS2
    # order = values.argsort()[::-1]
    # values = values[order]
    # vectors = vectors[:, order]
    # # 'Scree' plot of values
    # plt.plot(values)
    # plt.show()
    # # Cumulative percentage of values
    # plt.plot(100 * values.cumsum() / sum(values))
    # plt.ylim((0, 110))
    # plt.show()
    # # projections = np.dot(X_scaled, vectors)
    # # # projections = data @ eigenvectors

    # plt.scatter(projections[:,0], projections[:,1], c=projections[:, 2], s=1)
    # plt.xlim(-100, 100)
    # plt.ylim(-100, 100)
    # plt.show()


    



    #DS2
    # Determining the **correlation** matrix
    # data = chr.iloc[:,2:].to_numpy() #.as_matrix()
    # print(data.mean(axis=0))
    # centered = data - data.mean(axis=0)
    # print(centered.std(axis=0))
    # zscores = centered / centered.std(axis=0)
    # correlation_matrix = np.dot(zscores.T, zscores) / (zscores.shape[0] - 1)

    # # Determining principal components
    # eigenvalues, eigenvectors = np.linalg.eig(correlation_matrix)
    # order = eigenvalues.argsort()[::-1]
    # eigenvalues = eigenvalues[order]
    # eigenvectors = eigenvectors[:, order]
    # # 'Scree' plot of eigenvalues
    # plt.plot(eigenvalues)
    # plt.show()
    # # Cumulative percentage of eigenvalues
    # plt.plot(100 * eigenvalues.cumsum() / sum(eigenvalues))
    # plt.ylim((0, 110))
    # plt.show()




def main():
    file_path = 'D:/Hanze_Groningen/STAGE/UMAP/chr17.tsv'
    chr = pd.read_csv(file_path, sep='\t')
    chr = chr.dropna()
    # print(chr.iloc[1:,1:])
    make_PCA(chr)



if __name__ == '__main__':
    main()





