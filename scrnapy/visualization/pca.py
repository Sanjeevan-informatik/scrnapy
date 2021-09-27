import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import numbers
import numpy as np
import pandas as pd
import warnings




def pca_graph():

	pc_op = PCA()


	data = pd.read_csv(r"C:\Users\maky\Desktop\SRNAPY\data\test_10X\output\data.csv")  

	index_data= data["0"]
	del data['0']

	data = data.set_index([pd.Index(index_data)])

	data_pcs = pc_op.fit_transform(data)
	plt.plot(data_pcs[:,0], data_pcs[:,1], 'o', color='black')

	# Saving the figure.
	plt.savefig(r"C:\Users\maky\Desktop\SRNAPY\data\test_10X\output\output.jpg")

	return "the PCA graph is successfully created"



resutlt = pca_graph()

print(resutlt)