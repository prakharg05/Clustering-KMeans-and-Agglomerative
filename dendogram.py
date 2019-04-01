import numpy as np
import plotly.plotly as py
import plotly.figure_factory as ff
from Bio import SeqIO
import pickle
import copy
import plotly.io as pio
from scipy.cluster.hierarchy import linkage
scores={}
data={}


def group():
	c=0
	name=[]
	for record in SeqIO.parse("vertebrate.txt","fasta"):
		name.append(record.id)
		data[c]=("".join((list(record))))
		c+=1


	global scores
	scores_r=pickle.load(open("./scores.pkl","rb"))
	scores=copy.deepcopy(scores_r)

	X=np.zeros((len(data),len(data)))
	for t1,t2 in scores:
		X[t1][t2]=-1*scores[(tuple([t1,t2]))]


	dendro = ff.create_dendrogram(X,labels=name,linkagefun=lambda x: linkage(x, 'average'))
	dendro['layout'].update({'width':800, 'height':500})
	py.plot(dendro, filename='group_dendrogram')
	pio.write_image(dendro, 'dendo.jpg')

def complete():
	c=0
	name=[]
	for record in SeqIO.parse("vertebrate.txt","fasta"):
		name.append(record.id)
		data[c]=("".join((list(record))))
		c+=1


	global scores
	scores_r=pickle.load(open("./scores.pkl","rb"))
	scores=copy.deepcopy(scores_r)

	X=np.zeros((len(data),len(data)))
	for t1,t2 in scores:
		X[t1][t2]=-1*scores[(tuple([t1,t2]))]


	dendro = ff.create_dendrogram(X,labels=name,linkagefun=lambda x: linkage(x, 'complete'))
	dendro['layout'].update({'width':800, 'height':500})
	py.plot(dendro, filename='complete_dendrogram')
	# pio.write_image(dendro, 'dendo.jpg')

def single():
	c=0
	name=[]
	for record in SeqIO.parse("vertebrate.txt","fasta"):
		name.append(record.id)
		data[c]=("".join((list(record))))
		c+=1


	global scores
	scores_r=pickle.load(open("./scores.pkl","rb"))
	scores=copy.deepcopy(scores_r)

	X=np.zeros((len(data),len(data)))
	for t1,t2 in scores:
		X[t1][t2]=-1*scores[(tuple([t1,t2]))]


	dendro = ff.create_dendrogram(X,labels=name,linkagefun=lambda x: linkage(x, 'single'))
	dendro['layout'].update({'width':800, 'height':500})
	py.plot(dendro, filename='single_dendrogram')
	# pio.write_image(dendro, 'dendo.jpg')

single()
complete()
group()