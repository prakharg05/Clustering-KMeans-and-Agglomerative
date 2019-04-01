currentClusters=[] # cluster of the most upper level
timeStamp=0
class Node: 
# a tree with bottom up nodes used to save different cluster. When two nodes merge together both left and right link in the node points towards
# nodes merged to from th higher level node
	
	# cluster=[]


	def __init__(self,data=None,left=None,right=None):
		global timeStamp
		self.cluster=[]
		if(left==None and right==None):
			# print(data)
			self.cluster.append(data)
			self.timestamp=timeStamp

		else:
			self.left=left
			self.right=right
			self.combineCluster()
			timeStamp+=1
			self.timestamp=timeStamp
		



	def combineCluster(self):
		self.cluster.extend(self.left.cluster)
		self.cluster.extend(self.right.cluster)

	def nodeData(self):
		return self.cluster

	def getTimeStamp(self):
		return self.timestamp

	def getChildren(self):
		return self.left,self.right




from Bio import SeqIO # required libraries imported
from Bio import pairwise2 as pw
import numpy as np
import copy 

data={} # dictionary used to save data 
nodeData={} # dict of objects(node) refering to sequence number
import pickle
scores={} # dictionary used to save scores i.e. similarity matrix

# Function which finds the single linkage distance

def singledist(clust1,clust2):
	m=[]
	m.append(clust1.cluster[0])
	m.append(clust2.cluster[0])
	dist_min=-1*scores[tuple(sorted(m))]

	for i in clust1.cluster:
		
		for j in clust2.cluster:
			l=[]
			l.append(i)
			l.append(j)
			distance=-1*scores[tuple(sorted(l))]
			
			if(distance<dist_min):
				dist_min=distance

	return dist_min

# Function which finds the complete linkage distance

def completedist(clust1,clust2):
	m=[]
	m.append(clust1.cluster[0])
	m.append(clust2.cluster[0])
	dist_max=-1*scores[tuple(sorted(m))]

	for i in clust1.cluster:
		
		for j in clust2.cluster:
			l=[]
			l.append(i)
			l.append(j)
			distance=-1*scores[tuple(sorted(l))]
			
			if(distance>dist_max):
				dist_max=distance

	return dist_max

# Function which finds the average linkage distance

def averagedist(clust1,clust2):
	totaldist=0
	count=0
	l1=len(clust1.cluster)
	l2=len(clust2.cluster)
	for i in clust1.cluster:
		for j in clust2.cluster:
			l=[]
			l.append(i)
			l.append(j)
			distance=-1*scores[tuple(sorted(l))]
			totaldist+=distance
			count+=1
	return totaldist/(l1*l2)

# function which finds two clusters with minimum average linkage distance

def findAverageMinDistCluster():

	clust1=currentClusters[1]
	clust2=currentClusters[0]
	minDist=averagedist(currentClusters[0],currentClusters[1])# runtime edgecase ,node objects get passed
	for a in currentClusters:
		for b in currentClusters:

			if a!=b:
				distance=averagedist(a,b)# node objects get passed
				if(distance<minDist):
					minDist=distance
					clust1=a
					clust2=b

	return clust1,clust2

# function which finds two clusters with minimum single linkage distance

def findSingleMinDistCluster():

	clust1=currentClusters[1]
	clust2=currentClusters[0]
	minDist=singledist(currentClusters[0],currentClusters[1])# runtime edgecase ,node objects get passed
	for a in currentClusters:
		for b in currentClusters:

			if a!=b:
				distance=singledist(a,b)# node objects get passed
				if(distance<minDist):
					minDist=distance
					clust1=a
					clust2=b

	return clust1,clust2

# function which finds two clusters with minimum complete linkage distance


def findCompleteMinDistCluster():

	clust1=currentClusters[1]
	clust2=currentClusters[0]
	minDist=completedist(currentClusters[0],currentClusters[1])# runtime edgecase ,node objects get passed
	for a in currentClusters:
		for b in currentClusters:

			if a!=b:
				distance=completedist(a,b)# node objects get passed
				if(distance<minDist):
					minDist=distance
					clust1=a
					clust2=b

	return clust1,clust2

#function to get n different clusters

def getClusters():
	clusters=[]
	clusters=list((currentClusters))
	return clusters




#main function which calls al other functions to achieve agglomerative clustering

def main():
	c=0
	l=[]
	for record in SeqIO.parse("vertebrate.txt","fasta"):
		l=[]
		data[c]=("".join((list(record))))
		nodeData[c]=Node(data=c)
		# print(nodeData[c].cluster)
		currentClusters.append(nodeData[c])
		c+=1
	global scores

	scores_r=pickle.load(open("./scores.pkl","rb")) #a pickle file is used to save similarity matrix between different sequences
	scores=copy.deepcopy(scores_r)


	opt = int(input("Enter 1 for single link \n 2 for complete link \n 3 for average link\n"))
	opt2 = int(input("Enter number of clusters needed"))

	for i in currentClusters:
		# pass
		print(i.cluster)
	if(opt==1):
		while(len(currentClusters)!=1):

			clust1,clust2=findSingleMinDistCluster()
			currentClusters.remove(clust1)      # to find clusters using single linkage
			currentClusters.remove(clust2)
			clust=Node(left=clust1,right=clust2)
			currentClusters.append(clust)
			print("Merging :",sorted(clust1.cluster),"   ",sorted(clust2.cluster))
			if(len(currentClusters)==opt2):
				final = getClusters()
		root=currentClusters[0]
	elif(opt==2):		# to find clusters using complete linkage
		while(len(currentClusters)!=1):

			clust1,clust2=findCompleteMinDistCluster()
			currentClusters.remove(clust1)
			currentClusters.remove(clust2)
			clust=Node(left=clust1,right=clust2)
			currentClusters.append(clust)
			print("Merging :",sorted(clust1.cluster),"   ",sorted(clust2.cluster))
			if(len(currentClusters)==opt2):
				final = getClusters()

		root=currentClusters[0]
	elif(opt==3):			# to find clusters using average linkage
		while(len(currentClusters)!=1):

			clust1,clust2=findSingleMinDistCluster()
			currentClusters.remove(clust1)
			currentClusters.remove(clust2)
			clust=Node(left=clust1,right=clust2)
			currentClusters.append(clust)
			print("Merging :",sorted(clust1.cluster),"   ",sorted(clust2.cluster))
			if(len(currentClusters)==opt2):
				final = getClusters()

		root=currentClusters[0]


	print()
	print(opt2," clusters are:")

	for clust in final:
		print(sorted(clust.cluster))
		print()

if __name__=="__main__":
	main()


