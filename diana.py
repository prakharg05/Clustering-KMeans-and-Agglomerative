import copy
from Bio import SeqIO
from Bio import pairwise2 as pw
import numpy as np
import pickle
currentClusters=[]
scores={}
data={}

class Node: # A tree data structure used to save different clusters
# WHen a cluster splits into two parent link in node is used to point towards it parent


	def __init__(self,data,parent=None):
		self.cluster=data
		self.parent=parent

# A function which finds average distance between a point and a given clusters. It ignores itself if present in cluster


def averageDistance(clustList,start):
	distance=0
	c=0
	if(len(clustList)==1):
		return -99999
	for j in clustList:
		l=[]
		c+=1
		if start!=j:
			l.append(start)
			l.append(j)
			distance=distance+(-1*scores[tuple(sorted(l))])

	return distance/c

# this function uses averageDistance function in order to find the point which would act as splinter group

def maxDistantPoint(clust):
	distMax=-99999
	val=-1
	for i in clust.cluster:
		distance=averageDistance(clust.cluster,i)
		if(distance>distMax):
			distance=distMax
			val=i
	return val

# def mindist(l,num):
# 	mindistance=999999
# 	for i in l:
# 		distance=averageDistance(l,num)
# 		if(mindistance>distance):
# 			mindistance=distance
# 	return i,mindistance

# this function takes a cluster then finds the splinter point using max distant point
# then using average distance function it calculates average distance of a point between in a non splinter group to other points in the same group
# If this distance is greater than the points of splinter group
# Then this functions adds this point to the splinter group


def split(clust):
	l=copy.deepcopy(clust.cluster)
	left=[]
	right=[]
	temp=[]
	num=maxDistantPoint(clust)
	left.append(num)
	l.remove(num)
	for i in l:
		interDistance=averageDistance(left,i)
		intraDistance=averageDistance(l,i)
		if(intraDistance>interDistance):
			left.append(i)
			l.remove(i)

	return left,l


# this function finds the cluster with maximum diameter i.e. the cluster which has to be split into two clusters

	
def maxDiameter(currentClusters):
	maxD=-99999
	clust=currentClusters[0]
	for node in currentClusters:
		for i in node.cluster:
			for j in node.cluster:
				l=[]
				if i!=j:
					l.append(i)
					l.append(j)
					distance=-1*scores[tuple(sorted(l))]
					if(distance>maxD):
						maxD=distance
						clust=node

	return clust

#this functions return a list containing n clusters

def getClusters():
	clusters=[]
	clusters=list((currentClusters))
	return clusters

#Main functions calls all other function in order to achieve divisive analysis 

def main():
	c=0
	l=[]
	
	for record in SeqIO.parse("vertebrate.txt","fasta"):
		data[c]=("".join((list(record))))
		l.append(c)
		c+=1


	global scores
	scores_r=pickle.load(open("./scores.pkl","rb"))
	scores=copy.deepcopy(scores_r)
	opt2 = int(input("Enter number of clusters needed"))

	root=Node(l,None)
	currentClusters.append(root)

	while(len(currentClusters)!=len(data)):
		node = maxDiameter(currentClusters)
		left,right=split(node)
		currentClusters.remove(node)
		leftNode=Node(left,node)
		rightNode=Node(right,node)
		print("Spliting "," ",sorted(node.cluster) , "--->", sorted(left),"   ",sorted(right),"  ",len(currentClusters))
		if(len(currentClusters)==opt2):
			final = getClusters()
		currentClusters.append(leftNode)
		currentClusters.append(rightNode)


	print()
	print(opt2,"clusters are : ")
	print()
	for clust in final:
		print(sorted(clust.cluster))
		print()





if __name__ == '__main__':
	main()


