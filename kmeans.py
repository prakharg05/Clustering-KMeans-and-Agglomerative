from Bio import SeqIO
from Bio import pairwise2 as pw  # required libraries imported
import numpy as np 
import copy
import pickle
data={} # dictionary used to save data 
scores={} # dictionary used to save scores i.e. similarity matrix

def proximity(k1,k2): #this returns the similarity between two sequences pointed by k1 and k2. k1 and k2 are used as key to point to a prticular sequence

	l=[]
	l.append(k1)
	l.append(k2)
	score=scores[tuple(sorted(l))]
	return score
# Initialize function initializes centroids for the first run. Initially sample data is used in order to remove outliers.
# Then a random point is taken to be first centroid. Then for each succesive initial centroid a point which is farthest from it is selected.
# In this way we obtain a set of initial centroids.

def initialize(k): 
	randomSample = np.random.choice(len(data),31,replace=False)
	l=list(randomSample)

	randomPoint = l[0]
	l.remove(randomPoint)

	centroids=[]
	centroids.append(randomPoint)
	centroid =-1
	count=1
	minval=99999
	while(count<k):
		count+=1
		minval=99999
		for point in centroids:
			for secondpoint in l:
				simmilarity=proximity(point,secondpoint)
				if simmilarity < minval:
					minval=simmilarity
					centroid=secondpoint
		centroids.append(centroid)
		l.remove(centroid)
	
	return centroids

# This function takes clusters and finds its centroid. 
def centroidUpdate(cluster):# CLUSTER is a list

	centroid=[]
	maxval=0
	for candidate in cluster:
		val=0
		for element in cluster:
			if(candidate!=element):
				val+=proximity(candidate,element)
		
		if(val>maxval):
			maxval=val
			centroid=candidate

	return centroid
# Given a centroid this function determines all the points which would be in the cluster represented by this centroid

def assignPoints(centroids): #centroid is list of keys
	
	cluster={}
	for centroid in centroids:
		temp=[]
		temp.append(int(centroid))

		cluster[centroid]=tuple(temp)

	t_key=cluster.keys()
	clust=0
	for key in data:
		maxval=-99999
		for centroid in centroids:

			if(proximity(centroid,key)>maxval):
				maxval=proximity(centroid,key)
				clust=centroid


		l=(list(cluster[clust]))
		l.append(key)
		cluster[clust]=tuple(set(tuple(l)))

	return cluster

# This functions first calls initialize function to initialize centroids
# then it uses termination control function in order to determine when to stop
# while the termination control function returns true while loops run and calls assignPoints function and centroidUpdate function



def kMeans(k):

	centroids=initialize(k)
	c=0
	clusters={}
	previousCluster={}
	while(terminationCondition(previousCluster,clusters)):

		previousCluster=copy.deepcopy(clusters)
		clusters=assignPoints(centroids)

		newCentroids=[]	
		for cluster in clusters:
			newCentroids.append(centroidUpdate(clusters[cluster]))
		centroids=copy.deepcopy(newCentroids)

		c+=1

	return clusters
#This function uses hard condition and checks if any of the points change their clusters or not and determines when to stop.

def terminationCondition(prev,current):
	prev_key=sorted(list(prev.keys()))
	current_key=sorted(list(current.keys()))
	if(prev_key==[]):
		return True
		
	counter=0

	if(prev_key != current_key):
		return True

	for key in prev_key:
		for element in prev[key]:
			if( element in current[key]):
				pass
			else:
				counter+=1
				
	if(counter>0):
		return True

	return False


def main():
	c=0

	for record in SeqIO.parse("vertebrate.txt","fasta"):
		data[c]=("".join((list(record))))
		c+=1
	global scores

	scores_r=pickle.load(open("./scores.pkl","rb")) #a pickle file is used to save similarity matrix between different sequences
	scores=copy.deepcopy(scores_r)

	opt2 = int(input("Enter number of clusters needed"))

	cluster=kMeans(opt2)
	for key in cluster:
		print(key," : ",cluster[key])




if __name__ == '__main__':
	main()


print("\nEND\n")