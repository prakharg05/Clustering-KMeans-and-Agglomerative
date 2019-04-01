from Bio import SeqIO
from Bio import pairwise2 as pw
from google.colab import files

scores={}
data={}

def proximity(k1,k2):

	t_key=k1,k2
	key=tuple(sorted(t_key))

	if(key in scores):
		# print("breast")
		pass
	
	else:
		scores[key]=pw.align.globalms(data[k1],data[k2],2,-1,-0.5,-0.3)[0][2]


def main():
	c=0

	for record in SeqIO.parse("vertebrate.txt","fasta"):
		data[c]=("".join((list(record))))
		c+=1

	for i in data:
		for j in data:
			print(i,j)
			proximity(i,j)

	pickle.dump(scores,open("./scores.pkl","wb"))
	files.download("./scores.pkl")

if __name__=="__main__":
	main()


