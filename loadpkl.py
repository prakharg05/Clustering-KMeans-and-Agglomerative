import pickle

scores=pickle.load(open("./scores.pkl","rb"))

c=0
for i in scores:
	print(i,scores[i])
	c+=1

	if c==35:
		break