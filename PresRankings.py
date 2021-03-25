#!python3
#
# From
# http://en.wikipedia.org/wiki/Historical_rankings_of_Presidents_of_the_United_States
# Historical rankings of Presidents of the United States - Wikipedia, the free encyclopedia
#
# PrefVote.py would have to be modified to make unchosen candidates neutral
# instead of lower in preference than the rest
#
# Does a modification of the Borda count and also calculates standard deviations
# to get a picture of the range of opinion

from math import sqrt
from functools import cmp_to_key

Data = (
	("George Washington",2,2,3,3,4,4,4,3,2,3,1,4,1,2,2,4,3,2,2,2,1),
	("John Adams",9,10,9,14,10,14,12,14,11,16,13,12,13,13,17,17,12,15,19,14,14),
	("Thomas Jefferson",5,5,4,5,2,3,5,4,4,7,4,5,4,4,7,5,4,5,7,5,5),
	("James Madison",14,12,14,17,9,8,9,10,17,18,15,9,17,15,20,6,14,13,17,12,7),
	("James Monroe",12,18,15,16,15,11,15,13,15,14,16,8,16,21,14,7,13,16,13,18,8),
	("John Quincy Adams",11,13,16,19,17,16,17,18,18,19,20,17,25,16,19,19,20,22,21,23,18),
	("Andrew Jackson",6,6,7,6,13,9,11,8,5,13,6,13,10,14,13,14,9,9,18,15,19),
	("Martin Van Buren",15,17,20,18,21,21,22,21,21,30,23,24,27,40,31,23,27,25,34,27,25),
	("William Henry Harrison",None,None,None,38,26,35,28,35,None,37,None,36,None,39,39,35,None,39,28,42,39),
	("John Tyler",22,25,28,29,34,33,34,34,32,36,34,37,35,31,35,37,37,36,39,37,37),
	("James K. Polk",10,8,12,11,12,13,14,11,9,12,10,11,9,9,12,12,16,19,14,20,12),
	("Zachary Taylor",25,24,27,28,29,34,33,29,29,28,31,34,33,28,29,33,33,33,31,35,30),
	("Millard Fillmore",24,26,29,31,32,32,35,36,31,35,35,38,36,33,37,38,35,37,37,38,38),
	("Franklin Pierce",27,28,31,35,35,36,37,37,33,39,37,39,38,41,40,40,39,40,41,41,40),
	("James Buchanan",26,29,33,36,37,38,39,40,38,41,39,41,40,42,42,42,40,43,43,43,43),
	("Abraham Lincoln",1,1,1,1,3,2,2,1,1,1,2,2,2,1,1,3,2,1,1,1,3),
	("Andrew Johnson",19,23,32,32,38,39,40,39,37,40,36,42,37,24,41,43,36,41,42,40,44),
	("Ulysses S. Grant",28,30,35,30,36,37,38,38,34,33,32,35,29,18,23,26,29,28,22,21,24),
	("Rutherford B. Hayes",13,14,22,22,22,23,24,25,23,26,22,27,24,27,33,31,30,30,32,29,32),
	("James A. Garfield",None,None,None,33,25,30,26,30,None,29,None,33,None,34,28,27,None,31,29,34,28),
	("Chester A. Arthur",17,21,26,24,24,26,27,28,26,32,26,30,26,22,32,25,32,32,35,31,34),
	("Grover Cleveland",8,11,17,13,18,17,19,16,13,17,12,20,12,19,21,20,21,23,23,24,23),
	("Benjamin Harrison",21,20,23,25,31,29,30,31,19,31,27,32,30,30,30,34,34,29,30,32,35),
	("William McKinley",18,15,18,10,19,19,18,17,16,15,14,19,14,17,16,21,17,21,16,19,20),
	("Theodore Roosevelt",7,7,5,4,5,5,3,5,6,4,5,3,5,5,4,2,5,4,4,4,4),
	("William Howard Taft",16,16,19,20,20,20,21,20,22,24,19,21,20,29,24,24,25,20,24,22,22),
	("Woodrow Wilson",4,4,6,7,6,6,6,6,7,6,11,6,11,10,9,8,6,10,11,11,11),
	("Warren G. Harding",29,31,36,37,39,40,41,41,39,38,37,40,39,35,38,41,38,42,40,39,41),
	("Calvin Coolidge",23,27,30,27,30,31,36,33,30,27,25,29,23,26,26,29,28,27,27,28,31),
	("Herbert Hoover",20,19,21,21,27,28,29,24,35,34,29,31,31,36,34,36,26,38,36,36,36),
	("Franklin D. Roosevelt",3,3,2,2,1,1,1,2,3,2,3,1,3,3,3,1,1,3,3,3,2),
	("Harry S. Truman",None,9,8,8,7,7,7,7,8,5,7,7,7,7,5,9,7,6,6,6,9),
	("Dwight D. Eisenhower",None,22,11,9,11,12,8,9,10,9,9,10,8,6,8,10,10,7,5,7,6),
	("John F. Kennedy",None,None,13,15,8,10,10,15,12,8,18,14,15,11,6,11,15,14,8,16,10),
	("Lyndon B. Johnson",None,None,10,12,14,15,13,12,14,10,17,15,18,12,11,16,11,12,10,10,16),
	("Richard Nixon",None,None,34,34,28,25,23,32,36,25,33,26,32,38,27,30,23,34,28,33,29),
	("Gerald Ford",None,None,24,23,23,27,32,27,28,23,28,28,28,25,22,28,24,24,25,25,27),
	("Jimmy Carter",None,None,25,26,33,24,25,19,27,22,30,25,34,32,25,32,18,26,26,26,26),
	("Ronald Reagan",None,None,None,None,16,22,20,26,25,11,8,16,6,8,10,18,8,11,9,9,13),
	("George H. W. Bush",None,None,None,None,None,18,31,22,24,20,21,22,21,20,18,22,22,17,20,17,21),
	("Bill Clinton",None,None,None,None,None,None,16,23,20,21,24,18,22,23,15,13,19,8,15,13,15),
	("George W. Bush",None,None,None,None,None,None,None,None,None,None,None,23,19,37,36,39,31,35,33,30,33),
	("Barack Obama",None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,15,None,18,12,8,17),
	("Donald Trump",None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,44,42)
)

Presidents = tuple([d[0] for d in Data])
NumPress = len(Presidents)

NumRankers = len(Data[0]) - 1
NumRanked = NumRankers*[0]

for d in Data:
	for k, r in enumerate(d[1:]):
		if r != None: NumRanked[k] += 1

RelRanks = []
for d in Data:
	rrline = []
	for k,r in enumerate(d[1:]):
		if r != None:
			rlr = 1 - 2*(r-1.0)/(NumRanked[k]-1.0)
		else:
			rlr = None
		rrline.append(rlr)
	RelRanks.append(tuple(rrline))
RelRanks = tuple(RelRanks)

PresNbrs = NumPress*[0]
PresAvgs = NumPress*[0.0]

for kp in range(NumPress):
	for k in range(NumRankers):
		r = RelRanks[kp][k]
		if r != None:
			PresNbrs[kp] += 1
			PresAvgs[kp] += r

for kp in range(NumPress):
	PresAvgs[kp] /= PresNbrs[kp]

PresStds = NumPress*[0.0]

for kp in range(NumPress):
	for k in range(NumRankers):
		r = RelRanks[kp][k]
		if r != None:
			PresStds[kp] += (r - PresAvgs[kp])**2

for kp in range(NumPress):
	PresStds[kp] = sqrt(PresStds[kp]/PresNbrs[kp])

PresStats = list(zip(Presidents,PresNbrs,PresAvgs,PresStds))
for PS in PresStats: print(PS)

print()

PresStats.sort(key=lambda x: (-x[2],x[0]))
for PS in PresStats: print(PS)

print()

PrefMat = [NumPress*[0] for k in range(NumPress)]
for kr in range(NumRankers):
	for k1 in range(NumPress):
		r1 = Data[k1][kr+1]
		if r1 != None:
			for k2 in range(NumPress):
				r2 = Data[k2][kr+1]
				if r2 != None:
					if r1 < r2:
						PrefMat[k1][k2] += 1

def cmp(x,y):
	if x > y: return 1
	elif x < y: return -1
	else: return 0

def SchulzeOrdering(BPMat,i,j):
	rc = - cmp(BPMat[i][j],BPMat[j][i])
	if rc != 0: return rc
	return cmp(i,j)

def SchulzePrefOrder(PrefMat):
	n = len(PrefMat)
	
	# Beatpaths
	BPMat = [n*[0] for k in range(n)]
	
	# Variant of Floyd-Warshall algorithm
	for i in range(n):
		for j in range(n):
			if j != i:
				if PrefMat[i][j] > PrefMat[j][i]:
					BPMat[i][j] = PrefMat[i][j]
				else:
					BPMat[i][j] = 0
	
	for i in range(n):
		for j in range(n):
			if j != i:
				for k in range(n):
					if i != k and j != k:
						BPMat[j][k] = \
							max(BPMat[j][k], min(BPMat[j][i],BPMat[i][k]))
	
	Ordering = list(range(n))
	Ordering.sort(key=cmp_to_key(lambda i,j: SchulzeOrdering(BPMat,i,j)))
	return Ordering

SPO = SchulzePrefOrder(PrefMat)
PresSP = tuple((Presidents[k] for k in SPO))

PresBC = [p[0] for p in PresStats]

for p in zip(PresBC,PresSP): print(p)