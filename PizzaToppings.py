#!python3
#
# Pizza Toppings, from
#
# AMS :: Feature Column from the AMS
# http://www.ams.org/publicoutreach/feature-column/fcarc-voting-decision
#
# Math Alive -- Voting & Social Choice -- Lab 1
# http://web.math.princeton.edu/math_alive/6/Lab1.shtml

from PrefVote import *

# AMS version, with
xa = "Sausage"
xb = "Anchovies"
xc = "Peppers"
xd = "Artichoke"
xe = "Mushrooms"

# Math Alive version:
# A = Molson
# B = Killians
# C = Guinness
# D = Samuel Adams
# E = Meister Brau

Ballots = (
	(18,(xa,xd,xe,xc,xb)),
	(12,(xb,xe,xd,xc,xa)),
	(10,(xc,xb,xe,xd,xa)),
	(9,(xd,xc,xe,xb,xa)),
	(4,(xe,xb,xd,xc,xa)),
	(2,(xe,xc,xd,xb,xa))
)

# Prefix, infix, suffix for making formatting for something else
# Here, BBCode
pfx = '[tr][td]'
ifx = '[/td][td]'
sfx = '[/td][/tr]'
def DumpBallots(Blts):
	for Blt in Blts:
		num = Blt[0]
		vts = Blt[1]
		print(pfx + ifx.join([str(num)] + list(vts)) + sfx)

BBox = BallotBox(Ballots)

print("Candidates:",)
Cands = BBox.Candidates()
NumCands = len(Cands)
for Cand in Cands: print(Cand, " ",)
print()
print()

print("Top One (First Past the Post, Plurality):")
res = TopOne(BBox)
for r in res: print(r)
print()

print("Top-Two Runoff:")
reslist = TopTwoRunoff(BBox)
for k, res in enumerate(reslist):
	print("Round", k+1)
	for r in res: print(r)
print		

print("Sequential Runoff:")
reslist = SequentialRunoff(BBox)
for k, res in enumerate(reslist):
	print("Round", k+1)
	for r in res: print(r)
print()

print("Borda Count:")
res = Borda(BBox)
for r in res: print(r)
print()

print("Condorcet Matrix:")
PrefMat = BBox.CondorcetMatrix()
for PMRow in PrefMat:
	for PMVal in PMRow:
		print(PMVal, '\t', end='')
	print()
print()

print("Condorcet Winner:")
print(CondorcetWinner(BBox))
print()