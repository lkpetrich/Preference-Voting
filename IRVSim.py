#!/usr/bin/env python
#
# For simulating Instant Runoff Voting
#
# Source: http://talkfreethought.org/showthread.php?9749-How-will-Instant-Runoff-Voting-in-Maine-effect-voting-reform
# Blahface #11
#
# Voters are an even spectrum between 0 and 100 inclusive

from PrefVote import *

candsets = ((26,40), (26,40,77), (16,26,40,77), (16,26,40,46,77))

def sortcands(cands,voter):
	cds = list(cands)
	cds.sort(lambda a,b: cmp(abs(a-voter),abs(b-voter)))
	return tuple(cds)

def makeballots(cands):
	ballots = [sortcands(cands,voter) for voter in xrange(0,101,1)]
	return tuple(ballots)

def makebbox(cands):
	return BallotBox(MakeWeighted(makeballots(cands)))

print
for cands in candsets:
	print "Candidates:", cands
	BBox = makebbox(cands)
	print "Top One"
	print TopOne(BBox)
	print "Top-Two Runoff"
	for round in TopTwoRunoff(BBox): print round
	print "Sequential Runoff"
	for round in SequentialRunoff(BBox): print round
	print "Schulze"
	print Schulze(BBox)
	print "Borda"
	print Borda(BBox)
	
	print
