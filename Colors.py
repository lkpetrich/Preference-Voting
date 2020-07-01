#!/usr/bin/env python
#
# Example from
# Voting Systems and the Condorcet Paradox | Infinite Series - YouTube
# https://www.youtube.com/watch?v=HoAnYQZrNrQ
# from
# Feature Column from the AMS
# http://www.ams.org/samplings/feature-column/fcarc-voting-decision

from PrefVote import *

# Similar to the beer ballots in PrefVote.py
ColorBallotSource = (
	("Green","Blue","Purple","Red","Orange"),
	(
		(18, (1,5,4,2,3)), (12, (5,1,4,3,2)), (10, (5,2,1,4,3)),
		(9, (5,4,2,1,3)), (4, (5,2,4,3,1)), (2,(5,4,2,3,1))
	)
)
ColorBallots = CWLPrefNumsToOrder(ColorBallotSource)

DumpAll(ColorBallots)
