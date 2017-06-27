#!/usr/bin/env python
#
# For doing various preference-based voting algorithms
# From individual preferences, it finds overall preferences with various algorithms
#
# MakeWeighted(Ballots) -- unweighted to weighted, all weights 1
#
# BallotBox objects -- init with weighted ballots:
#   list of (weight, (candidates in prefence order))
#
#
# The arg BBox is a ballot-box object
#
# The functions return a sorted list (candidate,count)
# downward by count, then upward by candidate, unless specified otherwise
#
# TotalCount(BBox): returns the sum of the weights
#
# AllVotes(BBox): count of all votes together: adds up all the weights
#
# TopOne(BBox): count the top one of each ballot
#
# MinusBottom(BBox): count the unselected ones
# or the bottom one of each ballot
#
# TopNum(BBox, Num): count the top Num of each ballot
# for simulating approval voting
#
# Borda(BBox): uses the Borda count:
# counting with ((num candidates) + 1 - rank)
#
# ModBorda(BBox): uses a modification of the Borda count:
# counting with ((num votes in ballot) + 1 - rank)
#
# CumulBorda(BBox): like ModBorda, but divides each ballot's counts
# by that ballot's total count, thus simulating cumulative voting
#
# MajorityJudgment(BBox): sorts the candidates
# by the weighted medians of their ranks
#
#
# These runoff methods return a list of results of rounds,
# each of them like the previous functions' output
#
# TopNumList(BBox): returns a list of all TopNum round results
# Bucklin's method cuts off when some candidate gets a majority
#
# TopTwoRunoff(BBox, DoRound=TopOne): repeats with the top two candidates,
# or all the candidates who got the most votes
#
# SequentialRunoff(BBox, ThresholdFunc=min, DoRound=TopOne):
# repeats without the candidates who got the fewest votes,
# continues until all the candidates are used up.
# With its defaults, it is Instant Runoff Voting
# With DoRound = Borda or ModBorda, it's Baldwin
# With that DoRound and also ThresholdFunc = Average, it's Nanson
# With DoRound = MinusBottom, it's Coombs
#
# SingleTransferableVote(BBox, Num): fill N seats (N is second arg)
# in a manner much like sequential runoff / IRV
# Appends a list of each round's winners or losers to that round's results
# + or -, then the list
# Also adds a fake round that lists the winners
#
# CompPairOutcomesSTV(BBox, Num): fill N seats (N is second arg)
# Comparison of Pairs of Outcomes by Single Transferable Vote
# Finds all the possible outcomes and does pairwise comparisons,
# creating a Condorcet matrix for them
# One then uses some Condorcet method on it to find the winners
#
# Schulze STV not implemented because of my difficulty in
# understanding the algorithm
#
#
# Condorcet methods
#
# CondorcetWinner(BBox)
# Returns (True,winner) if there is one, (False,) otherwise
#
# Schulze(BBox)
# Schulze's beatpath method
# Returns simple list of candidates from winners to losers
#
# Copeland(BBox)
# Copeland's pairwise-aggregation method
# Returns sorted list of (candidate, score)
#
# Minimax(BBox, Which): second arg is which sort:
# "wins" -- winning votes
# "marg" -- margins
# "oppo" -- pairwise opposition
# Returns sorted list of (candidate, score)
#
# KemenyYoung(BBox)
# The Kemeny-Young method
# It goes as n! for n candidates,
# so it may need to use simulated annealing or something similar
# for a large number of candidates.
# Returns simple list of candidates from winners to losers
#
# Dodgson(BBox)
# Dodgson's method
# It goes as n! for n candidates
# Finds all permutations of the ballot preferences,
# then finds the Condorcet winner (if any) for each permutation
# Returns a list of (candidate, distance)
# where distance is the smallest permutation distance
# that made that candidate a Condorcet winner.
#
# Permutation distance =
# Minimum number of interchanges from ascending order to that permutation
# I find that value to be (length) - (number of cycles)
#
# RankedPairs(BBox)
# Tideman's ranked-pairs method
# Like Kemeny-Young, but with simple hill-climbing optimization
# Returns simple list of candidates from winners to losers
#
# Maximum Affirmed Majorities and Maximum Majority Voting
# seem identical to it
#
# Maximal lotteries - Wikipedia - https://en.wikipedia.org/wiki/Maximal_lotteries
# [1503.00694] Consistent Probabilistic Social Choice - https://arxiv.org/abs/1503.00694
# http://econweb.ucsd.edu/~jsobel/172aw02/notes9.pdf
# Find (row # beat col #) - (col # beat row #) matrix:
# C = Condorcet matrix
# D = C - transpose(C)
# Find a probability vector p for the candidates:
# p.D > 0
# For solving by linear programming:
# Maximize w with p.D >= w
# Minimize v with p.D <= v
# where sum of p = 1 and p >= 0
#
# MaximalSet(BBox, Type)
# Returns the Condorcet maximal set
# Type = "Smith": each member beats all those outside of it
# Type = "Schwartz": union of all sets of candidates
# where each member beats or ties all those outside of it,
# but which has no proper subset with that property
#
# MaximalSetSequence(BBox, Type)
# Returns a sequence of maximal sets found,
# with the previous ones removed from the ballots

# Turns preference numbering of candidates into an order
# Assumes 1-based numbers
def CandPrefNumsToOrder(Cands, PrefNums):
	Ballot = len(Cands)*[None]
	for c,k in zip(Cands,PrefNums):
		Ballot[k-1] = c
	return tuple(Ballot)

# Does so on an object with (candidates, list of (weight,numbering))
def CWLPrefNumsToOrder(CWL):
	Cands = CWL[0]
	WLs = CWL[1]
	Ballots = [(WL[0], CandPrefNumsToOrder(Cands, WL[1])) for WL in WLs]
	return tuple(Ballots)


# If one has unweighted ballots, make them weighted
def MakeWeighted(Ballots):
	return tuple(((1,b) for b in Ballots))


# Setup for ballots in BallotBox
def BallotSetup(Ballots):
	return tuple([(b[0],tuple(b[1])) for b in Ballots if len(b[1]) > 0])

# These accessors lazy initing of calculated quantities:
# Candidates -- sorted list of candidates
# Condorcet matrix

def Candidates(self):
	if self._Candidates == None:
		CandSet = set()
		for Ballot in self.Ballots:
			CandSet |= set(Ballot[1])
		CandList = list(CandSet)
		CandList.sort()
		self._Candidates = tuple(CandList)
	
	return self._Candidates

def CondorcetMatrix(self):
	if self._CondorcetMatrix == None:
		Cands = self.Candidates()
		NumCands = len(Cands)
		self._CondorcetMatrix = [NumCands*[0] for k in xrange(NumCands)]
		
		CandIndices = {}
		for k, Cand in enumerate(Cands):
			CandIndices[Cand] = k
		
		for Ballot in self.Ballots:
			# Ranks -- ranked candidates get from at least 1 to NumCands
			# (max rank)
			# Unranked ones get 0
			Weight = Ballot[0]
			Votes = Ballot[1]
			Ranks = NumCands*[0]
			for k, Vote in enumerate(Votes):
				Ranks[CandIndices[Vote]] = NumCands - k
			
			# First index: winner
			# Second index: loser
			for k1 in xrange(NumCands):
				for k2 in xrange(NumCands):
					if Ranks[k1] > Ranks[k2]:
						self._CondorcetMatrix[k1][k2] += Weight
	
	return self._CondorcetMatrix

# Container for ballots
# Members:
# Candidates
# Ballots: list of (weight, list of candidates in preference order)
class BallotBox:
	
	# Members:
	def __init__(self, Ballots):
		self.Ballots = BallotSetup(Ballots)
		# Lazy init:
		self._Candidates = None
		self._CondorcetMatrix = None

BallotBox.Candidates = Candidates
BallotBox.CondorcetMatrix = CondorcetMatrix


# Ballot processing
# Function arg: ballot member. Returns whether or not to keep the member
def ProcessBallots(BBox, ProcessFunction):
	NewBallots = []
	for Ballot in BBox.Ballots:
		Weight = Ballot[0]
		Votes = Ballot[1]
		NewVotes = []
		for Vote in Votes:
			if ProcessFunction(Vote):
				NewVotes.append(Vote)
		NewBallots.append((Weight,tuple(NewVotes)))
	
	return BallotBox(NewBallots)

def KeepCandidates(BBox, Cands):
	return ProcessBallots(BBox, lambda bc: bc in Cands)

def RemoveCandidates(BBox, Cands):
	return ProcessBallots(BBox, lambda bc: bc not in Cands)


# Vote-count functions

# In case we want to do fractions:
# Arg: BallotBox object
def TotalCount(BBox):
	Count = 0
	for Ballot in BBox.Ballots:
		Count += Ballot[0]
	return Count


# a and b are (Cand, Count)
# Sort reverse by count, then forward by cand
def CCSortFunction(a,b):
	rc = - cmp(a[1],b[1])
	if rc != 0: return rc
	return cmp(a[0],b[0])

# Generic count function
# Args:
# BallotBox object
# Count function, with args
#   Counts: dict with candidate-name keys
#   Weight: from ballot
#   Votes: list of candidates in preference order)
#   Cands: list of candidates, in case that might be needed by the count
def BBoxCounter(BBox, CountFunction):
	Cands = BBox.Candidates()
	Counts = {}
	for Cand in Cands:
		Counts[Cand] = 0
	
	for Ballot in BBox.Ballots:
		CountFunction(Counts, Ballot[0], Ballot[1], Cands)
	
	CountList = [Counts[Cand] for Cand in Cands]
	CandCounts =  zip(Cands,CountList)
	CandCounts.sort(CCSortFunction)
	return tuple(CandCounts)


# Simple counting

def AllVotesFunction(Counts, Weight, Votes, Cands):
	for Vote in Votes:
		Counts[Vote] += Weight

def AllVotes(BBox):
	return BBoxCounter(BBox, AllVotesFunction)

# http://en.wikipedia.org/wiki/First-past-the-post_voting

def TopOneFunction(Counts, Weight, Votes, Cands):
	if len(Votes) >= 1:
		Counts[Votes[0]] += Weight

def TopOne(BBox):
	return BBoxCounter(BBox, TopOneFunction)

def MinusBottomFunction(Counts, Weight, Votes, Cands):
	Unselected = tuple(set(Cands) - set(Votes))
	if len(Unselected) > 0:
		for Vote in Unselected:
			Counts[Vote] -= Weight
	elif len(Votes) >= 1:
		Counts[Votes[-1]] -= Weight

def MinusBottom(BBox):
	return BBoxCounter(BBox, MinusBottomFunction)

# http://en.wikipedia.org/wiki/Bucklin_voting

def TopNumFunction(Counts, Weight, Votes, Cands, Num):
	NumAdj = min(len(Votes),Num)
	for k in xrange(NumAdj):
		Counts[Votes[k]] += Weight

def TopNum(BBox,Num):
	return BBoxCounter(BBox, lambda Cn,Wt,Vt,Ca: TopNumFunction(Cn,Wt,Vt,Ca,Num))

def TopNumList(BBox):
	Cands = BBox.Candidates()
	NumCands = len(Cands)
	return tuple((TopNum(BBox,k+1) for k in xrange(NumCands)))


# Rank counting

# http://en.wikipedia.org/wiki/Borda_count

def BordaFunction(Counts, Weight, Votes, Cands):
	MaxNum = len(Cands)
	for k,Vote in enumerate(Votes):
		Counts[Vote] += Weight*(MaxNum - k)

def Borda(BBox):
	return BBoxCounter(BBox, BordaFunction)

def ModBordaFunction(Counts, Weight, Votes, Cands):
	MaxNum = len(Votes)
	for k,Vote in enumerate(Votes):
		Counts[Vote] += Weight*(MaxNum - k)

def ModBorda(BBox):
	return BBoxCounter(BBox, ModBordaFunction)

# http://en.wikipedia.org/wiki/Cumulative_voting
def CumulBordaFunction(Counts, Weight, Votes, Cands):
	MaxNum = len(Votes)
	IndCount = {}
	for k,Vote in enumerate(Votes):
		if Vote not in IndCount: IndCount[Vote] = float(0)
		IndCount[Vote] += Weight*(MaxNum - k)
	ICTotal = float(0)
	for k,Vote in enumerate(Votes):
		ICTotal += IndCount[Vote]
	for k,Vote in enumerate(Votes):
		Counts[Vote] += (Weight*IndCount[Vote])/ICTotal

def CumulBorda(BBox):
	return BBoxCounter(BBox, CumulBordaFunction)

# https://en.wikipedia.org/wiki/Majority_judgment
def WeightedMedian(wtrnks):
	# Sort by ranks and find cumulative weights
	wrsort = list(wtrnks)
	wrsort.sort(lambda a,b: cmp(a[0],b[0]))
	
	totwt = 0
	cwtrks = []
	for rk, wt in wrsort:
		cwr = (rk, totwt)
		cwtrks.append(cwr)
		totwt += wt
	
	hftot = 0.5*totwt
	
	# Search through to find contributors to the median:
	# those that make the boolean variable addmdn true. 
	# Add up their contributions in rksum (sum of ranks)
	# and rkcnt (count of ranks).
	rksum = 0.; rkcnt = 0
	n = len(cwtrks)
	for k in xrange(n):
		addmdn = cwtrks[k][1] <= hftot
		if addmdn:
			if k < n-1:
				addmdn = cwtrks[k+1][1] >= hftot
			else:
				addmdn = True
		
		if addmdn:
			rksum += cwtrks[k][0]; rkcnt += 1
	
	# Finally, the median 
	return rksum/rkcnt

def MajorityJudgment(BBox):
	Cands = BBox.Candidates()
	NumCands = len(Cands)
	
	CandRanks = []
	for Ballot in BBox.Ballots:
		wt, prefs = Ballot
		for k,pref in enumerate(prefs):
			CandRanks.append( (pref, NumCands-k, wt) )
		for Cand in Cands:
			if Cand not in prefs:
				CandRanks.append( (Cand, 0, wt) )
	
	CandRankLists = {}
	for Cand in Cands:
		CandRankLists[Cand] = []
	
	for CREntry in CandRanks:
		CandRankLists[CREntry[0]].append(CREntry[1:3])
	
	CandWMs = [(Cand, WeightedMedian(CandRankLists[Cand])) for Cand in Cands]
	
	CandWMs.sort(CCSortFunction)
	return tuple(CandWMs)


# Multistep methods

# http://en.wikipedia.org/wiki/Contingent_vote
# http://en.wikipedia.org/wiki/Two-round_system

def TopTwoRunoff(BBox, DoRound=TopOne):
	rounds = []
	res = DoRound(BBox)
	rounds.append(res)
	
	if len(res) == 0: return tuple(rounds)
	
	# Find top two or top ones with the same vote
	# Using res being sorted
	TopTwoCands = []
	MaxVotes = res[0][1]
	for r in res:
		if r[1] != MaxVotes:
			if len(TopTwoCands) >= 2: break
			MaxVotes = r[1]
		TopTwoCands.append(r[0])
	TopTwoCands.sort()
	TopTwoCands = tuple(TopTwoCands)
	
	NewBBox = KeepCandidates(BBox, TopTwoCands)
	res = DoRound(NewBBox)
	rounds.append(res)

	return tuple(rounds)

# http://en.wikipedia.org/wiki/Instant-runoff_voting
# http://en.wikipedia.org/wiki/Exhaustive_ballot
# http://en.wikipedia.org/wiki/Nanson%27s_method
# http://en.wikipedia.org/wiki/Coombs%27_method

def SequentialRunoff(BBox, ThresholdFunc=min, DoRound=TopOne):
	NewBBox = BBox
	rounds = []
	while len(NewBBox.Ballots) > 0:
		res = DoRound(NewBBox)
		rounds.append(res)
		
		BottomCands = []
		VoteThreshold = ThresholdFunc([r[1] for r in res])
		for r in res:
			if r[1] <= VoteThreshold:
				BottomCands.append(r[0])	
		BottomCands.sort()
		BottomCands = tuple(BottomCands)
		
		NewBBox = RemoveCandidates(NewBBox, BottomCands)

	return tuple(rounds)

def Average(lst):
	return float(sum(lst))/float(len(lst))

# http://en.wikipedia.org/wiki/Single_transferable_vote

def SingleTransferableVote(BBox, Num):
	# The Hagenbach-Bischoff quota,
	# a fractional version of the Droop quota
	# http://en.wikipedia.org/wiki/Hagenbach-Bischoff_quota
	Quota = TotalCount(BBox)/(Num + 1.0)
	NewBBox = BBox
	rounds = []
	Winners = []
	while len(NewBBox.Ballots) > 0 and len(Winners) < Num:
		res = TopOne(NewBBox)
		
		# Check for winners
		# Find top candidates with the same Vote
		TopCands = []
		MaxVotes = res[0][1]
		for r in res:
			if r[1] != MaxVotes: break
			TopCands.append(r[0])
		TopCands.sort()
		TopCands = tuple(TopCands)
		
		ExcessVotes = MaxVotes - Quota
		if ExcessVotes >= 0:
			# Adjust weights to use winners' excess votes
			NewBallots = []
			for Weight, Votes in NewBBox.Ballots:
				if Votes[0] in TopCands:
					Weight = (Weight * ExcessVotes) / MaxVotes
				NewBallots.append((Weight,Votes))
			NewBBox = BallotBox(NewBallots)
			NewBBox = RemoveCandidates(NewBBox, TopCands)
			Winners.extend(TopCands)
			rounds.append(tuple(list(res) + [('+',TopCands)]))
			continue	
		
		# Check for losers
		# Find bottom ones with the same vote
		BottomCands = []
		MinVotes = res[-1][1]
		for r in reversed(res):
			if r[1] != MinVotes: break
			BottomCands.append(r[0])	
		BottomCands.sort()
		BottomCands = tuple(BottomCands)
		
		NewBBox = RemoveCandidates(NewBBox, BottomCands)
		rounds.append(tuple(list(res) + [('-',BottomCands)]))
	
	rounds.append((tuple(Winners),))
	
	return tuple(rounds)

# http://en.wikipedia.org/wiki/CPO-STV

def CompPairOutcomesSTV(BBox, Num):
	Quota = TotalCount(BBox)/(Num + 1.0)
	Cands = BBox.Candidates()
	n = len(Cands)
	
	# Find all possible sets of Num candidates
	Posses = [[]]
	for m in xrange(Num):
		NewPosses = []
		for Poss in Posses:
			base = max(Poss)+1 if len(Poss) > 0 else 0
			for k in xrange(base,n-Num+m+1):
				NewPoss = Poss + [k]
				NewPosses.append(NewPoss)
		Posses = NewPosses
	
	CandPosses = [[Cands[p] for p in Poss] for Poss in Posses]
	
	ncp = len(CandPosses)
	CprMat = [ncp*[0] for k in xrange(ncp)]
	
	CS = set(Cands)
	for k1 in xrange(ncp-1):
		Cands1 = CandPosses[k1]
		CS1 = set(Cands1)
		for k2 in xrange(k1+1,ncp):
			Cands2 = CandPosses[k2]
			CS2 = set(Cands2)
			InBoth = CS1 & CS2
			InNeither = CS - (CS1 | CS2)
			
			# Eliminate candidates in neither possibilities
			NewBBox = RemoveCandidates(BBox, InNeither)
			
			TotalOverThreshold = 0
			# Eliminate over-threshold candidates in both possibilities
			while True:
				res = TopOne(NewBBox)
				VotesOverThreshold = {}
				for r in res:
					if r[0] in InBoth and r[1] >= Quota:
						VotesOverThreshold[r[0]] = r[1]
				CandsOverThreshold = VotesOverThreshold.keys()
				NumOverThreshold = len(CandsOverThreshold)
				if NumOverThreshold <= 0: break
				TotalOverThreshold += NumOverThreshold
				
				NewBallots = []
				for Weight, Votes in NewBBox.Ballots:
					if Votes[0] in CandsOverThreshold:
						TotalVotes = VotesOverThreshold[Votes[0]]
						Weight = (Weight * (TotalVotes - Quota)) / TotalVotes
					NewBallots.append((Weight,Votes))
				NewBBox = BallotBox(NewBallots)
				NewBBox = RemoveCandidates(NewBBox, CandsOverThreshold)
			
			Total1 = Total2 = Quota*TotalOverThreshold
			for r in res:
				if r[0] in Cands1:
					Total1 += r[1]
				if r[0] in Cands2:
					Total2 += r[1]
			CprMat[k1][k2] = Total1
			CprMat[k2][k1] = Total2
	
	return (CandPosses, CprMat)

# http://en.wikipedia.org/wiki/Schulze_STV
# Not implemented because of the difficulty of understanding the algorithm
	

# Condorcet methods

# http://en.wikipedia.org/wiki/Condorcet_method

def Condorcet(BBox,CondorcetCandPM):
	Cands = BBox.Candidates()
	PrefMat = BBox.CondorcetMatrix()
	return CondorcetCandPM(Cands,PrefMat)

def CondorcetCandPMPrefOrder(Cands,PrefMat,PrefOrderFunction):
	PrefOrder = PrefOrderFunction(PrefMat)
	return tuple((Cands[PO] for PO in PrefOrder))

def CondorcetPrefOrder(BBox,PrefOrderFunction):
	return Condorcet(BBox, lambda Cands, PrefMat: \
		CondorcetCandPMPrefOrder(Cands,PrefMat,PrefOrderFunction))

def CondorcetCandPMPrefCount(Cands,PrefMat,PrefCountFunction):
	PrefCount = PrefCountFunction(PrefMat)
	CandCounts = zip(Cands, PrefCount)
	CandCounts.sort(CCSortFunction)
	return tuple(CandCounts)

def CondorcetPrefCount(BBox,PrefCountFunction):
	return Condorcet(BBox, lambda Cands, PrefMat: \
		CondorcetCandPMPrefCount(Cands,PrefMat,PrefCountFunction))


def CondorcetWinnerIndex(PrefMat):
	n = len(PrefMat)
	for i in xrange(n):
		wix = i
		for j in xrange(n):
			if j != i and PrefMat[i][j] <= PrefMat[j][i]:
				wix = None
				break
		if wix != None: return wix
	return None

def CondorcetWinner(BBox):
	Cands = BBox.Candidates()
	PrefMat = BBox.CondorcetMatrix()
	wix = CondorcetWinnerIndex(PrefMat)
	if wix != None:
		return (True,Cands[wix])
	else:
		return (False,)
	

# http://en.wikipedia.org/wiki/Schulze_method

def SchulzeOrdering(BPMat,i,j):
	rc = - cmp(BPMat[i][j],BPMat[j][i])
	if rc != 0: return rc
	return cmp(i,j)

def SchulzePrefOrder(PrefMat):
	n = len(PrefMat)
	
	# Beatpaths
	BPMat = [n*[0] for k in xrange(n)]
	
	# Variant of Floyd-Warshall algorithm
	for i in xrange(n):
		for j in xrange(n):
			if j != i:
				if PrefMat[i][j] > PrefMat[j][i]:
					BPMat[i][j] = PrefMat[i][j]
				else:
					BPMat[i][j] = 0
	
	for i in xrange(n):
		for j in xrange(n):
			if j != i:
				for k in xrange(n):
					if i != k and j != k:
						BPMat[j][k] = \
							max(BPMat[j][k], min(BPMat[j][i],BPMat[i][k]))
	
	Ordering = range(n)
	Ordering.sort(lambda i,j: SchulzeOrdering(BPMat,i,j))
	return Ordering

def Schulze(BBox): return CondorcetPrefOrder(BBox,SchulzePrefOrder)

# http://en.wikipedia.org/wiki/Copeland%27s_method

def CopelandPrefCount(PrefMat):
	n = len(PrefMat)
	
	Count = n*[0]
	for k1 in xrange(n):
		for k2 in xrange(n):
			Count[k1] += cmp(PrefMat[k1][k2],PrefMat[k2][k1])
	
	return Count

def Copeland(BBox): return CondorcetPrefCount(BBox,CopelandPrefCount)

# http://en.wikipedia.org/wiki/Minimax_condorcet

def MinimaxPrefCount(PrefMat,Which):
	n = len(PrefMat)
	
	Scores = n*[0]
	for k in xrange(n):
		Score = None
		for kx in xrange(n):
			if kx == k: continue
			
			if Which == "wins":
				if PrefMat[kx][k] > PrefMat[k][kx]:
					NewScore = PrefMat[kx][k]
				else:
					NewScore = 0
			elif Which == "marg":
				NewScore = PrefMat[kx][k] - PrefMat[k][kx]
			elif Which == "oppo":
				NewScore = PrefMat[kx][k]
			
			if Score == None:
				Score = NewScore
			elif NewScore > Score:
				Score = NewScore
		
		Scores[k] = - Score if n > 1 else 0
	
	return Scores

def Minimax(BBox,Which):
	return CondorcetPrefCount(BBox,lambda PM: MinimaxPrefCount(PM,Which))

# Permutation generator
# Uses in-place algorithm from
# http://en.wikipedia.org/wiki/Permutation#Generation_in_lexicographic_order
#
# Its arg is the data to permute
#
# It returns the permutations in iterator fashion,
# starting with ascending sort end ending with descending sort

def Permutations(data):
	# Set up the initial permutation
	perm = list(data)
	perm.sort()
	n = len(perm)

	while True:
		# Emit each permutation:
		yield perm
		
		# Find the next permutation
		# Quit if not possible
		ix0 = None
		for k in xrange(n-1):
			if perm[k] < perm[k+1]:
				ix0 = k
		if ix0 == None: break
		
		ix1 = ix0+1
		for k in xrange(ix0+1,n):
			if perm[ix0] < perm[k]:
				ix1 = k
		
		temp = perm[ix0]
		perm[ix0] = perm[ix1]
		perm[ix1] = temp
		
		for k in xrange(ix0+1,n):
			kr = (ix0+n) - k
			if kr <= k: break
			temp = perm[k]
			perm[k] = perm[kr]
			perm[kr] = temp

# http://en.wikipedia.org/wiki/Kemeny%E2%80%93Young_method

def KemenyYoungPrefOrder(PrefMat):
	n = len(PrefMat)
	
	# Find all permutations
	startperm = range(n)
	bestperm = startperm
	bestscore = None
	
	for perm in Permutations(startperm):
		score = 0
		for k1 in xrange(n-1):
			for k2 in xrange(k1+1,n):
				score += PrefMat[perm[k1]][perm[k2]]
		if bestscore == None:
			bestperm = tuple(perm)
			bestscore = score
		elif score > bestscore:
			bestperm = tuple(perm)
			bestscore = score
	
	return bestperm

def KemenyYoung(BBox): return CondorcetPrefOrder(BBox,KemenyYoungPrefOrder)


# This finds out what disjoint cycles a permutation has
def PermCycles(perm):
	n = len(perm)
	# Tag each element of the permutation with
	# whether it's available for including in a cycle
	Avail = n*[True]
	Cycles = []
	
	for i in xrange(n):
		if not Avail[i]: continue
		# Start a cycle
		CycleBegin = i
		Cycle = [CycleBegin]
		Avail[CycleBegin] = False
		NextIndex = i
		# Find the next element of a cycle
		# until one comes to the first element
		while True:
			NextIndex = perm[NextIndex]
			if NextIndex == CycleBegin: break
			Cycle.append(NextIndex)
			Avail[NextIndex] = False
		# Cycle complete
		Cycles.append(tuple(Cycle))
	
	return tuple(Cycles)


# https://en.wikipedia.org/wiki/Dodgson's_method

def CDLSortFunc(a,b):
	rc = cmp(a[1],b[1])
	if rc != 0: return rc
	return cmp(a[0],b[0])

def Dodgson(BBox):
	# Get the ballots
	Ballots = BBox.Ballots
	Cands = BBox.Candidates()
	n = len(Cands)
	
	# Fill out the ballots to full length
	# for doing permutations
	ExtendedBallots = []
	for wt, prefs in Ballots:
		extprefs = [(True,pref) for pref in prefs] + (n-len(prefs))*[(False,)]
		ExtendedBallots.append((wt,extprefs))
	
	# Permutation distance for each candidate:
	# initialize to null value
	CandDists = {}
	for Cand in Cands: CandDists[Cand] = n
	
	# Go through the permutations
	startperm = range(n);
	for perm in Permutations(startperm):
		# Get the permutation distance
		pcs = PermCycles(perm)
		permdist = len(perm) - len(pcs)
		
		# Permute the ballots
		PermBallots = []
		for wt, prefs in ExtendedBallots:
			extpermprefs = [prefs[p] for p in perm]
			permprefs = [pref[1] for pref in extpermprefs if pref[0]]
			PermBallots.append((wt,permprefs))
		PermBBox = BallotBox(PermBallots)
		
		# Find the Condorcet winner (if any) and use it
		CWin = CondorcetWinner(PermBBox)
		if CWin[0]:
			Cand = CWin[1]
			if permdist < CandDists[Cand]:
				CandDists[Cand] = permdist
	
	CDList = [(Cand, CandDists[Cand]) for Cand in Cands]
	CDList.sort(CDLSortFunc)
	return tuple(CDList)


# http://en.wikipedia.org/wiki/Ranked_pairs
# Uses cyclicity test in
# http://www.cs.hmc.edu/~keller/courses/cs60/s98/examples/acyclic/

def RankedPairsCompare(a,b):
	rc = - cmp(a[2],b[2])
	if rc != 0: return rc
	
	rc = cmp(a[3],b[3])
	return rc

def RankedPairsPrefOrder(PrefMat):
	n = len(PrefMat)

	# Find the beat margins and sort them
	BeatMargins = []
	for k1 in xrange(n):
		for k2 in xrange(n):
			if k2 != k1:
				bmg = (k1, k2, PrefMat[k1][k2],PrefMat[k2][k1])
				BeatMargins.append(bmg)
	BeatMargins.sort(RankedPairsCompare)
	
	# Find the beat list of what beats what
	BeatList = []
	for bmg in BeatMargins:
		# Get a pair of (start index, end index)
		bd = tuple(bmg[:2])
		
		# The simplest test for cycles -- reverse direction
		rvbd = (bd[1],bd[0])
		if rvbd in BeatList: continue
		
		# Nontrivial test: remove all the nodes that are
		# only on left sides or on right sides
		# Repeat until it is no longer possible to remove any more
		# Are any left?
		BLX = BeatList + [bd]
		while True:
			# What is present on each side?
			NumLeft = n*[0]
			NumRight = n*[0]
			for bdx in BLX:
				NumLeft[bdx[0]] += 1
				NumRight[bdx[1]] += 1
			
			# What is present on only one side?
			RemoveLeft = []
			RemoveRight = []
			for k in xrange(n):
				if NumLeft[k] > 0 and NumRight[k] == 0:
					RemoveLeft.append(k)
				if NumLeft[k] == 0 and NumRight[k] > 0:
					RemoveRight.append(k)
			
			# Can one go any further?
			if len(RemoveLeft) + len(RemoveRight) == 0: break
			
			# Remove!
			BLY = []
			for bdx in BLX:
				if bdx[0] in RemoveLeft: continue
				if bdx[1] in RemoveRight: continue
				BLY.append(bdx)
			BLX = BLY
		
		# If no cycles found, then add the new one
		if len(BLX) == 0:
			BeatList.append(bd)
	
	CandIxs = range(n)
	CandIxs.sort(lambda a,b: -1 if (a,b) in BeatList else (0 if a == b else 1))
	return CandIxs

def RankedPairs(BBox): return CondorcetPrefOrder(BBox,RankedPairsPrefOrder)


# https://en.wikipedia.org/wiki/Maximal_lotteries


# Use linear programming

# Solve a system of linear equations with the Gauss-Jordan algorithm
# The input has form (matrix) (vector, vector, ...)
def SolveLinearEquations(mat):
	n = len(mat)
	nx = len(mat[0])
	workmat = [[float(mat[i][j]) for j in xrange(nx)] for i in xrange(n)]
	
	# Do forward substitution
	for icol in xrange(n):
		# Necessary to exchange rows
		# to bring a nonzero value into position?
		# Return None if singular
		if workmat[icol][icol] == 0:
			ipvt = None
			for i in xrange(icol+1,n):
				if workmat[i][icol] != 0:
					ipvt = i
					break
			if ipvt == None: return None
			temp = workmat[icol]
			workmat[icol] = workmat[ipvt]
			workmat[ipvt] = temp
		# Make diagonal 1:
		wmicol = workmat[icol]
		dgvrecip = 1/wmicol[icol]
		for i in xrange(icol,nx):
			wmicol[i] *= dgvrecip
		# Forward substitute:
		for i in xrange(icol+1,n):
			wmi = workmat[i]
			elimval = wmi[icol]
			for j in xrange(icol,nx):
				wmi[j] -= elimval*wmicol[j]
	
	# Do back substitution
	for icol in xrange(n-1,0,-1):
		wmicol = workmat[icol]
		for i in xrange(icol):
			wmi = workmat[i]
			elimval = wmi[icol]
			for j in xrange(icol,nx):
				wmi[j] -= elimval*wmicol[j]
	
	# Done!
	return [[workmat[i][j] for j in xrange(n,nx)] for i in xrange(n)]

def LP_SimplexMethodStep(basvars, tableau):
	nrows = len(tableau)
	ncols = len(tableau[0])
	
	# Find pivot column
	pvcol = -1
	objf = tableau[0]
	for k in xrange(ncols-1):
		newval = objf[k]
		if pvcol < 0:
			pvcol = k
			val = newval
		elif newval < val:
			pvcol = k
			val = newval
	
	# Check for imperfect numerical cancellation
	if val >= 0: return (False,"Complete")	
	
	# Find pivot row
	pvrow = -1
	for k in xrange(1,nrows):
		dvsr = tableau[k][pvcol]
		if dvsr > 0:
			newval = tableau[k][-1]/dvsr
			if newval >= 0:
				if pvrow < 0:
					pvrow = k
					val = newval
				elif newval < val:
					pvrow = k
					val = newval
	if pvrow < 0: return (False, "Unsolvable - Improvement Step")
	
	basvars[pvrow] = pvcol
	
	xrow = tableau[pvrow]
	xrc = xrow[pvcol]
	for k in xrange(ncols):
		xrow[k] /= xrc
	for j in xrange(nrows):
		if j == pvrow: continue
		yrow = tableau[j]
		yrc = yrow[pvcol]
		for k in xrange(ncols):
			yrow[k] -= yrc*xrow[k]
	
	# Fix the numerical values
	for k in xrange(nrows):
		tableau[j][pvcol] = 1. if k == pvrow else 0.
	
	# Make approximate cancellations exact
	objf = tableau[0]
	objfmax = max(max(objf),-min(objf))
	eps = (1e-12)*objfmax
	for k in xrange(ncols-2):
		if abs(objf[k]) <= eps: objf[k] = 0.
	
	return (True,"Successful Step")

def LP_SimplexMethodSeq(basvars, tableau):
	while True:
		donext, res = LP_SimplexMethodStep(basvars, tableau)
		if not donext:
			return res

def LP_SimplexMethod(tableau_):
	nrows = len(tableau_)
	ncols = len(tableau_[0])
	tableau = [[float(tableau_[i][j]) for j in xrange(ncols)] for i in xrange(nrows)]
	
	basvars = nrows*[-1]
	basvars[0] = ncols-2
	allzero = False
	for i in xrange(ncols-2):
		numnz = 0
		for j in xrange(1,nrows):
			val = tableau[j][i]
			if val != 0:
				numnz += 1
				nzrow = j
		if numnz == 0:
			allzero = True
			break
		elif numnz == 1:
			if basvars[nzrow] < 0: basvars[nzrow] = i
	
	if allzero: return (False, "Unsolvable - Some Variables Unspecified")
	
	rc = LP_SimplexMethodSeq(basvars, tableau)
	if rc != "Complete": return (False,rc)
	
	# Find the nonzero variables, all but the final slack variable
	basvars = [b for b in basvars[1:] if b >= 0]
	objf = tableau[0]
	nzvrixs = [k for k in xrange(ncols-2) if objf[k] == 0]
	nzvrdiff = list(set(nzvrixs) - set(basvars))
	nzvrdiff.sort()
	nzvrixs = basvars + nzvrdiff
	if len(nzvrixs) < (nrows-1): return (False, "Unsolvable")
	
	# Handle degenerate cases by using the first of the variables
	nzvrixs = nzvrixs[:nrows-1]
	
	# Solve a system of linear equations by
	# Gaussian elimination with pivoting
	slmat = [nrows*[0] for k in xrange(nrows-1)]
	for j in xrange(nrows-1):
		slmrow = slmat[j]
		tblrow = tableau[j+1]
		for k in xrange(nrows-1):
			slmrow[k] = tblrow[nzvrixs[k]]
	for k in xrange(nrows-1):
		slmat[k][-1] = tableau[k+1][-1]
	res = SolveLinearEquations(slmat)
	
	vals = (ncols-1)*[0]
	for k in xrange(nrows-1):
		vals[nzvrixs[k]] = res[k][0]
	
	# Solve for the final slack variable
	vsum = objf[ncols-1]
	for k in xrange(ncols-2):
		vsum -= objf[k]*vals[k]
	vals[ncols-2] = vsum/objf[ncols-2]
	
	return (True, vals)


def MaximalLotteriesPrefCount_Simplex(PrefMat):
	n = len(PrefMat)
	
	if n == 0: return ()
	
	tableau = [(2*n+3)*[0] for k in xrange(n+2)]
	
	for i in xrange(n):
		for j in xrange(n):
			tableau[i+2][j] = PrefMat[i][j] - PrefMat[j][i]
	
	for i in xrange(n):
		tableau[1][i] = 1
		tableau[i+2][n] = 1
		tableau[i+2][i+n+1] = 1
	
	tableau[0][n] = -1
	tableau[0][2*n+1] = 1
	tableau[1][2*n+2] = 1
	
	solved, res = LP_SimplexMethod(tableau)
	
	# Failure mode: return all zeros
	if solved:
		Scores = res[:n]
	else:
		Scores = n*[0]
	
	return tuple(Scores)


def MultMatVec(mat,vec):
	n = len(vec)
	
	res = n*[0]
	for i in xrange(n):
		mr = mat[i]
		s = 0
		for j in xrange(n):
			s += mr[j]*vec[j]
		res[i] = s

	return res

def LineCalc(pt,dst,dir):
	n = len(pt)
	return [pt[i] + dst*dir[i] for i in xrange(n)]

def FixScores(Scores):
	NewScores = [max(s,0) for s in Scores]
	NSSR = 1.0/sum(NewScores)
	return [NSSR*s for s in NewScores]

# How many iterations and divisions for this algorithm

# How many score-finding iterations
mlitrrpt = 16
# How many iterations to improve a score
mlitrpos = 16
# How many iterations to find a good direction from a score value
mlitrdir = 16
# How many divisions in a line search along a direction
mldvsrch = 16

def MaximalLotteriesPrefCount_RandomSearch(PrefMat):
	n = len(PrefMat)
	
	Scores = n*[0]
	
	# Since I do not have a "real" linear-programming solver
	# available to me and callable from Python, I will write my own,
	# a stochastic line-search algorithm
	
	if n == 0: return Scores
	
	if n == 1:
		Scores[0] = 1
		return Scores
	
	# Preference matrix minus its transpose
	PMAS = [[PrefMat[j][i] - PrefMat[i][j] for j in xrange(n)] for i in xrange(n)]
	
	# Initial point: place in the center of the allowed region
	# It may or may not be allowed by the preference-matrix constraints
	SCAvg = 1.0/n
	Scores = n*[SCAvg]
	
	# Find initial constraint values:
	CstrVals = MultMatVec(PMAS,Scores)
	# Value to maximize
	CstrMin = min(CstrVals)
	
	# Intended final values
	NewScores = Scores
	NewCstrMin = CstrMin
	
	import random
	
	# Try several overall runs
	for iterrpt in xrange(mlitrrpt):

		# Start a point
		NewPosScores = Scores
		NewPosCstrMin = CstrMin
		
		# Iterate over advancement steps
		for iterpos in xrange(mlitrpos):
			
			# Start a direction search
			NewDirScores = NewPosScores
			NewDirCstrMin = NewPosCstrMin
		
			# Iterate over directions
			for iterdir in xrange(mlitrdir):
				
				# Find a random direction
				# Be sure that no component is zero
				Dir = [random.gauss(0,1) for i in xrange(n)]
				DirAvg = sum(Dir)/n
				Dir = [Dir[i] - DirAvg for i in xrange(n)]
				for d in Dir:
					if d == 0: continue
				
				# Find the interval that is allowable
				# by the scores being nonnegative and adding up to 1
				PMin = None
				PMax = None
				for i in xrange(n):
					if Dir[i] >= 0:
						NMin = - NewPosScores[i]/Dir[i]
						NMax = (1 - NewPosScores[i])/Dir[i]
					else:
						NMin = (1 - NewPosScores[i])/Dir[i]
						NMax = - NewPosScores[i]/Dir[i]
					
					if PMin == None: PMin = NMin
					elif NMin > PMin: PMin = NMin
					
					if PMax == None: PMax = NMax
					elif NMax < PMax: PMax = NMax
					
				if PMin == None: continue
				if PMax == None: continue
			
				# Do a line search between these points
				LSDivs = mldvsrch
				LSDRcp = 1./LSDivs
				for i in xrange(0,LSDivs+1):
					Dst = PMin + (PMax - PMin)*LSDRcp*i
					NewSrchScores = FixScores(LineCalc(NewPosScores, Dst, Dir))
					CstrVals = MultMatVec(PMAS,NewSrchScores)
					NewSrchCstrMin = min(CstrVals)
					
					# Improved?
					if NewSrchCstrMin > NewDirCstrMin:
						NewDirScores = NewSrchScores
						NewDirCstrMin = NewSrchCstrMin
			
			# Improved direction-search result?
			# If so, then accept it as the next point
			if NewDirCstrMin > NewPosCstrMin:
				NewPosScores = NewDirScores
				NewPosCstrMin = NewDirCstrMin
		
		# Improved final point?
		if NewPosCstrMin > NewCstrMin:
			NewScores = NewPosScores
			NewCstrMin = NewPosCstrMin

	return NewScores

def MaximalLotteriesPrefCount(PrefMat):
	Scores = MaximalLotteriesPrefCount_Simplex(PrefMat)
		
	# Did the simplex method return a nonsensical result?
	ResOK = True
	eps = 1e-12
	smin = - eps
	smax = 1 + eps
	for s in Scores:
		if s < smin or s > smax:
			ResOK = False
			break

	# Did the simplex method fail?
	# Failure mode: returning all 0's
	if ResOK:
		if len(Scores) > 0:
			if max(Scores) == 0: ResOK = False
		
	if ResOK: return Scores
		
	# Trying the random-line-search interior-point method instead
	return MaximalLotteriesPrefCount_RandomSearch(PrefMat)


def MaximalLotteries(BBox):
	return CondorcetPrefCount(BBox, MaximalLotteriesPrefCount)


# Find the Smith and Schwartz sets
# From http://wiki.electorama.com/wiki/Maximal_elements_algorithms

def PrefMatRelations(PrefMat, Type):
	n = len(PrefMat)
	Relations = [n*[False] for k in xrange(n)]
	
	for k1 in xrange(n):
		for k2 in xrange(n):
			if k2 != k1:
				if Type == "Schwartz":
					Relations[k1][k2] = PrefMat[k1][k2] > PrefMat[k2][k1]
				elif Type == "Smith":
					Relations[k1][k2] = PrefMat[k1][k2] >= PrefMat[k2][k1]
	
	return Relations


def FloydWarshallMaximal(Relations):
	n = len(Relations)
	InMaximal = n*[True]
	
	# Init HasPath for relations with length 1
	HasPath = [n*[False] for k in xrange(n)]
	for k1 in xrange(n):
		for k2 in xrange(n):
			if k2 != k1:
				HasPath[k1][k2] = Relations[k1][k2]
	
	# Consider paths with intermediate nodes from 1 to k
	for k in xrange(n):
		for k1 in xrange(n):
			if k1 != k:
				for k2 in xrange(n):
					if k2 != k and k2 != k1:
						if HasPath[k1][k] and HasPath[k][k2]:
							HasPath[k1][k2] = True
	
	# Candidates with paths to them but none to complete a cycle,
	# they are not in the maximal set
	for k1 in xrange(n):
		for k2 in xrange(n):
			if k2 != k1:
				if HasPath[k2][k1] and not HasPath[k1][k2]:
					InMaximal[k1] = False
	
	return InMaximal

def MaximalSet(BBox, Type):
	Cands = BBox.Candidates()
	PrefMat = BBox.CondorcetMatrix()
	
	Rels = PrefMatRelations(PrefMat, Type)
	InMax = FloydWarshallMaximal(Rels)
	
	return tuple((cm[0] for cm in zip(Cands,InMax) if cm[1]))

def MaximalSetSequence(BBox, Type):
	MaximalSets = []
	NewBBox = BBox
	
	while True:
		MaxSet = MaximalSet(NewBBox, Type)
		if len(MaxSet) == 0: break
		MaximalSets.append(MaxSet)
		NewBBox = RemoveCandidates(NewBBox, MaxSet)
	
	return tuple(MaximalSets)


# For debugging

def DumpSingleAlgorithm(Algorithm,BBox):

	# For reference
	print "Schulze and Tideman:"
	print Schulze(BBox)
	print RankedPairs(BBox)
	print
	
	print "Algorithm to Test:"
	res = Algorithm(BBox)
	print res
	print


def DumpAll(Ballots, DoNFactorial=True):
	BBox = BallotBox(Ballots)

	# DumpSingleAlgorithm(MajorityJudgment,BBox); return
	
	print "Candidates:",
	Cands = BBox.Candidates()
	NumCands = len(Cands)
	for Cand in Cands: print Cand, " ",
	print
	print
	
	print "Total Count:", TotalCount(BBox)
	print
	
	print "All Votes:"
	res = AllVotes(BBox)
	for r in res: print r
	print
	
	print "Top One (First Past the Post, Plurality):"
	res = TopOne(BBox)
	for r in res: print r
	print
	
	print "All Top N's (Bucklin):"
	reslist = TopNumList(BBox)
	for k, res in enumerate(reslist):
		print "N =", k+1
		for r in res: print r
	print
	
	print "Borda Count:"
	res = Borda(BBox)
	for r in res: print r
	print
	
	print "Modified Borda Count:"
	res = ModBorda(BBox)
	for r in res: print r
	print
	
	print "Cumulative Borda Count:"
	res = CumulBorda(BBox)
	for r in res: print r
	print
	
	print "Majority Judgment:"
	res = MajorityJudgment(BBox)
	for r in res: print r
	print
	
	print "Top-Two Runoff:"
	reslist = TopTwoRunoff(BBox)
	for k, res in enumerate(reslist):
		print "Round", k+1
		for r in res: print r
	print		

	print "Sequential Runoff:"
	reslist = SequentialRunoff(BBox)
	for k, res in enumerate(reslist):
		print "Round", k+1
		for r in res: print r
	print
	
	print "Baldwin:"
	reslist = SequentialRunoff(BBox,min,Borda)
	for k, res in enumerate(reslist):
		print "Round", k+1
		for r in res: print r
	print
	
	print "Nanson:"
	reslist = SequentialRunoff(BBox,Average,Borda)
	for k, res in enumerate(reslist):
		print "Round", k+1
		for r in res: print r
	print
	
	print "Coombs:"
	reslist = SequentialRunoff(BBox,min,MinusBottom)
	for k, res in enumerate(reslist):
		print "Round", k+1
		for r in res: print r
	print

	for k in xrange(min(3,NumCands)):
		Num = k+1
		print "Single Transferable Vote:", Num
		reslist = SingleTransferableVote(BBox,Num)
		for k, res in enumerate(reslist):
			print "Round", k+1
			for r in res: print r
		print

	for k in xrange(min(3,NumCands)):
		Num = k+1
		print "Comparison of Pairs of Outcomes by STV:", Num
		CandSets, PrefMat = CompPairOutcomesSTV(BBox,Num)
		CandSetsOrdered = CondorcetCandPMPrefOrder(CandSets, PrefMat, SchulzePrefOrder)
		if len(CandSetsOrdered) > 0:
			print CandSetsOrdered[0]
		print
	
	print "Condorcet Matrix:"
	PrefMat = BBox.CondorcetMatrix()
	for PMRow in PrefMat:
		for PMVal in PMRow:
			print PMVal, '\t',
		print
	print
	
	print "Condorcet Winner:"
	print CondorcetWinner(BBox)
	print
	
	print "Schulze Beatpath:"
	print Schulze(BBox)
	print
	
	print "Copeland Pairwise Aggregation:"
	res = Copeland(BBox)
	for r in res: print r
	print
	
	for Which in ("wins", "marg", "oppo"):
		print "Minimax:", Which
		res = Minimax(BBox,Which)
		for r in res: print r
		print
	
	print "Kemeny-Young:"
	if DoNFactorial:
		print KemenyYoung(BBox)
	else:
		print "Skipped because it is O(n!)"
	print
	
	print "Dodgson:"
	if DoNFactorial:
		print Dodgson(BBox)
	else:
		print "Skipped because it is O(n!)"
	print
	
	print "Tideman Ranked Pairs:"
	print RankedPairs(BBox)
	print
	
	print "Maximal Lotteries:"
	res = MaximalLotteries(BBox)
	for r in res: print r
	print
	
	print "Maximal sets: Smith and Schwartz sets:"
	MaxSetSmith = MaximalSet(BBox,"Smith")
	MaxSetSchwartz = MaximalSet(BBox,"Schwartz")
	BBoxSmith = KeepCandidates(BBox, MaxSetSmith)
	BBoxSchwartz = KeepCandidates(BBox, MaxSetSchwartz)
	print MaxSetSmith
	print MaxSetSchwartz
	print
	
	print "Borda on Smith Set:"
	res = Borda(BBoxSmith)
	for r in res: print r
	print
	
	print "Borda on Schwartz Set:"
	res = Borda(BBoxSchwartz)
	for r in res: print r
	print
	
	print "IRV on Smith Set:"
	reslist = SequentialRunoff(BBoxSmith)
	for k, res in enumerate(reslist):
		print "Round", k+1
		for r in res: print r
	print
	
	print "IRV on Schwartz Set:"
	reslist = SequentialRunoff(BBoxSchwartz)
	for k, res in enumerate(reslist):
		print "Round", k+1
		for r in res: print r
	print
	
	print "Sequence of maximal sets: Smith and Schwartz sets:"
	print MaximalSetSequence(BBox,"Smith")
	print MaximalSetSequence(BBox,"Schwartz")
	print
	
	print
	
	
# Debug:
if __name__ == "__main__":
	
	# Edge cases:
	ZeroBallots = ()
	OneBallot = ((1, ("The Dear Leader",)),)
	TwoBallots = ((3, ("T Dee", "T Dum")), (2, ("T Dum", "T Dee")))
	
	# The beer and distilled-spirits examples are from
	# http://web.math.princeton.edu/math_alive/Voting/Lab1.shtml
	
	BeerBallotSource = (
		("Killians", "Molson", "Samuel Adams", "Guinness", "Meister Brau"),
		(
			(18, (5,1,2,4,3)), (12, (1,5,3,4,2)), (10, (2,5,4,1,3)),
			(9, (4,5,1,2,3)), (4, (2,5,3,4,1)), (2, (4,5,3,2,1))
		)
	)
	BeerBallots = CWLPrefNumsToOrder(BeerBallotSource)
	
	DistilledBallotSource = (
		("Absolut", "Jim Beam", "Jose Cuervo", "Jack Daniels", "So. Comfort"),
		(
			(5, (1,2,3,4,5)), (2, (5,1,2,3,4)), (3, (5,2,1,3,4)),
			(3, (5,2,3,1,4)), (4, (5,2,3,4,1))
		)
	)
	DistilledBallots = CWLPrefNumsToOrder(DistilledBallotSource)
	
	
	# Several Wikipedia articles use these ones
	
	TennesseeCapitalBallots = (
		(0.42, ("Memphis","Nashville","Chattanooga","Knoxville")),
		(0.26, ("Nashville","Chattanooga","Knoxville","Memphis")),
		(0.15, ("Chattanooga","Knoxville","Nashville","Memphis")),
		(0.17, ("Knoxville","Chattanooga","Nashville","Memphis"))
	)
	
	ABCDSBallots = (
		(25, ('Andrea',)),
		(34, ('Carter','Brad','Delilah')),
		(7, ('Brad','Delilah')),
		(8, ('Delilah','Brad')),
		(5, ('Delilah','Scott')),
		(21, ('Scott', 'Delilah'))
	)
	
	# My one: ice-cream flavors
	
	ICFBallots = (('Vanilla', 'Neapolitan', 'Chocolate', 'Cookie Dough'),
		('Rocky Road', 'Mint CC', 'Chocolate', 'Vanilla'),
		('Strawberry', 'Cherry'),
		('Cookie Dough', 'Chocolate'),
		('Cherry', 'Mint CC', 'Rocky Road'),
		('Chocolate', 'Strawberry', 'Vanilla', 'Neapolitan'),
		('Vanilla', 'Cookie Dough', 'Rocky Road'),
		('Cherry', 'Strawberry'),
		('Neapolitan', 'Rocky Road', 'Mint CC'),
		('Cherry', 'Vanilla', 'Strawberry'),
		('Chocolate', 'Rocky Road', 'Mint CC'),
		('Chocolate', 'Vanilla', 'Strawberry'),
		('Cherry', 'Mint CC'),
		('Strawberry', 'Mint CC')
		)
	
	print "*** Empty Ballot Box ***"
	DumpAll(ZeroBallots)
	
	print "*** One-Party Election ***"
	DumpAll(OneBallot)
	
	print "*** Two-Party Election ***"
	DumpAll(TwoBallots)
	
	print "*** Beers ***"
	DumpAll(BeerBallots)
	
	print "*** Distilled Spirits ***"
	DumpAll(DistilledBallots)
	
	print "*** Tennessee Capital ***"
	DumpAll(TennesseeCapitalBallots)
	
	print "*** ABCDS ***"
	DumpAll(ABCDSBallots)
	
	print "*** Ice-Cream Flavors ***"
	DumpAll(MakeWeighted(ICFBallots))