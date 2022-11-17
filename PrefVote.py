#!python3
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
# Dowdall(BBox): like Borda, uses a different weighting:
# 1, 1/2, 1/3, 1/4, ...
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
# WinnerDropping(BBox, DoRound, GetWinner): returns a list of
# DoRound results, with the winner of each one dropped from the count
# Needs GetWinner for getting the winner from each round's results
#
# TopTwoRunoff(BBox, DoFirstRound=TopOne, DoSecondRound=TopOne):
# The first round selects the top two candidates,
# and in the second round, they go head to head.
#
# For STAR voting, use Borda or some other ranks-to-ratings method
# for the first round.
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
# CPOSTV(BBox, Num): fill N seats (N is second arg)
# Comparison of Pairs of Outcomes by Single Transferable Vote
# Finds all the possible outcomes and does pairwise comparisons,
# creating a Condorcet matrix for them
# This code uses the Schulze method to find the winning set
#
# SchulzeSTV(BBox, Num): fill N seats (N is second arg)
# Schulze STV compares sets of candidates that differ
# by only one candidate, creating a Condorcet matrix for them
# This code uses the Schulze method to find the winning set
#
#
# Condorcet methods
#
# CondorcetWinnerIndex(PrefMat)
# Returns index if there is one, None otherwise
#
# CondorcetWinner(BBox)
# Returns (True,winner) if there is one, (False,) otherwise
#
# CondorcetLoserIndex(PrefMat)
# Returns index if there is one, None otherwise
#
# CondorcetLoser(BBox)
# Returns (True,loser) if there is one, (False,) otherwise
#
# CondorcetBorda(BBox)
# Does a version of the Borda count with the Condorcet matrix
#
# CondorcetFlipBorda(BBox)
# Does a flipped version, with minus the transpose,
# to rank by losing instead of winning
#
# CondorcetCmbnBorda(BBox)
# Combines both of the previous ratings for sorting, winning then losing
#
# CondorcetWithFallback(BBox, fallback=Borda, settype="All")
# Black's method. It tries to find the Condorcet winner,
# and if it fails to do so, then falls back on another method
# evaluated on some subset of the candidates: "All", "Smith", "Schwartz".
# Returns (True, Condorcet winner) or (False, fallback-method output)
#
# CondorcetSequentialRunoff(BBox, CWn=False, CLs=False,
#   ThresholdFunc=min, DoRound=TopOne):
# Like SequentialRunoff, with the addition of
# CWn: use the Condorcet winner
# CLs: drop the Condorcet loser
#
# Copeland(BBox)
# Copeland's pairwise-aggregation method
# Returns sorted list of (candidate, score)
#
# Schulze(BBox)
# Schulze's beatpath method
# Returns sorted list of (candidate, score)
# The score is calculated from the beatpath-order matrix in Copeland fashion
#
# Minimax(BBox, Which): second arg is which sort:
# "wins" -- winning votes
# "marg" -- margins
# "oppo" -- pairwise opposition
# Returns sorted list of (candidate, score)
#
# RankedPairs(BBox, Which="oppo")
# Tideman's ranked-pairs method
# Like Kemeny-Young, but with simple hill-climbing optimization
# Returns sorted list of (candidate, score)
# The score is calculated from the pair-ranking matrix in Copeland fashion
#
# Which is like for Minimax.
# Tideman's original method has "oppo" (the default)
# Maximize Affirmed Majorities has "wins"
#
# Maximum Affirmed Majorities and Maximum Majority Voting
# seem identical to it
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
# Makes a sequence of maximal sets by removing them
# from the candidates as it goes
#
# MaximalSetSequentialRunoff(BBox, Type, CWn=False, CLs=False,
#   ThresholdFunc=min, DoRound=TopOne):
# Like SequentialRunoff, with the addition of the type of the maximal set:
# "Smith" or "Schwartz".
#
# Extra Methods
#
# DescendingSolidCoalitions(BBox)
# Descending Solid Coalitions
# Equivalent to Descending Acquiescing Coalitions
# because there are no tied preferences here
# Returns its winner(s)

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
		self._CondorcetMatrix = [NumCands*[0] for k in range(NumCands)]
		
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
			for k1 in range(NumCands):
				for k2 in range(NumCands):
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

# Threads negation over a list if need be
def negate(x):
	if hasattr(x,'__iter__'):
		return tuple(map(negate,x))
	else:
		return -x

# x is (Cand, Count / List of Counts)
# Sort reverse by count, then forward by cand
def CCSortFunction(x):
	return (negate(x[1]),x[0])

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
	CandCounts = list(zip(Cands,CountList))
	CandCounts.sort(key=CCSortFunction)
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
	for k in range(NumAdj):
		Counts[Votes[k]] += Weight

def TopNum(BBox,Num):
	return BBoxCounter(BBox, lambda Cn,Wt,Vt,Ca: TopNumFunction(Cn,Wt,Vt,Ca,Num))

def TopNumList(BBox):
	Cands = BBox.Candidates()
	NumCands = len(Cands)
	return tuple((TopNum(BBox,k+1) for k in range(NumCands)))


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

def DowdallFunction(Counts, Weight, Votes, Cands):
	MaxNum = len(Cands)
	for k,Vote in enumerate(Votes):
		Counts[Vote] += Weight/float(k+1)

def Dowdall(BBox):
	return BBoxCounter(BBox, DowdallFunction)

# https://en.wikipedia.org/wiki/Majority_judgment
def WeightedMedian(wtrnks):
	# Sort by ranks and find cumulative weights
	wrsort = list(wtrnks)
	wrsort.sort(key=lambda x: x[0])
	
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
	for k in range(n):
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
	
	CandWMs.sort(key=CCSortFunction)
	return tuple(CandWMs)


# Multistep methods

def WinnerDropping(BBox, DoRound, GetWinner):
	NewBBox = BBox
	rounds = []
	while len(NewBBox.Ballots) > 0:
		res = DoRound(NewBBox)
		rounds.append(res)
		NewBBox = RemoveCandidates(NewBBox, [GetWinner(res)])
	return tuple(rounds)

# http://en.wikipedia.org/wiki/Contingent_vote
# http://en.wikipedia.org/wiki/Two-round_system

def TopTwoRunoff(BBox, DoFirstRound=TopOne, DoSecondRound=TopOne):
	rounds = []
	res = DoFirstRound(BBox)
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
	res = DoSecondRound(NewBBox)
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
	Losers = []
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
		
		for BotCand in BottomCands:
			Losers.append(BotCand)
		
		NewBBox = RemoveCandidates(NewBBox, BottomCands)
		rounds.append(tuple(list(res) + [('-',BottomCands)]))
		
		if len(Winners) < Num:
			for loser in reversed(Losers):
				Winners.append(loser)
				if len(Winners) >= Num: break
		
	rounds.append((tuple(Winners),))
	
	return tuple(rounds)

# http://en.wikipedia.org/wiki/CPO-STV

# Subsets with length sslen
# of a range of integers from 0 to rnglen-1
def RangeSubsets(rnglen, sslen):
	RSS = [[]]
	for ns in range(sslen):
		NewRSS = []
		for SS in RSS:
			base = max(SS)+1 if len(SS) > 0 else 0
			for k in range(base,rnglen-sslen+ns+1):
				NewSS = SS + [k]
				NewRSS.append(NewSS)
		RSS = NewRSS
	return RSS

# Subsets with length sslen
# of a list
def ListSubsets(lst, sslen):
	rss = RangeSubsets(len(lst),sslen)
	return [[lst[n] for n in ss] for ss in rss]

# The main event
# Returns sets of candidates and a Condorcet matrix for them
def CPOSTVMatrix(BBox, Num):
	Quota = TotalCount(BBox)/(Num + 1.0)
	Cands = BBox.Candidates()
	
	# Find all possible sets of Num candidates
	CandPosses = ListSubsets(Cands,Num)
	
	ncp = len(CandPosses)
	CprMat = [ncp*[0] for k in range(ncp)]
	
	CS = set(Cands)
	for k1 in range(ncp-1):
		Cands1 = CandPosses[k1]
		CS1 = set(Cands1)
		for k2 in range(k1+1,ncp):
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

# Much like CPO-STV, but comparing sets of candidates that differ by only one member

# Returns sets of candidates and a Condorcet matrix for them
def SchulzeSTVMatrix(BBox,Num):
	Cands = BBox.Candidates()

	# Find all possible sets of Num candidates
	CandPosses = ListSubsets(Cands,Num)
	
	# Turn into sets for convenience
	CandPossSets = list(map(set,CandPosses))
	
	# Find the score of each set relative to each other set
	ncp = len(CandPosses)
	ScoreMat = [ncp*[0] for k in range(ncp)]
	for i in range(ncp):
		CPOrig = CandPosses[i]
		for j in range(ncp):
			scnddiff = list(CandPossSets[j] - CandPossSets[i])
			if len(scnddiff) != 1: continue
			sdval = scnddiff[0]
			ncpogt = Num*[0]
			for Ballot in BBox.Ballots:
				wt, prefs = Ballot
				if sdval not in prefs: continue
				sdvix = prefs.index(sdval)
				for k, cpoval in enumerate(CPOrig):
					if cpoval not in prefs: continue
					cpovix = prefs.index(cpoval)
					if cpovix < sdvix:
						ncpogt[k] += wt
			ScoreMat[i][j] = min(ncpogt)

	return (CandPosses, ScoreMat)


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
	CandCounts = list(zip(Cands, PrefCount))
	CandCounts.sort(key=CCSortFunction)
	return tuple(CandCounts)

def CondorcetPrefCount(BBox,PrefCountFunction):
	return Condorcet(BBox, lambda Cands, PrefMat: \
		CondorcetCandPMPrefCount(Cands,PrefMat,PrefCountFunction))

def CondorcetWinnerIndex(PrefMat):
	n = len(PrefMat)
	for i in range(n):
		wix = i
		for j in range(n):
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
	
def CondorcetLoserIndex(PrefMat):
	n = len(PrefMat)
	for i in range(n):
		lix = i
		for j in range(n):
			if j != i and PrefMat[i][j] >= PrefMat[j][i]:
				lix = None
				break
		if lix != None: return lix
	return None

def CondorcetLoser(BBox):
	Cands = BBox.Candidates()
	PrefMat = BBox.CondorcetMatrix()
	lix = CondorcetLoserIndex(PrefMat)
	if lix != None:
		return (True,Cands[lix])
	else:
		return (False,)

def CondorcetBordaCount(PrefMat):
	return tuple(map(sum,PrefMat))

def CondorcetBorda(BBox):
	return CondorcetPrefCount(BBox,CondorcetBordaCount)

def CondorcetFlipBordaCount(PrefMat):
	n = len(PrefMat)
	cnts = n*[0]
	for i in range(n):
		cnt = 0
		for j in range(n):
			cnt += PrefMat[j][i]
		cnts[i] = - cnt
	return tuple(cnts)

def CondorcetFlipBorda(BBox):
	return CondorcetPrefCount(BBox,CondorcetFlipBordaCount)

def CondorcetCmbnBordaCount(PrefMat):
	cnt1 = CondorcetBordaCount(PrefMat)
	cnt2 = CondorcetFlipBordaCount(PrefMat)
	return tuple(zip(cnt1,cnt2))

def CondorcetCmbnBorda(BBox):
	return CondorcetPrefCount(BBox,CondorcetCmbnBordaCount)

# Black's method
def CondorcetWithFallback(BBox, fallback=Borda, settype="All"):
	cond = CondorcetWinner(BBox)
	if cond[0]:
		return cond
	else:
		if settype == "All":
			fbblts = BBox
		else:
			fbblts = KeepCandidates(BBox, MaximalSet(BBox,"Smith"))
		return (False, fallback(fbblts))

def CondorcetSequentialRunoff(BBox, CWn=False, CLs=False, \
	ThresholdFunc=min, DoRound=TopOne):
	NewBBox = BBox
	rounds = []
	while len(NewBBox.Ballots) > 0:
		if CWn:
			condwin = CondorcetWinner(NewBBox)
			if condwin[0]:
				rounds.append((1,condwin[1]))
				break
		
		if CLs:
			condlose = CondorcetLoser(NewBBox)
			if condlose[0]:
				rounds.append((-1,condlose[1]))
				NewBBox = RemoveCandidates(NewBBox, [condlose[1]])
				continue
		
		res = DoRound(NewBBox)
		rounds.append((0,res))
		
		BottomCands = []
		VoteThreshold = ThresholdFunc([r[1] for r in res])
		for r in res:
			if r[1] <= VoteThreshold:
				BottomCands.append(r[0])	
		BottomCands.sort()
		BottomCands = tuple(BottomCands)
		
		NewBBox = RemoveCandidates(NewBBox, BottomCands)

	return tuple(rounds)


# http://en.wikipedia.org/wiki/Copeland%27s_method

def compare(x,y):
	if x < y: return -1
	elif x > y: return 1
	else: return 0

def CopelandPrefCount(PrefMat):
	n = len(PrefMat)
	
	Count = n*[0]
	for k1 in range(n):
		for k2 in range(n):
			Count[k1] += compare(PrefMat[k1][k2],PrefMat[k2][k1])
	
	return Count

def Copeland(BBox): return CondorcetPrefCount(BBox,CopelandPrefCount)


# http://en.wikipedia.org/wiki/Schulze_method

def OrderingMatrix(n, ordf):
	return tuple(( tuple(( ordf(i,j) for j in range(n) )) for i in range(n) ))

def ListOrderFromMatrix(ordmat):
	n = len(ordmat)
	ixskeys = []
	for i in range(n):
		keyval = 0
		for j in range(n):
			omdiff = ordmat[i][j] - ordmat[j][i]
			if omdiff > 0:
				keyval += 1
			elif omdiff < 0:
				keyval -= 1
		ixskeys.append( (i,keyval) )
	
	ixskeys.sort(key=lambda x: x[1])
	return tuple( (x[0] for x in ixskeys) )

def BeatpathMatrix(PrefMat):
	n = len(PrefMat)
	
	# Initial matrix
	BPMat = [n*[0] for k in range(n)]
	
	for i in range(n):
		for j in range(n):
			if j != i:
				if PrefMat[i][j] > PrefMat[j][i]:
					BPMat[i][j] = PrefMat[i][j]
			# Otherwise zero - BPMat initialized to that
	
	# Variant of the Floyd-Warshall algorithm
	for i in range(n):
		for j in range(n):
			if j != i:
				for k in range(n):
					if i != k and j != k:
						BPMat[j][k] = \
							max(BPMat[j][k], min(BPMat[j][i],BPMat[i][k]))
	
	return BPMat

def SchulzePrefCount(PrefMat):
	
	# Beatpaths
	BPMat = BeatpathMatrix(PrefMat)
	
	# Win-Lose Counts
	return CopelandPrefCount(BPMat)

def SchulzePrefOrder(PrefMat):
	n = len(PrefMat)

	counts = SchulzePrefCount(PrefMat)
	
	withnums = list(zip(range(n),counts))
	
	withnums.sort(key=lambda x: x[1], reverse=True)
	
	return tuple([x[0] for x in withnums])

def Schulze(BBox): return CondorcetPrefCount(BBox,SchulzePrefCount)


# Sorters for multiseat methods mentioned earlier

def SortWithSchulze(CandPrefMat):
	Cands, PrefMat = CandPrefMat
	CandsOrdered = CondorcetCandPMPrefOrder(Cands, PrefMat, SchulzePrefOrder)
	if len(CandsOrdered) > 0:
		return CandsOrdered[0]
	else:
		return None

def CPOSTV(BBox, Num): return SortWithSchulze(CPOSTVMatrix(BBox, Num))

def SchulzeSTV(BBox, Num): return SortWithSchulze(SchulzeSTVMatrix(BBox, Num))


# http://en.wikipedia.org/wiki/Minimax_condorcet

def ModifiedPrefMat(PrefMat,Which):
	n = len(PrefMat)
	ModMat = [n*[0] for i in range(n)]
	
	if Which == "wins":
		for i in range(n):
			for j in range(n):
				if PrefMat[i][j] > PrefMat[j][i]:
					ModMat[i][j] = PrefMat[i][j]
				else:
					ModMat[i][j] = 0
	elif Which == "marg":
		for i in range(n):
			for j in range(n):
				ModMat[i][j] = PrefMat[i][j] - PrefMat[j][i]
	elif Which == "oppo":
		for i in range(n):
			for j in range(n):
				ModMat[i][j] = PrefMat[i][j]
	
	return ModMat
		

def MinimaxPrefCount(PrefMat,Which):
	ModMat = ModifiedPrefMat(PrefMat,Which)

	n = len(PrefMat)
	Scores = n*[0]
	for k in range(n):
		Score = None
		for kx in range(n):
			if kx == k: continue
			NewScore = ModMat[kx][k]			
			if Score == None:
				Score = NewScore
			elif NewScore > Score:
				Score = NewScore
		
		Scores[k] = - Score if n > 1 else 0
	
	return Scores

def Minimax(BBox,Which):
	return CondorcetPrefCount(BBox,lambda PM: MinimaxPrefCount(PM,Which))


# http://en.wikipedia.org/wiki/Ranked_pairs
# Uses cyclicity test in
# http://www.cs.hmc.edu/~keller/courses/cs60/s98/examples/acyclic/

# Uses the choices for the modified Condorcet preference matrix
# that are in the minimax method

# The testers are objects with two methods:
# __init__
# __call__ on BeatList, BLCand
# BeatList: sequence type of tuples of pairs of candidates
# BLCand: tuple of a pair of candidates to test
# Returns whether or not the added candidate makes the graph acyclic

class RPTester_RemoveEnds:

	def __init__(self, NCands):
		self.n = NCands
		
	def __call__(self, BeatList, BLCand):
		# Nontrivial test: remove all the nodes that are
		# only on left sides or on right sides
		# Repeat until it is no longer possible to remove any more
		# Are any left?
		BLX = list(BeatList) + [BLCand]
		
		while True:
			# What is present on each side?
			NumLeft = self.n*[0]
			NumRight = self.n*[0]
			for bdx in BLX:
				NumLeft[bdx[0]] += 1
				NumRight[bdx[1]] += 1
			
			# What is present on only one side?
			RemoveLeft = []
			RemoveRight = []
			for k in range(self.n):
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
		return len(BLX) == 0


class RPTester_DepthFirst:
	def __init__(self,NCands):
		# Depth-first setup:
		# Use iteration and an explicit stack rather than recursion
		# In the search stack:
		#   Candidate
		#   Position in list of next candidates
		# List of next candidates for each candidate
		self.CandSeq = NCands*[0]
		self.CandSeqPos = NCands*[0]
		self.NextCand = [[] for i in range(NCands)]
		
		# Force initialization
		self.ix = -1
		
	def __call__(self, BeatList, BLCand):
		(bd0, bd1) = BLCand
		
		# Starting?
		if self.ix < 0:
			self.ix = 0
			self.CandSeq[self.ix] = 0
			self.CandSeqPos[self.ix] = 0
			self.NextCand[bd0].append(bd1)
			return True
		
		# Continuing
		
		# Initial state:
		# Use final candidate of pair
		self.ix = 0
		self.CandSeq[self.ix] = bd1
		self.CandSeqPos[self.ix] = 0
		
		# Add it unless we find that it causes a cycle				
		while True:
			# Check on whether one has found a cycle
			cd = self.CandSeq[self.ix]
			if cd == bd0:
				# Found a cycle
				return False
			
			# Where to go next?
			nxtcds = self.NextCand[cd]
			nix = self.CandSeqPos[self.ix]
			if nix < len(nxtcds):
				# Advance to next candidate
				self.ix += 1
				self.CandSeq[self.ix] = nxtcds[nix]
				self.CandSeqPos[self.ix] = 0
			else:
				# Back up
				# If not possible, then the search is finished
				if self.ix <= 0: break
				self.ix -= 1
				self.CandSeqPos[self.ix] += 1
		
		# The pair does not create a cycle,
		# so we go ahead with it
		self.NextCand[bd0].append(bd1)	
		return True


class RPTester_AllNext:
	def __init__(self,NCands):
		self.Next = [set() for i in range(NCands)]
		
	def __call__(self, BeatList, BLCand):
		(bd0, bd1) = BLCand
		
		# Found a cycle?
		nx1 = self.Next[bd1]
		if bd0 in nx1: return False
		
		# If not, then update the first one's next ones
		nx0 = self.Next[bd0]
		nx0 |= nx1
		nx0.add(bd1)
		return True

# Takes the Condorcet matrix, the comparison algorithm,
# and the tester for being acyclic.

def RankedPairsCompare(x):
	return (x[2],-x[3])

def RankedPairsPrefCount(PrefMat,Which,Tester):
	ModMat = ModifiedPrefMat(PrefMat,Which)
	
	# Find the beat margins and sort them
	n = len(PrefMat)
	BeatMargins = []
	for k1 in range(n):
		for k2 in range(n):
			if k2 != k1:
				bmg = (k1, k2, ModMat[k1][k2], ModMat[k2][k1])
				BeatMargins.append(bmg)
	BeatMargins.sort(key=RankedPairsCompare)
	
	# Find the beat list of what beats what
	BeatList = set()
	for bmg in BeatMargins:
		# Initial and final candidates of a pair
		(bd0, bd1) = bmg[:2]
		bd = (bd0, bd1)
		
		# Is the reverse direction present?
		rvbd = (bd1, bd0)
		if rvbd in BeatList: continue
		
		if Tester(BeatList,bd):
			BeatList.add(bd)
	
	# The signs for for getting the right sort order
	OrdMat = [n*[0] for i in range(n)]
	for bd0, bd1 in BeatList:
		OrdMat[bd0][bd1] = -1
		OrdMat[bd1][bd0] = 1
	
	return CopelandPrefCount(OrdMat)

# Global variable for which ranked-pairs tester to use
# RPTester = RPTester_RemoveEnds
# RPTester = RPTester_DepthFirst
RPTester = RPTester_AllNext

def RankedPairs(BBox, Which="oppo"):
	NCands = len(BBox.Candidates())
	Tester = RPTester(NCands)
	return CondorcetPrefCount(BBox, \
		lambda BB: RankedPairsPrefCount(BB,Which,Tester))


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
		for k in range(n-1):
			if perm[k] < perm[k+1]:
				ix0 = k
		if ix0 == None: break
		
		ix1 = ix0+1
		for k in range(ix0+1,n):
			if perm[ix0] < perm[k]:
				ix1 = k
		
		temp = perm[ix0]
		perm[ix0] = perm[ix1]
		perm[ix1] = temp
		
		for k in range(ix0+1,n):
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
		for k1 in range(n-1):
			for k2 in range(k1+1,n):
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
	
	for i in range(n):
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

def CDLSortFunc(x):
	return (x[1],x[0])

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
	CDList.sort(key=CDLSortFunc)
	return tuple(CDList)


# https://en.wikipedia.org/wiki/Maximal_lotteries

# Use linear programming

# Solve a system of linear equations with the Gauss-Jordan algorithm
# The input has form (matrix) (vector, vector, ...)
def SolveLinearEquations(mat):
	n = len(mat)
	nx = len(mat[0])
	workmat = [[float(mat[i][j]) for j in range(nx)] for i in range(n)]
	
	# Do forward substitution
	for icol in range(n):
		# Necessary to exchange rows
		# to bring a nonzero value into position?
		# Return None if singular
		if workmat[icol][icol] == 0:
			ipvt = None
			for i in range(icol+1,n):
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
		for i in range(icol,nx):
			wmicol[i] *= dgvrecip
		# Forward substitute:
		for i in range(icol+1,n):
			wmi = workmat[i]
			elimval = wmi[icol]
			for j in range(icol,nx):
				wmi[j] -= elimval*wmicol[j]
	
	# Do back substitution
	for icol in range(n-1,0,-1):
		wmicol = workmat[icol]
		for i in range(icol):
			wmi = workmat[i]
			elimval = wmi[icol]
			for j in range(icol,nx):
				wmi[j] -= elimval*wmicol[j]
	
	# Done!
	return [[workmat[i][j] for j in range(n,nx)] for i in range(n)]

def LP_SimplexMethodStep(basvars, tableau):
	nrows = len(tableau)
	ncols = len(tableau[0])
	
	# Find pivot column
	pvcol = -1
	objf = tableau[0]
	for k in range(ncols-1):
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
	for k in range(1,nrows):
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
	for k in range(ncols):
		xrow[k] /= xrc
	for j in range(nrows):
		if j == pvrow: continue
		yrow = tableau[j]
		yrc = yrow[pvcol]
		for k in range(ncols):
			yrow[k] -= yrc*xrow[k]
	
	# Fix the numerical values
	for k in range(nrows):
		tableau[j][pvcol] = 1. if k == pvrow else 0.
	
	# Make approximate cancellations exact
	objf = tableau[0]
	objfmax = max(max(objf),-min(objf))
	eps = (1e-12)*objfmax
	for k in range(ncols-2):
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
	tableau = [[float(tableau_[i][j]) for j in range(ncols)] for i in range(nrows)]
	
	basvars = nrows*[-1]
	basvars[0] = ncols-2
	allzero = False
	for i in range(ncols-2):
		numnz = 0
		for j in range(1,nrows):
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
	nzvrixs = [k for k in range(ncols-2) if objf[k] == 0]
	nzvrdiff = list(set(nzvrixs) - set(basvars))
	nzvrdiff.sort()
	nzvrixs = basvars + nzvrdiff
	if len(nzvrixs) < (nrows-1): return (False, "Unsolvable")
	
	# Handle degenerate cases by using the first of the variables
	nzvrixs = nzvrixs[:nrows-1]
	
	# Solve a system of linear equations by
	# Gaussian elimination with pivoting
	slmat = [nrows*[0] for k in range(nrows-1)]
	for j in range(nrows-1):
		slmrow = slmat[j]
		tblrow = tableau[j+1]
		for k in range(nrows-1):
			slmrow[k] = tblrow[nzvrixs[k]]
	for k in range(nrows-1):
		slmat[k][-1] = tableau[k+1][-1]
	res = SolveLinearEquations(slmat)
	
	vals = (ncols-1)*[0]
	for k in range(nrows-1):
		vals[nzvrixs[k]] = res[k][0]
	
	# Solve for the final slack variable
	vsum = objf[ncols-1]
	for k in range(ncols-2):
		vsum -= objf[k]*vals[k]
	vals[ncols-2] = vsum/objf[ncols-2]
	
	return (True, vals)


def MaximalLotteriesPrefCount_Simplex(PrefMat):
	n = len(PrefMat)
	
	if n == 0: return ()
	
	tableau = [(2*n+3)*[0] for k in range(n+2)]
	
	for i in range(n):
		for j in range(n):
			tableau[i+2][j] = PrefMat[i][j] - PrefMat[j][i]
	
	for i in range(n):
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
	for i in range(n):
		mr = mat[i]
		s = 0
		for j in range(n):
			s += mr[j]*vec[j]
		res[i] = s

	return res

def LineCalc(pt,dst,dir):
	n = len(pt)
	return [pt[i] + dst*dir[i] for i in range(n)]

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
	PMAS = [[PrefMat[j][i] - PrefMat[i][j] for j in range(n)] for i in range(n)]
	
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
	for iterrpt in range(mlitrrpt):

		# Start a point
		NewPosScores = Scores
		NewPosCstrMin = CstrMin
		
		# Iterate over advancement steps
		for iterpos in range(mlitrpos):
			
			# Start a direction search
			NewDirScores = NewPosScores
			NewDirCstrMin = NewPosCstrMin
		
			# Iterate over directions
			for iterdir in range(mlitrdir):
				
				# Find a random direction
				# Be sure that no component is zero
				Dir = [random.gauss(0,1) for i in range(n)]
				DirAvg = sum(Dir)/n
				Dir = [Dir[i] - DirAvg for i in range(n)]
				for d in Dir:
					if d == 0: continue
				
				# Find the interval that is allowable
				# by the scores being nonnegative and adding up to 1
				PMin = None
				PMax = None
				for i in range(n):
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
				for i in range(0,LSDivs+1):
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
# From https://electowiki.org/wiki/Maximal_elements_algorithms

def MaximalSetIndices(PrefMat, Type):
	n = len(PrefMat)
	
	# Init HasPath for relations with length 1
	HasPath = [n*[False] for k in range(n)]
	for k1 in range(n):
		for k2 in range(n):
			if k2 != k1:
				if Type == "Schwartz":
					HasPath[k1][k2] = PrefMat[k1][k2] > PrefMat[k2][k1]
				elif Type == "Smith":
					HasPath[k1][k2] = PrefMat[k1][k2] >= PrefMat[k2][k1]
	
	# Consider paths with intermediate nodes from 1 to k
	for k in range(n):
		for k1 in range(n):
			if k1 != k:
				for k2 in range(n):
					if k2 != k and k2 != k1:
						if HasPath[k1][k] and HasPath[k][k2]:
							HasPath[k1][k2] = True
	
	# Candidates with paths to them but none to complete a cycle,
	# they are not in the maximal set
	InMaximal = n*[True]
	for k1 in range(n):
		for k2 in range(n):
			if k2 != k1:
				if HasPath[k2][k1] and not HasPath[k1][k2]:
					InMaximal[k1] = False
	
	# Find indices
	return tuple((im[0] for im in zip(range(n),InMaximal) if im[1]))

def MaximalSet(BBox, Type):
	Cands = BBox.Candidates()
	PrefMat = BBox.CondorcetMatrix()
	
	MSIxs = MaximalSetIndices(PrefMat, Type)
	
	return tuple((Cands[ix] for ix in MSIxs))

def MaximalSetSequenceIndices(PrefMat, Type):
	MSIxSeq = []
	
	# BaseSeq is for containing the original indices,
	# because MaximalSetIndices finds them relative to
	# its input's size
	n = len(PrefMat)
	BaseSeq = list(range(n))
	PrevPrefMat = PrefMat
	
	while True:
		# Get the set's indices and translate them into the original ones
		MSIxs = MaximalSetIndices(PrevPrefMat, Type)
		MSIxSeq.append(tuple((BaseSeq[i] for i in MSIxs)))
		
		# Find remaining indices and trim down
		# the base indices and the preference matrix with them
		OrigIxs = set(range(len(PrevPrefMat)))
		SubIxs = set(MSIxs)
		RmIxsSet = OrigIxs - SubIxs
		RmIxs = list(RmIxsSet)
		RmIxs.sort()
		BaseSeq = [BaseSeq[i] for i in RmIxs]
		NewPrefMat = [[PrevPrefMat[i][j] for j in RmIxs] for i in RmIxs]
		# Any more left to do?
		if len(NewPrefMat) == 0: break
		PrevPrefMat = NewPrefMat
	
	return tuple(MSIxSeq)

def MaximalSetSequence(BBox, Type):
	
	Cands = BBox.Candidates()
	MSIxSeq = MaximalSetSequenceIndices(BBox.CondorcetMatrix(), Type)
	return tuple((tuple((Cands[ix] for ix in MSIxs)) for MSIxs in MSIxSeq))

# Tideman alternative method
# https://en.wikipedia.org/wiki/Tideman_alternative _method
#
def MaximalSetSequentialRunoff(BBox, Type, ThresholdFunc=min, DoRound=TopOne):
	NewBBox = BBox
	rounds = []
	while len(NewBBox.Ballots) > 0:
		MaxSet = MaximalSet(NewBBox, Type)
		if len(MaxSet) == 0: break
		NewBBox = KeepCandidates(NewBBox, MaxSet)
		
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

# All the subsets of a Container
# Returns as tuple of tuples
def Subsets(ctnr):
	ssts = [[]]
	for x in ctnr:
		nwssts = []
		for st in ssts:
			nwssts.append(st)
			nwssts.append(st+[x])
		ssts = nwssts
	return tuple(map(tuple,ssts))
	
# Descending Solid Coalitions
# From https://electowiki.org/wiki/Maximal_elements_algorithms

def DescendingSolidCoalitions(BBox):
	sbsts = Subsets(BBox.Candidates())
	# Nonempty ones only
	sbsts = [s for s in sbsts if len(s) > 0]
	sbcts = []
	for sb in sbsts:
		sbln = len(sb)
		sbct = 0
		for Ballot in BBox.Ballots:
			Weight = Ballot[0]
			Votes = Ballot[1]
			if len(Votes) < sbln: continue
			vtc = list(Votes[:sbln])
			vtc.sort()
			if tuple(vtc) == sb:
				sbct += Weight
		sbcts.append(sbct)
	sbcx = list(zip(sbsts,sbcts))
	sbcx.sort(key=lambda x: -x[1])
	
	candset = set(BBox.Candidates())
	winset = set(BBox.Candidates())
	for sb,ct in sbcx:
		sbc = candset - set(sb)
		nwwinset = winset - sbc
		if len(nwwinset) > 0:
			winset = nwwinset
	winlst = list(winset)
	winlst.sort()
	return tuple(winlst)

#
# For debugging
#

def DumpSingleAlgorithm(Algorithm,BBox):

	# For reference
	print("Schulze and Borda:")
	print(Schulze(BBox))
	print(Borda(BBox))
	print()
	
	print("Algorithm to Test:")
	res = Algorithm(BBox)
	print(res)
	print()


def DumpMultiWinnerAlgorithms(BBox):

	Cands = BBox.Candidates()
	NumCands = len(Cands)

	for k in range(min(3,NumCands)):
		Num = k+1
		print("Single Transferable Vote:", Num)
		reslist = SingleTransferableVote(BBox,Num)
		# Only show the final result
		print(reslist[-1])
		print()

	for k in range(min(3,NumCands)):
		Num = k+1
		print("CPO by STV:", Num)
		print(CPOSTV(BBox,Num))
		print()

	for k in range(min(3,NumCands)):
		Num = k+1
		print("Schulze STV:", Num)
		print(SchulzeSTV(BBox,Num))
		print()

# import time

def DumpAll(Ballots, DoNFactorial=True):
	BBox = BallotBox(Ballots)
	
	DumpSingleAlgorithm(lambda b: RankedPairs(b),BBox); return
	# DumpMultiWinnerAlgorithms(BBox); return
	
	print("Candidates:",)
	Cands = BBox.Candidates()
	NumCands = len(Cands)
	for Cand in Cands: print(Cand, " ",)
	print()
	print()
	
	print("Total Count:", TotalCount(BBox))
	print()
	
	print("All Votes:")
	res = AllVotes(BBox)
	for r in res: print(r)
	print()
	
	print("Top One (First Past the Post, Plurality):")
	res = TopOne(BBox)
	for r in res: print(r)
	print()
	
	print("All Top N's (Bucklin):")
	reslist = TopNumList(BBox)
	for k, res in enumerate(reslist):
		print("N =", k+1)
		for r in res: print(r)
	print()
	
	print("Borda Count:")
	res = Borda(BBox)
	for r in res: print(r)
	print()
	
	print("Modified Borda Count:")
	res = ModBorda(BBox)
	for r in res: print(r)
	print()
	
	print("Cumulative Borda Count:")
	res = CumulBorda(BBox)
	for r in res: print(r)
	print()
	
	print("Dowdall System:")
	res = Dowdall(BBox)
	for r in res: print(r)
	print()
	
	print("Majority Judgment:")
	res = MajorityJudgment(BBox)
	for r in res: print(r)
	print()
	
	print("Top-Two Runoff:")
	reslist = TopTwoRunoff(BBox)
	for k, res in enumerate(reslist):
		print("Round", k+1)
		for r in res: print(r)
	print()
	
	print("STAR Voting:")
	reslist = TopTwoRunoff(BBox,Borda)
	for k, res in enumerate(reslist):
		print("Round", k+1)
		for r in res: print(r)
	print()

	print("Sequential Runoff:")
	reslist = SequentialRunoff(BBox)
	for k, res in enumerate(reslist):
		print("Round", k+1)
		for r in res: print(r)
	print()

	print("Winner Dropping for SR:")
	reslistlist = WinnerDropping(BBox, SequentialRunoff, lambda r: r[-1][0][0])
	for k, reslist in enumerate(reslistlist):
		print("Overall Round", k+1)
		for l, res in enumerate(reslist):
			print("SR Round", l+1)
			for r in res: print(r)
	print()
	
	print("Baldwin:")
	reslist = SequentialRunoff(BBox,min,Borda)
	for k, res in enumerate(reslist):
		print("Round", k+1)
		for r in res: print(r)
	print()
	
	print("Nanson:")
	reslist = SequentialRunoff(BBox,Average,Borda)
	for k, res in enumerate(reslist):
		print("Round", k+1)
		for r in res: print(r)
	print()
	
	print("Coombs:")
	reslist = SequentialRunoff(BBox,min,MinusBottom)
	for k, res in enumerate(reslist):
		print("Round", k+1)
		for r in res: print(r)
	print()

	for k in range(min(3,NumCands)):
		Num = k+1
		print("Single Transferable Vote:", Num)
		reslist = SingleTransferableVote(BBox,Num)
		for k, res in enumerate(reslist):
			print("Round", k+1)
			for r in res: print(r)
		print()

	for k in range(min(3,NumCands)):
		Num = k+1
		print("CPO by STV:", Num)
		print(CPOSTV(BBox,Num))
		print()

	for k in range(min(3,NumCands)):
		Num = k+1
		print("Schulze STV:", Num)
		print(SchulzeSTV(BBox,Num))
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
	
	print("Condorcet Loser:")
	print(CondorcetLoser(BBox))
	print()
	
	print("Condorcet Combined Borda (plain + flipped):")
	res = CondorcetCmbnBorda(BBox)
	for r in res: print(r)
	print()
	
	print("Black: Condorcet with Fallback")
	res = CondorcetWithFallback(BBox)
	print(res[0])
	if res[0]:
		print(res[1])
	else:
		for r in res[1]: print(r)
	print()
	
	print("Black: Condorcet with Fallback (Sequential Runoff, Smith set)")
	res = CondorcetWithFallback(BBox,SequentialRunoff,"Smith")
	print(res[0])
	if res[0]:
		print(res[1])
	else:
		for k, rx in enumerate(res[1]):
			print("Round", k+1)
			for r in rx: print(r)
	print()
	
	print("Condorcet SR using C Winner:")
	reslist = CondorcetSequentialRunoff(BBox,True)
	for k, res in enumerate(reslist):
		print("Round", k+1)
		for r in res: print(r)
	print()
	
	print("Condorcet SR dropping C Loser:")
	reslist = CondorcetSequentialRunoff(BBox,False,True)
	for k, res in enumerate(reslist):
		print("Round", k+1)
		for r in res: print(r)
	print()
	
	print("Copeland Pairwise Aggregation:")
	res = Copeland(BBox)
	for r in res: print(r)
	print()
		
	print("Schulze Beatpath:")
	res = Schulze(BBox)
	for r in res: print(r)
	print()
	
	for Which in ("wins", "marg", "oppo"):
		print("Minimax:", Which)
		res = Minimax(BBox,Which)
		for r in res: print(r)
		print()
	
	print("Tideman Ranked Pairs:")
	res = RankedPairs(BBox)
	for r in res: print(r)
	print()
	
	print("Maximize Affirmed Majorities:")
	res = RankedPairs(BBox,"wins")
	for r in res: print(r)
	print()
	
	print("Kemeny-Young:")
	if DoNFactorial:
		print(KemenyYoung(BBox))
	else:
		print("Skipped because it is O(n!)")
	print()
	
	print("Dodgson:")
	if DoNFactorial:
		res = Dodgson(BBox)
		for r in res: print(r)
	else:
		print("Skipped because it is O(n!)")
	print()
	
	print("Maximal Lotteries:")
	res = MaximalLotteries(BBox)
	for r in res: print(r)
	print()
	
	print("Maximal sets: Smith and Schwartz sets:")
	MaxSetSmith = MaximalSet(BBox,"Smith")
	MaxSetSchwartz = MaximalSet(BBox,"Schwartz")
	BBoxSmith = KeepCandidates(BBox, MaxSetSmith)
	BBoxSchwartz = KeepCandidates(BBox, MaxSetSchwartz)
	print(MaxSetSmith)
	print(MaxSetSchwartz)
	print()
	
	print("Borda on Smith Set:")
	res = Borda(BBoxSmith)
	for r in res: print(r)
	print()
	
	print("Borda on Schwartz Set:")
	res = Borda(BBoxSchwartz)
	for r in res: print(r)
	print()
	
	print("IRV on Smith Set:")
	reslist = SequentialRunoff(BBoxSmith)
	for k, res in enumerate(reslist):
		print("Round", k+1)
		for r in res: print(r)
	print()
	
	print("IRV on Schwartz Set:")
	reslist = SequentialRunoff(BBoxSchwartz)
	for k, res in enumerate(reslist):
		print("Round", k+1)
		for r in res: print(r)
	print()
	
	print("Sequence of maximal sets: Smith and Schwartz sets:")
	print(MaximalSetSequence(BBox,"Smith"))
	print(MaximalSetSequence(BBox,"Schwartz"))
	print()

	print("Smith-Set Sequential Runoff:")
	reslist = MaximalSetSequentialRunoff(BBox,"Smith")
	for k, res in enumerate(reslist):
		print("Round", k+1)
		for r in res: print(r)
	print()
		
	print("Descending Solid Coalitions")
	res = DescendingSolidCoalitions(BBox)
	for r in res: print(r)
	print()


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
	
	print("*** Empty Ballot Box ***")
	DumpAll(ZeroBallots)
	
	print("*** One-Party Election ***")
	DumpAll(OneBallot)
	
	print("*** Two-Party Election ***")
	DumpAll(TwoBallots)
	
	print("*** Beers ***")
	DumpAll(BeerBallots)
	
	print("*** Distilled Spirits ***")
	DumpAll(DistilledBallots)
	
	print("*** Tennessee Capital ***")
	DumpAll(TennesseeCapitalBallots)
	
	print("*** ABCDS ***")
	DumpAll(ABCDSBallots)
	
	print("*** Ice-Cream Flavors ***")
	DumpAll(MakeWeighted(ICFBallots))