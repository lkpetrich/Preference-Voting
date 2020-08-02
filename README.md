# Preference-Voting
For counting votes in preference or ranked-choice voting. Large number of algorithms implemented.

PrefVote.py -- in Python. Instructions on how to use it are in the file itself.
LunchMeals.py -- in Python. A simple example of preference voting.
Preference-Voting Algorithms.nb -- in Mathematica.

It works from lists of weighted ballots, with format (weight, (list of candidates in preference order)),
but it has a convenience function, MakeWeighted, for adding a weight value to unweighted ballots (default 1).

It implements several vote-counting methods, in several categories of methods.

All-or-nothing votes:
- Top One: first-past-the-post or plurality voting
- Top N: approval voting with N votes
- Minus Bottom: lowest-ranked and unranked (used in some other methods)

Rated or scored votes:
- Borda: rating starts at # candidates, goes down with rank
- Modified Borda: rating starts at # voted-for candidates, goes down with rank
- Cumulative Borda: simulates cumulative voting by dividing modified Borda by # voted-for candidates
- Dowdall: rating is 1, 1/2, 1/3, ...
- Majority Judgment: sorts by weighted medians of ranks

Multiple-round methods:
- Winner Dropping: finds a winner with some method, drops that candidate, then repeats
- Top-Two Runoff: does an initial count, then one with the top two winners. Default: top one, STAR: Borda
- Sequential Runoff: counts, drops lowest scorer(s), repeats
- Variations like Nanson's, Bucklin's, and Coombs's methods
- Single Transferable Vote (STV): like sequential runoff, but drops winners and reweights their ballots
- CPO-STV: Comparison of Pairs of Outcomes by Single Transferable Vote
- Schulze STV: Like above, but only comparing outcomes that differ by one candidate

Condorcet methods:
- Calculate the Condorcet matrix: preferences as one-on-one contests -- virtual round robin
- Condorcet winner and loser, if present
- Condorcet-Borda: variant of the Borda count implemented with the Condorcet matrix. Has winner, loser, and combined scores
- Condorcet-IRV; optionally search for Condorcet winner and/or loser in sequential runoff
- Condorcet with Fallback (Black's method): uses the Condorcet winner, or else some other method. Can use the Smith or Schwartz sets
- Schulze beatpath
- Copeland: reduces Condorcet matrix to winner, tie, loser entries, then adds them up
- Minimax: finds candidate with the smallest of its maximum defeats
- Kemeny-Young: goes through all permutations of candidate rankings, finds best-performance one
- Dodgson: goes through all permutations of ballot orders, finds closest one to identity with a Condorcet winner
- Tideman ranked pairs: assembles the best-performing ranking of candidates from the highest-ranked pairs
- Maximize Affirmed Majorities: ranked pairs with different sorting
- Maximal lotteries: generalizes Condorcet winner with the help of linear programming

Maximal sets:
- Smith: smallest set of candidates that beat all others
- Schwartz: union of sets of candidates that beat or tie all others, with no proper subset with that property

Extra Methods:
- Descending Solid Coalitions
