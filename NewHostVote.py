#!python3
#
# Does preference voting on web hosts, working from various reviews of them

from PrefVote import *

SourceNames = []
SourceURLs = []
Ballots = []

# Shared Hosting Reviews - http://webhostingreviews.us/shared/
b = ("InMotion", "WebHostingHub", "iPage", "FatCow", "JustHost")
SourceNames.append("Shared Hosting Reviews")
SourceURLs.append("http://webhostingreviews.us/shared/")
Ballots.append(b)

# Web Hosting Reviews - Hosting Reviews - Hosting Review - Best Web Hosting - http://www.hostingreview.com/
b = ("HostGator", "iPage", "BlueHost", "FatCow", "HostMonster")
SourceNames.append("Web Hosting Reviews - Hosting Reviews - Hosting Review - Best Web Hosting")
SourceURLs.append("http://www.hostingreview.com/")
Ballots.append(b)

# Web Hosting Reviews - http://www.webhostingreviewed.net/
SourceNames.append("Web Hosting Reviews")
SourceURLs.append("http://www.webhostingreviewed.net/")
b = ("BlueHost", "HostGator", "iPage")
Ballots.append(b)

# Web Hosting - Reviews of Domain Hosting Services - http://www.web-hosting-reviews.org/
b = ("BlueHost", "HostMonster", "StartLogic", "LunarPages", "JumpLine", "IPowerWeb", \
	"PowWeb", "HostRocket", "Website Source", "Apollo Hosting")
SourceNames.append("Web Hosting - Reviews of Domain Hosting Services")
SourceURLs.append("http://www.web-hosting-reviews.org/")
Ballots.append(b)

# Web Hosting Reviews - Thousands of web hosting USER reviews. - http://www.webhostingreviews.com/
b = ("InMotion", "HostGator", "WebzPro", "AwardSpace", "WebHostingHub")
SourceNames.append("Web Hosting Reviews - Thousands of web hosting USER reviews.")
SourceURLs.append("http://www.webhostingreviews.com/")
Ballots.append(b)

# Web Hosting Review 2012 | Best Web Hosting | Cheap Web Hosting - TopTenREVIEWS - http://web-hosting-review.toptenreviews.com/
b = ("JustHost", "GoDaddy", "iPage", "Arvixe", "MidPhase", "HostGator", "AN Hosting", \
	"Network Solutions", "Global", "WestHost", "FatCow", "BlueHost", "GreenGeeks", "LunarPages", \
	"myHosting", "HostWay", "DotEasy", "HostRocket", "DreamHost", "HostMonster", "Yahoo", \
	"1and1", "CoolHandle", "WebHostingPad", "Register.com", "JumpLine", "StartLogic", "InMotion", \
	"NetFirms", "FortuneCity", "Web.com", "Blue Domino", "PowWeb", "WebHero", "Gate.com", \
	"iPOWER", "Easy CGI", "Aplus")
SourceNames.append("Web Hosting Review 2012 | Best Web Hosting | Cheap Web Hosting - TopTenREVIEWS")
SourceURLs.append("http://web-hosting-review.toptenreviews.com/")
Ballots.append(b)

# Web Hosting Reviews You Can Trust - http://www.whoishostingthis.com/hosting-reviews/?sort=rating
b = ("InMotion", "BlueHost", "HostGator")
SourceNames.append("Web Hosting Reviews You Can Trust")
SourceURLs.append("http://www.whoishostingthis.com/hosting-reviews/?sort=rating")
Ballots.append(b)

# Best Web Hosting Services - Top 10 Website Hosts 2012 - http://webhostinggeeks.com/
b = ("InMotion", "WebHostingHub", "WebHostingPad", "FatCow", "iPage", "GreenGeeks", "JustHost", \
	"HostGator", "BlueHost", "HostMonster")
SourceNames.append("Best Web Hosting Services - Top 10 Website Hosts 2012")
SourceURLs.append("http://webhostinggeeks.com/")
Ballots.append(b)


# Remove all the hosts with too few counts
Counts = {}
for b in Ballots:
	for h in b:
		if h not in Counts: Counts[h] = 0
		Counts[h] += 1

AcceptedHosts = tuple([h for h in Counts if Counts[h] >= 3])

FrequentHostBallots = []
for b in Ballots:
	bnew = []
	for h in b:
		if h in AcceptedHosts: bnew.append(h)
	FrequentHostBallots.append(bnew)


BBoxAll = BallotBox(MakeWeighted(Ballots))
BBoxFrq = BallotBox(MakeWeighted(FrequentHostBallots))

def BrdTrim(x): return x[0]
def SRWinner(x): return x[-1][0][0]

def Eval(BBox):
	print()
	r1 = map(SRWinner, WinnerDropping(BBox, SequentialRunoff, SRWinner))
	r2 = map(BrdTrim,Borda(BBox))
	r3 = Schulze(BBox)
	res = zip(r1,r2, r3)
	for r in res: print(r)

print("('IRV Seq','Borda','Schulze')")
Eval(BBoxAll)
Eval(BBoxFrq)