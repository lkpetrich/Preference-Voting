#!/usr/bin/env python
#
# Example from
# http://freethoughtblogs.com/dispatches/2013/11/14/time-to-rewrite-the-constitution/
#
# jd142
# 47
# Alice loves pizza, likes dim sum, doesn't mind tapas and hates hamburgers
# Bob loves hamburgers, likes dim sum, doesn't mind tapas, and hates pizza
# Charlie loves hamburgers, likes tapas, doesn't mind dim sum and hates pizza
# Dawn loves dim sum, likes pizza, doesn't mind tapas and hates hamburgers
# Eve loves tapas, likes dim sum, doesn't mind pizza and hates hamburger

from PrefVote import *

bsrc = []

# Alice
b = ("pizza", "dim sum", "tapas", "hamburgers")
bsrc.append(b)

# Bob
b = ("hamburgers", "dim sum", "tapas", "pizza")
bsrc.append(b)

# Charlie
b = ("hamburgers", "tapas", "dim sum", "pizza")
bsrc.append(b)

# Dawn
b = ("dim sum", "pizza", "tapas", "hamburgers")
bsrc.append(b)

# Eve
b = ("tapas", "dim sum", "pizza", "hamburgers")
bsrc.append(b)

DumpAll(MakeWeighted(bsrc))
