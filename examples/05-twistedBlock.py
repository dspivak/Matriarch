# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 10:21:53 2015

@author: Tristan Giesa, David I. Spivak, Ravi Jagadeesan
"""

import matriarch

mySeq = "AAAPPY"
myChain = matriarch.chain(mySeq)


def parabola(t):
    return [-t * t, 0, 10 * t]


Rout = matriarch.Ray([18, 0, 0], [0, 0, 0])
W = matriarch.buildAxisTwister(parabola, Rout)
twistedBlock = matriarch.twist(myChain, W)
matriarch.fileOut(twistedBlock, "05-twistedBlock.pdb")
