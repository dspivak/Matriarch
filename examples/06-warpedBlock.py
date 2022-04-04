# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 10:21:53 2015

@author: Tristan Giesa, David I. Spivak, Ravi Jagadeesan
"""

import matriarch

mySeq = "AAAPPY"
myChain = matriarch.chain(mySeq)


def F(x):
    return [-x[2], x[1], x[0] - 25]


def FPrime(x):
    return [-x[2], x[1], x[0]]


g = [F, FPrime]
myChain_newAxis = matriarch.moveOrbs(myChain, g)


def parabola(t):
    return [-t * t, 0, 10 * t]


Rout = matriarch.Ray([18, 0, 0], [0, 0, 0])
W = matriarch.buildAxisTwister(parabola, Rout)
warpedBlock = matriarch.twist(myChain_newAxis, W)
matriarch.fileOut(warpedBlock, "06-warpedBlock.pdb")
