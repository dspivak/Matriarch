# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 10:21:53 2015

@author: Tristan Giesa, David I. Spivak, Ravi Jagadeesan
"""

import matriarch

mySeq = "AAAPPY"
myChain = matriarch.chain(mySeq)


def T(x):
    return [-x[2], x[1], x[0] - 25]


def TPrime(x):
    return [-x[2], x[1], x[0]]


g = [T, TPrime]
movedChain = matriarch.moveOrbs(myChain, g)
movedChain.fileOut("04-movedChain.pdb")
