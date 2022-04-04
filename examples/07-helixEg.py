# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 10:21:53 2015

@author: Tristan Giesa, David I. Spivak, Ravi Jagadeesan
"""

import matriarch
from math import *

mySeqA = "AAAAAAAAAAAAAAAAAAAA"
myChainA = matriarch.chain(mySeqA)


def parameterizedHelix(t):
    return [4 * cos(2 * pi * t), -4 * sin(2 * pi * t), 8 * t]


W = matriarch.buildAxisTwister(parameterizedHelix)
myHelix = matriarch.twist(myChainA, W)
matriarch.fileOut(myHelix, "07-helixEg.pdb")
