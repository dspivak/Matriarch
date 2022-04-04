# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 10:21:53 2015

@author: Tristan Giesa, David I. Spivak, Ravi Jagadeesan
"""

from matriarch import *
from math import *


def triangleparameter(side, pitch, length, smoothingFactor):
    iterations = int(length / (3 * side)) + 2
    const = sqrt(3) / 6

    def loop(Z0):
        return [
            [side * 2 * const, 0, Z0],
            [-side * const, side / 2.0, Z0 + pitch / 3.0],
            [-side * const, -side / 2.0, Z0 + pitch * 2 / 3.0],
        ]

    PList = []
    for n in range(0, iterations):
        PList.extend(loop(pitch * n))
    Rout = Ray([0, 0, 0], [0, 0, 1])
    provideTheta = []
    return smoothedPieceWiseLinear(PList, Rout, provideTheta, smoothingFactor)


def triangle(a, side, pitch, length, smoothingFactor):
    return twist(a, triangleparameter(side, pitch, length, smoothingFactor))


aminoLength = 3.4

actualSeq = "TNVIIEGNVTLGHRVKIGTGCVIKNSVIGDDCEISP"
Mult = 3
Side = 9 * aminoLength
Pitch = 7
SmthFact = 0.33
totalLen = Mult * aminoLength * len(actualSeq)
myChain = attachSeries(chain(actualSeq), Mult)
myTriangle = triangle(myChain, Side, Pitch, totalLen, SmthFact)

fileOut(myTriangle, "12-triangleHelix.pdb")
