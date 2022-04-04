# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 10:21:53 2015

@author: Tristan Giesa, David I. Spivak, Ravi Jagadeesan
"""

from matriarch import *
from math import *

seq1 = "AAAAGGPGGYGGPGGAAAA"
a = chain(seq1)
ser = attachSeries(a, 25)


def SphSprl(k, Rout):
    def curve(t):
        return [k * sin(t) * cos(20 * t), k * sin(t) * sin(20 * t), k * cos(t)]

    tmax = pi
    Thetaspec = 0
    return buildAxisTwister(curve, Rout, Thetaspec, tmax)


Rout = Ray([0, 0, 0], [0, 0, 0])
k1 = 1  # first try
W1 = SphSprl(k1, Rout)

contourLength = length(ser)
lengthOfCurveWithKEquals1 = W1[0].length
knew = contourLength / lengthOfCurveWithKEquals1

SphSprlTwister = SphSprl(knew + 0.001, Rout)
output = twist(ser, SphSprlTwister)

fileOut(output, "13-sphericalSpiral.pdb")
