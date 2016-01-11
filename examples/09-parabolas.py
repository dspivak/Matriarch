# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 10:21:53 2015

@author: Tristan Giesa, David I. Spivak, Ravi Jagadeesan
"""

from matriarch import *

def parabola(t):
    return [-t*t,0,10*t]

mySeq = 'AAAAAAAAAAAPPPPPPPPPPP'
myChain = chain(mySeq)
Rout = Ray([-1,0,0], [0,0,0])
W = buildAxisTwister(parabola,Rout)
pos = twist(myChain,W)
rev = reverseOrbs(myChain)
twi = twist(rev,W)
neg = reverseOrbs(twi)
t = attach(neg,pos)
fileOut(t, '09-parab.pdb')
fileOut(pos, '09-parab-pos.pdb')
fileOut(rev, '09-parab-rev.pdb')
fileOut(twi, '09-parab-twi.pdb')
fileOut(neg, '09-parab-neg.pdb')