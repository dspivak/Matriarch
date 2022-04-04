# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 10:21:53 2015

@author: Tristan Giesa, David I. Spivak, Ravi Jagadeesan
"""

from matriarch import *
from math import *

mySeqA = "AAAAAA"
myChainA = chain(mySeqA)
myPadded = pad(myChainA, 3)
fileOut(myPadded, "08-myPadded.pdb")
