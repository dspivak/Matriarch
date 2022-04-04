# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 10:21:53 2015

@author: Tristan Giesa, David I. Spivak, Ravi Jagadeesan
"""

import matriarch

mySeqA = "AAAAAAAAAAAAAAAAAAAA"
myChainA = matriarch.chain(mySeqA)
myHelix = matriarch.helix(myChainA, 4, 8, "L")
myHelix.fileOut("02-myHelix.pdb")
