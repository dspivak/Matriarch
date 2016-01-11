# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 10:21:53 2015

@author: Tristan Giesa, David I. Spivak, Ravi Jagadeesan
"""

from matriarch import *

mySeqG = 'GGGGGGGGG'
mySeqA = 'AAA'
myChainG = chain(mySeqG)
myChainA = chain(mySeqA)
spacedBlock = space(myChainG, myChainA, 10)
array = makeArray(spacedBlock, 10, 15, 4, 6, True)
fileOut(array, '11-arraySpaced.pdb')
