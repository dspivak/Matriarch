# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 10:21:53 2015

@author: Tristan Giesa, David I. Spivak, Ravi Jagadeesan
"""

import matriarch
mySeq = 'AAAPPY'
myChain = matriarch.chain(mySeq)
myLongChain = matriarch.attachSeries(myChain,5)
myLongChain.fileOut('01-myLongChain.pdb')