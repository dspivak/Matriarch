# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 10:21:53 2015

@author: Tristan Giesa, David I. Spivak, Ravi Jagadeesan
"""

from matriarch import *
import math


def collagen(seq1, seq2):
    a1 = chain(seq1)
    a2 = chain(seq2)
    hel1 = helix(a1, 1.5, 9.5238, "L")
    hel2 = helix(a2, 1.5, 9.5238, "L")
    helhel1 = helix(hel1, 4, 85.5, "L")
    helhel2 = helix(hel2, 4, 85.5, "L")
    helhel1rot = shiftOrbs(rotateOrbs(helhel1, 2 * math.pi / 3), 2.8)
    helhel2rot = shiftOrbs(rotateOrbs(helhel2, 4 * math.pi / 3), 5.6)
    homodimer = overlay(helhel1, helhel1rot)
    output = overlay(homodimer, helhel2rot)
    return output


seq1 = "GFZGPKGTAGEZGKAGERGVZGPZGAVGPAGKDGEAGAQGAZGPAGPAGERGEQGPA"
seq2 = "GFZGPKGPSGDZGKZGEKGHPGLAGARGAZGPDGNNGAQGPZGPQGVQGGKGEQGPA"
collgn = collagen(seq1, seq2)
fileOut(collgn, "10-collgn.pdb")
