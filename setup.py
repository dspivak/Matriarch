import parseAminoAcids

from distutils.core import setup
setup(name="Matriarch",
	version = '1.0',
	description = 'A python library for materials archtecture',
	author = 'Ravi Jagadeesan, Tristan Giesa, David Spivak, Markus Buehler',
        url='http://www.web.mit.edu/matriarch',
        packages=["matriarch"]
	)

import os
try:
	os.remove("matriarch/aminoAcidData.py")
except:
	print("Failed to remove aminoAcidData.py. Continuing...")
try:
	os.remove("parseAminoAcids.pyc")
except:
	print("Failed to remove parseAminoAcids.pyc. Continuing...")