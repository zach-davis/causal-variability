"""
Generating queries automatically

import os
from PIL import Image
import numpy as np

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
"""

def query_generator(vals):
	base = Image.open('base.png')

	# cause
	X = Image.open('X/X_' + vals[0] + '.png')
	X.thumbnail((275,275))

	# effects
	img_size = (345,345)
	Y = Image.open('Y/Y_' + vals[1] + '.png')
	Y.thumbnail(img_size)
	Z = Image.open('Z/Z_' + vals[2] + '.png')
	Z.thumbnail(img_size)

	base.paste(X, (0,30), X)
	base.paste(Y, (372,0), Y)
	base.paste(Z, (372,62), Z)
	base.show()

query_generator(['query', 'normal', 'unknown'])










