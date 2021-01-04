"""
Generating queries automatically
"""

import os
from PIL import Image
import numpy as np

# move file into images/domain folder for this
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)


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

 #args are X=Y, Y=X1, Z=X2 in notation used for setup experiment slides
 #econ values: Low, High, Small
query_generator(['query', 'normal', 'low'])


# to save the image instead of opening them.
# make sure to set the file folder correctly (now inferences1)
def query_generator2(vals, infnr, domain):
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
    imname = domain+"inf"+str(infnr+1)+vals[0]+vals[1]+vals[2]+'.png'
    base.save("inferences1\\"+imname)
    


allinfs = np.array([['nonnormal', 'query', 'nonnormal'],
           ['nonnormal', 'query', 'unknown'],
           ['nonnormal', 'query', 'normal'],
           ['nonnormal', 'query', 'nonnormal'],
           ['nonnormal', 'query', 'unknown'],
           ['nonnormal', 'query', 'normal'],
           ['nonnormal', 'nonnormal', 'query'],
           ['nonnormal', 'unknown', 'query'],
           ['nonnormal', 'normal', 'query'],
           ['nonnormal', 'nonnormal', 'query'],
           ['nonnormal', 'unknown', 'query'],
           ['nonnormal', 'normal', 'query'],
           ['query', 'nonnormal', 'nonnormal'],
           ['query', 'nonnormal', 'unknown'],
           ['query', 'nonnormal', 'normal'],
           ['query', 'nonnormal', 'nonnormal'],
           ['query', 'nonnormal', 'unknown'],
           ['query', 'nonnormal', 'normal'],
           ['query', 'nonnormal', 'nonnormal'],
           ['query', 'unknown', 'nonnormal'],
           ['query', 'normal', 'nonnormal'],
           ['query', 'nonnormal', 'nonnormal'],
           ['query', 'unknown', 'nonnormal'],
           ['query', 'normal', 'nonnormal']
          ])


def allinfsvals(allinfs, domainvals):
    #changes nonnormal vals in allinfs to specified vals
    #allinfs list of lists defined above
    #domainvals the nonnormal values as [Y, X1, X2]
    A = allinfs
    for i in range(24):
        for j in range(3):
            if A[i,j]=='nonnormal':
                A[i,j] = domainvals[j]
    return A        
    
def all_query_generator(allinfs, domainvals, domain):
    #make all images based on matrix of inf vals
    A = allinfsvals(allinfs, domainvals)
    for i in range(24):
        query_generator2(A[i], i, domain)

domain = 'econ1' #domain and inference set nr (for if we do multiple sets with differing valences)
econvals1 = ['low', 'high', 'low'] 

all_query_generator(allinfs, econvals1, domain)
