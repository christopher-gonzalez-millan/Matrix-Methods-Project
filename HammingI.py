
from scipy.linalg import null_space
import numpy as np

def hamming(a):
    b = str(a)
    cw=list(b)
    k,m=len(cw),0
    while k+m+1>2**m:
        m+=1
    if k+m+1!=2**m: print("invalid number of bits to construct hamming code")
    """IDK why but this if command doesnt seem to be working properly for (11 bit)"""


    Ik= np.identity(k)
    PP=np.ones((k,m))
    li = False
    while (li == False):
        PP = np.ones((k, m))
        for x in range(0, m):
            z = 0
            for z in range(0, k - m):
                PP[np.random.randint(0, k - 1), x] = 0
                z += 1
        b = null_space(PP)
        li= np.array_equal(b,np.empty((m,0)))
    G=np.hstack((Ik,PP))
    """Here we create the Generator Matrix that will out put the original k message bits with additional m parity bits
       appended. We construct the functions for each parity bit by generating m linearly independent vectors that will contain
       all possible combinations of bit pairs as its rows.The first four columns will return the original 4 bits and thus
       are the I(k) matrix, the resulting columns will define the bits and will be the aforementioned linearly independent
       basis for all of the parity check bits"""
    x=np.asarray(cw,dtype=int) ##was attempting to make dtype boolean here, but the calculations wherent working properly
    ##G.dtype=bool## STRUGGLING WITH IMPLEMENTING BIT ALGEBRA by making all the array's dtype=bool
    bitsfortransmission = np.matmul(x,G)
    print(bitsfortransmission)

def decode(d):
    r = str(d)
    e = list(r)

    received= np.asarray(cw, dtype=int)
    Pcheck=np.hstack((np.matrix.transpose(PP),np.identity(m)))
    """We create the Parity Check Matrix, the orthogonal complement to the space defined by G. This is expedited using 
    the definition provided by _________"""
    parity=np.matmul(Pcheck,np.transpose(received))
    if parity == np.empty((m,0))
        print(received[0:3]) ####print first four elements of the input hamming code word
    else:
        w=Pcheck.index(np.transpose(parity))
        if received[w]==1: received[w]=0:
        print(received[0:3])
        elif received[w]==0: \received[w]=1
        print(received[0:3])











hamming(1101)
"""Idk how to fix leading decimal issue either for integers when we are putting them into the function,
tried to put other stuff into it but it hasn't been working"""

