"""Below we will construct a program to correct error in a [k,n]-Hamming code"""

from scipy.linalg import null_space
import numpy as np

a= input('Bit Message to be encoded: ')
b = str(a)
cw=list(b)
k,m=len(cw),0
while k+m+1>2**m:
    m+=1
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
x=np.asarray(cw,dtype=int)
"""Here we create the Generator Matrix that will out put the original k message bits with additional m parity bits
appended. We construct the functions for each parity bit by generating m linearly independent vectors that will contain
all possible combinations of bit pairs as its rows.The first four columns will return the original 4 bits and thus
are the I(k) matrix, the resulting columns will define the bits and will be the aforementioned linearly independent
basis for all of the parity check bits"""
print( 'Generator Matrix: \n ', G)
bitsfortransmission = np.matmul(x,G)
"""We now multiply a (1xk) matrix containing the code word by the (kxm) Generator Matrix to return the encoded message"""
for i in range(0, k + m):
    if bitsfortransmission[i] % 2 == 0:
        bitsfortransmission[i] = 0
    elif bitsfortransmission[i]%2 == 1:
        bitsfortransmission[i] = 1
print(bitsfortransmission)
"""This is a very round about way of implementing and/or bit algebra."""
print('Encoded message: ', bitsfortransmission)

###Introducing error to the code ###
amountoferrors= 1
for x in range(0,amountoferrors):
    bitsfortransmission[np.random.randint(0,k+m)]=np.random.randint(0,2)
"""Random error is introduced in the specified amount. Note: Hamming Code's can only correctly 
identify and correct up to 1 error"""

d = bitsfortransmission
e = d.tolist()
if len(e)!= (k+m):
    print('in/del present, cannot be corrected')
else:
    received=np.asarray(e,dtype=int)
    Pcheck=np.hstack((np.matrix.transpose(PP),np.identity(m)))
    """We create the Parity Check Matrix, the orthogonal complement to the space defined by G. This is expedited using 
    the definition provided by _________"""
    parity=np.matmul(Pcheck,np.transpose(received))
    """Multiplying the parity check matrix by the transposed of the """
    p= np.zeros((m,))
    for i in range(0,m):
        if parity[i]%2==0:
            p[i]=0
        elif parity[i]%2==1:
            p[i]=1
    """Once again, another round about way to implemnent and/or bit algebra"""


    if np.array_equal(p, np.zeros((m,))) == True :
        print('no error detected')
        print('message sent: ', received[0:k])
        """If the syndrome returns a zero vector, we know that it exists in the kernal of the matrix which by the 
        Fredholm Alternitive we know is the image of its orthogonal compliment, the Generator matrix and therefore is a
        valid codeword"""
    else:
        f = np.transpose(Pcheck)
        error = False
        errorlocation = 0
        for x in range(0, k + m):
            error = np.array_equal(p, f[x])
            if error == True : errorlocation = x
        """Here we located the position which the error occurred in transmission by comparing each column of the 
        Parity check matrix to the syndrome received"""
        if received[errorlocation]==0:
            received[errorlocation]=1
        elif received[errorlocation]==1:
            received[errorlocation]=0
        """Now that the location of the error has been identified, the bit at its location can be corrected"""
        print('error occured in ' , errorlocation + 1, 'position of received code')
        print('message sent: ', received[0:k])


