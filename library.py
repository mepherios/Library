# Library . py 

# -----------------------------------------------------------------------------
# Lagre filer med pickle

# Prerequisites
import pickle
# Lagrer verdi til fil-navn, navn = 'navn'+'.pickle'
# Bruk -    Skriv navn som 'navn'
def pickle_store(name,value):
    namepickle=name+'.pickle'
    pickling_on = open(namepickle, 'wb')
    pickle.dump(value,pickling_on)
    pickling_on.close()

# Bruk - Skriv navn som 'navn', Lagre verdien som variabel med navn for å bruke verdien.
def pickle_open(name):
    pickle_name=name+'.pickle'
    pickling_on = open(pickle_name,'rb')
    pickle_value = pickle.load(pickling_on)
    return pickle_value
# ----------------------------------------------------------------------------------------------------
# Returnerer punktet for verdien -1 i listen.
import numpy as np

a=np.array([1,3,5,7])
b=np.argwhere(a==5)

# Gir første kryssning å teste for som krysser med y verdien i f.
import numpy as np
'''
a=np.array([1,3,5,4])
b=np.array([2,4,6,8])

c=(np.isin(a,b,assume_unique=True))
e=np.nonzero(c)
# Gir indeks nummeret som er på kurven.
f=e[0]

print(e)

x=[]
y=np.array(x)
'''
#
# -------------------------------------------------------------------------------------------------------
# Wallet import format osv. https://github.com/crcarlo/btcwif
'''
import btcwif
priv=hex(3) # Fjern x for å få svaret. Og 50/50 for 0 i starten.

priv='01'
a=3
print(a)
wif=btcwif.privToWif(priv)
print(wif)
'''

# ------------------------------------------------------------------------------------------------------
# Ecc Add given Point a, curve and S.

# Print function that gives values in a list for testing.
# !/usr/bin/env python3

import collections
import random

EllipticCurve = collections.namedtuple('EllipticCurve', 'name p a b g n h')

curve = EllipticCurve(
    'secp256k1',
    # Field characteristic.
    p=0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f,
    # Curve coefficients.
    a=0,
    b=7,
    # Base point.
    g=(0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798,
       0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8),
    # Subgroup order.
    n=0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141,
    # Subgroup cofactor.
    h=1,
)


'''
curve = EllipticCurve(
    'test31',
    # Field characteristic.
    p=31,
    # Curve coefficients.
    a=0,
    b=7,
    # Base point.
    g=(20,
       3),
    # Subgroup order.
    n=21,
    # Subgroup cofactor.
    h=1,
)
'''

# Modular arithmetic ##########################################################

def inverse_mod(k, p):
    """Returns the inverse of k modulo p.
    This function returns the only integer x such that (x * k) % p == 1.
    k must be non-zero and p must be a prime.
    """
    if k == 0:
        raise ZeroDivisionError('division by zero')

    if k < 0:
        # k ** -1 = p - (-k) ** -1  (mod p)
        return p - inverse_mod(-k, p)

    # Extended Euclidean algorithm.
    s, old_s = 0, 1
    t, old_t = 1, 0
    r, old_r = p, k

    while r != 0:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t

    gcd, x, y = old_r, old_s, old_t

    assert gcd == 1
    assert (k * x) % p == 1

    return x % p


# Functions that work on curve points #########################################

def is_on_curve(point):
    """Returns True if the given point lies on the elliptic curve."""
    if point is None:
        # None represents the point at infinity.
        return True

    x, y = point

    return (y * y - x * x * x - curve.a * x - curve.b) % curve.p == 0


def point_neg(point):
    """Returns -point."""
    assert is_on_curve(point)

    if point is None:
        # -0 = 0
        return None

    x, y = point
    result = (x, -y % curve.p)

    assert is_on_curve(result)

    return result


def point_add(point1, point2):
    """Returns the result of point1 + point2 according to the group law."""
    assert is_on_curve(point1)
    assert is_on_curve(point2)

    if point1 is None:
        # 0 + point2 = point2
        return point2
    if point2 is None:
        # point1 + 0 = point1
        return point1

    x1, y1 = point1
    x2, y2 = point2

    if x1 == x2 and y1 != y2:
        # point1 + (-point1) = 0
        return None

    if x1 == x2:
        # This is the case point1 == point2.
        m = (3 * x1 * x1 + curve.a) * inverse_mod(2 * y1, curve.p)
    else:
        # This is the case point1 != point2.
        m = (y1 - y2) * inverse_mod(x1 - x2, curve.p)

    x3 = m * m - x1 - x2
    y3 = y1 + m * (x3 - x1)
    result = (x3 % curve.p,
              -y3 % curve.p)

    assert is_on_curve(result)

    return result



def point_add2(point1, point2):
    """Returns the result of point1 + point2 according to the group law."""
    assert is_on_curve(point1)
    assert is_on_curve(point2)

    if point1 is None:
        # 0 + point2 = point2
        return point2
    if point2 is None:
        # point1 + 0 = point1
        return point1

    x1, y1 = point1
    x2, y2 = point2

    if x1 == x2 and y1 != y2:
        # point1 + (-point1) = 0
        return None

    if x1 == x2:
        # This is the case point1 == point2.
        m = (3 * x1 * x1 + curve.a) * inverse_mod(2 * y1, curve.p)
    else:
        # This is the case point1 != point2.
        m = pow(((y1 - y2) * inverse_mod(x1 - x2, curve.p)),2,curve.p)

    x3 = pow((m * m - x1 - x2),2,curve.p)
    y3 = pow((y1 + m * (x3 - x1)),2,curve.p)
    result = (x3 % curve.p,
              -y3 % curve.p)

    assert is_on_curve(result)

    return result




def scalar_mult(k, point):
    """Returns k * point computed using the double and point_add algorithm."""
    assert is_on_curve(point)

    if k % curve.n == 0 or point is None:
        return None

    if k < 0:
        # k * point = -k * (-point)
        return scalar_mult(-k, point_neg(point))

    result = None
    addend = point

    while k:
        if k & 1:
            # Add.
            result = point_add(result, addend)

        # Double.
        addend = point_add(addend, addend)

        k >>= 1

    assert is_on_curve(result)

    return result


# Keypair generation and ECDHE ################################################

# Bruk denne til å lage alle verdiene.
def make_keypair():
    """Generates a random private-public key pair."""
    # private_key = random.randrange(1, curve.n)
    private_key = 3
    public_key = scalar_mult(private_key, curve.g)

    return private_key, public_key



# Inverts value given mpdulus
def inv(x,modulus):
    #inverse=pow(x,modulus-2,modulus)
    inverse=p-x
    return inverse

# makes invertion value of list given modulus
def list_swap(list,modulus):
    modulatedlist=[]
    for value in list:
        modulated=inv(value,modulus)
        modulatedlist.append(modulated)
    return modulatedlist

# Checks for duplicates in list1 vs list 2.
def checkIfDuplicates(listOfElems,setOfElems):
    ''' Check if given list contains any duplicates '''    
    for elem in listOfElems:
        if elem in setOfElems:
            #print('True')
            return True
        else:
            pass         
    #return False



