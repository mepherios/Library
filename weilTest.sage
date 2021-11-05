
#
# Finite field over extension field
F1.<a>=GF(2^7)
# Elliptic curve extension field, [algebraic constants] ax3,bx2,c,dx2,e
E1=EllipticCurve(F1,[0,0,1,1,1])
# Point one
P1=E1.random_point()

# 
F2.<b>=GF(2^28)
E2=EllipticCurve(F2,[0,0,1,1,1])
# Roots to compute homomorphism of F1 to F2
aa=F1.modulus().roots(F2)[0][0]
# Algebraic homomorphism for extending point one into complete extension field to compute the ate/weil pairing
phi=Hom(F1,F2)(aa)

# Point generation for point on curve 1
P1=E1.random_point()
# Point 2 to check if ate/weil pairing works.
R1=18*P1
# Extended P1, R1 to F2 and named R2,P2.
P2=E2(phi(P1.xy()[0]),phi(P1.xy()[1]))
R2=E2(phi(R1.xy()[0]),phi(R1.xy()[1]))

# E2 point of order of E1. ## ker(pi-q)
Q=145*E2.random_point()

# Point 1, point 2 of weil pairing. P2,R2 = points of E1, Q is pairing point.
alpha=P2.weil_pairing(Q,113)
beta=R2.weil_pairing(Q,113)
solution=beta.log(alpha)
if solution==18:
    print('Weil pairing is functional')
else:
    print("Weil pairing is non-functional")


