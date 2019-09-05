   kk=QQ;--Field of definition
  E=kk[e_0..e_3,SkewCommutative=>true];--The exterior algebra
  m=matrix{{e_0*e_1+e_2*e_3}};--the 1*1n map in Tate resolution of G
  betti(fm=res coker m)--the betti tabel of the resolution
          
 R=kk[x_0..x_3]--coordinate ring of P^3
 A=sum(4,i->x_i*sub(contract(E_i,fm.dd_4),R));-- The 16*35 matrix
B=sum(4,i->x_i*sub(contract(E_i,fm.dd_5),R)); -- The 35*64 matrix
C=A*B--the 16*64 matrix
E9=kk[e_0..e_9,SkewCommutative=>true]
x2=flatten entries basis(2,R)
M=sum(10,i->e_i*sub(contract(x2_i,C),E9))
      --the 16*64 matrix whose transpose is phi_0
M'=(syz transpose M)**E9^{-3};--the 16*16 syzygy matrix
M'+transpose M--not 0, M' is not skewsymmetric
 ---We reform M' to get skew-symmetric kernel
S=kk[a_0..a_119]
A=genericSkewMatrix(S,a_0,16)
SE=S**E9
AM=sub(A,SE)*sub(M,SE);
N=map(E9^1024,,sub(transpose(flatten AM//sub(vars S,SE)),E9));
sAM'=syz(N,DegreeLimit=>5)
Psi=map(E9^16,,sub(A, flatten transpose sAM'))
Psi+transpose Psi--is 0, is skewsymmetric



kk=QQ;
AR=kk[a_0..a_19][x_0..x_3];
F=sub(vars(kk[a_0..a_19]), AR)*transpose basis(3,AR);
Q=jacobian F;
J=jacobian transpose Q;--jacobian matrix
Q16=flatten Q**matrix{{x_0,x_1,x_2,x_3}};--16 cubic equations
Dj=det J;
f4=transpose jacobian matrix{{Dj}};--4 differentials of the Jacobian
al20=Q16|f4;--the 20 cubic forms
(C,m) = coefficients al20;
A=kk[a_0..a_19];
M=sub(m,A);--The 20*0 matrix
L=ideal random(A^1,A^{17:-1});--choose 17 random linear form to cut 
Mq=sub(M,A/L);
B=kk[a_17..a_19]--new ring with remaining variables
MB=sub(Mq,B)--the 20*20 matrix in remaining variables
time Dis=ideal det MB;--the discrimninant polynomial of degree 32
dim Dis, degree Dis--(2,32)
Delta=ideal jacobian gens Dis;--singular locus of the discriminant
dim Delta, degree Delta
     --the singular locus of degree 520 is 280 nodes and 120 cusps
time V2=minors(19,MB);-- used 3299.31 seconds
dim V2, degree V2--zero-dim'l consists of 400 points
isSubset(Delta,V2)
 --true means V2 is the singular locus of discriminant curve
