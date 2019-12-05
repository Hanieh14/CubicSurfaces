"DetRepOfCubicDiscriminant",
    	Version => "1.10", 
    	Date => "September, 2019",
    	Authors => {{Name => "Hanieh Keneshlou", 
		  Email => "hanieh.keneshlou@mis.mpg.de", 
		  HomePage => "https://keneshlou.wordpress.com/"},
		  {Name => "Dominic Bunnett", 
		  Email => "bunnett@math.tu-berlin.de", 
		  HomePage => "http://page.math.tu-berlin.de/~bunnett/"
		 "This Macaulay2 file contains supporting codes for the computational proofs 
                 in the paper 'Determinantal representations of the cubic discriminant' " }}


-------------------------------------------------------------------------
--Construction of the Tate resolution of the Null-Correspondence bundle--
-------------------------------------------------------------------------


kk=QQ;--Field of definition
E=kk[e_0..e_3,SkewCommutative=>true];
   --The exterior algebra
m=matrix{{e_0*e_1+e_2*e_3}};
   --the 1*1n map in Tate resolution of G
betti(fm=res coker m)
   --the betti table of the resolution

R=kk[x_0..x_3]--coordinate ring of P^3
A=sum(4,i->x_i*sub(contract(E_i,fm.dd_4),R));-- The 16*35 matrix
B=sum(4,i->x_i*sub(contract(E_i,fm.dd_5),R)); -- The 35*64 matrix
C=A*B--the 16*64 matrix

E9=kk[e_0..e_9,SkewCommutative=>true]
x2=flatten entries basis(2,R)
M=sum(10,i->e_i*sub(contract(x2_i,C),E9))
--the 16*64 matrix whose transpose is phi_0
M'=(syz transpose M)**E9^{-3};--the 16*16 syzygy matrix
M'+transpose M'--not 0, M' is not skewsymmetric

---We reform M' to get skew-symmetric kernel
S=kk[a_0..a_119]
A=genericSkewMatrix(S,a_0,16)
SE=S**E9
AM=sub(A,SE)*sub(M,SE);
N=map(E9^1024,,sub(transpose(flatten AM//sub(vars S,SE)),E9));
sAM'=syz(N,DegreeLimit=>5)
Psi=map(E9^16,,sub(A, flatten transpose sAM'))--the resultant matrix in Plucker coordinates
Psi+transpose Psi--is 0, is skewsymmetric



--Now here we apply the computation to that of an actual cubic
AR=kk[a_0..a_19][x_0..x_3];
F=sub(vars(kk[a_0..a_19]), AR)*transpose basis(3,AR);
--F=random(kk^1,kk^20)*transpose basis(3,AR)--test example numerically
Q=jacobian F;

cfs=(coefficients(transpose Q,Monomials=>flatten entries basis(2,AR)))_1
sb=subsets({0,1,2,3,4,5,6,7,8,9},4)

AE=kk[a_0..a_19,e_0..e_9]
Ie=apply(sb,l->sub(e_(l_0)*e_(l_1)*e_(l_2)*e_(l_3),AE));
de=matrix{apply(sb,l->sub(det cfs^l,AE))}
c=apply(16,k->de*(coefficients(sub(Psi^{k},AE),Monomials=>Ie))_1);
PsA=c_0||c_1||c_2||c_3||c_4||c_5||c_6||c_7||c_8||c_9||c_10||c_11||c_12||c_13||c_14||c_15;
PsA+transpose PsA
Ar=kk[a_0..a_19]
Pi=sub(PsA,Ar)--The 16*16 matrix in 20 coefficients of general cubic form

-----------------------------------------Example
--The Chow form of P^1 in P^3

 kk=QQ
 R=kk[x_0,x_1]--coordinate ring of P1
 B=basis(2,R)--H0(O(2))
 B3=basis(3,R)--W
 B5=basis(5,R)--H0(O(5))
 
 ER=kk[e_0..e_3, SkewCommutative=>true][x_0,x_1]

 E=kk[e_0..e_3, SkewCommutative=>true]--exterior algebra
 basis(2,E)--6 plucker coordinates
 
 T=sub(B,ER)
 P=sub(B5,ER)
 c=sub(B3,ER)
 use ER
 F=(sub(matrix{{e_0..e_3}},ER)*transpose c)_(0,0)
 y=sub(F*T//P,E)--phi0 in the Tate resolution
 phi0=mingens ker y--3*3 matrix defining the non-zero matrix in the chow complex 
 phi1=mingens ker phi0--phi_{-1}
 mingens ker phi1--phi_{-2} 

Rab=kk[a_0..a_3,b_0..b_3]
sylvester = matrix{{a_0,a_1,a_2,a_3,0,0},{0,a_0,a_1,a_2,a_3,0},{0,0,a_0,a_1,a_2,a_3},{b_0,b_1,b_2,b_3,0,0},{0,b_0,b_1,b_2,b_3,0},{0,0,b_0,b_1,b_2,b_3}}
   --this is the sylvester matrix for two binary cubics.

B = E**Rab
P = subsets({0,1,2,3},2)
m = matrix{{a_0..a_3},{b_0..b_3}}
I = ideal apply(P,i->e_(i_0)*e_(i_1)-det m_i) -- lists the minors
Q=B/I
tate = sub(sub(phi0,Q),Rab)
det tate == det sylvester


----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------

Nanson2.

--We build the 20*20 matrix described by Nanson, whose determinant gives
--the resultant of four quaternary quardic forms. We then carry out some basic checks
--on this form.
restart
kk=ZZ/32749;
A=kk[a_0..a_19]
AR=kk[a_0..a_19][x_0..x_3]

F=sub(vars(kk[a_0..a_19]), AR)*transpose basis(3,AR) 
    --a generic cubic equation

Q=jacobian F
Q16=flatten Q**matrix{{x_0,x_1,x_2,x_3}};
    --These are 16 cubic equations in the x_i

(C1,m1) = coefficients Q16
degree sub(m1_1_1,A)
    --one can check that all the coeffiecients are degree one (or zero) in the a_i

H=jacobian transpose Q
    --take the second derivatives to give te Hessian Matrix
time Dh=det H;
f4=transpose jacobian matrix{{Dh}};
    --4 differentials of the determinant of the Hessian matrix

(C2,m2) = coefficients f4;
degree sub(m2_0,A)
    --one can check  that the coefficients are of degree 4 in the a_i   
     
al20=Q16|f4;
    --We collect the 16 cubic equations and the differentials 

(C,m) = coefficients al20;
M=sub(m,A);
    --Nanson's 20*20 matrix.



------------------------------------------------------------------
--------Check the normal form locus of binodal cubics P^11 in P^19

L1 =sub(ideal(a_0, a_1, a_2, a_3, a_9, a_15, a_18, a_19),A)
    --normal form of binodal cubics with nodes at (1,0,0,0) and (0,0,0,1)
    
L=ideal random(A^1,A^{9:-1}) + L1
    --gives a random plane in P^11

Mq=sub(M,A/L1);
B=kk[a_4..a_8,a_10..a_14,a_16,a_17]
    --coordinate ring of P^2 with 3 remaining variables
    
MB=sub(Mq,B);
    --the 20*20 matrix in 3 remaining variables
    
time V2=minors(19,MB)
    --uses 124984 seconds.
    --we compute the V2 locus in these normal forms
    
V2==0
    --true: thus V2 contains these binodal cubics in normal form.
    
    
    
----------Now lets compare intersections of V2 with random planes
----------with no extra conditions:


L=ideal random(A^1,A^{17:-1});
    --choose 17 random linear form to cut out V_i's
    
Mq=sub(M,A/L);
B=kk[a_17..a_19]
    --coordinate ring of P^2 with 3 remaining variables
    
MB=sub(Mq,B)
    --the 20*20 matrix in 3 remaining variables

discr =  matrix{{det MB}};
    --discriminant curve on the plane (degree 32)
    
Delta = ideal jacobian  discr;
dim Delta, degree Delta
    --the singular locus is 0 dimensional and of degree 520

time V2random=minors(19,MB);
    -- used 3299.31 seconds

dim V2random, degree V2random
    --zero-dim'l consists of 400 points
    
isSubset(Delta,V2random)
    --true means V2 is the singular locus of discriminant curve
    
    
    
--------Now we say a little about the locus V3,
--------we can conclude that its codimesion is at least 3.


V3=ideal(det(submatrix(MB,{0..17},{0..17})));
dim(V3+V2)
    --(V3\cap P2) is empty
   
