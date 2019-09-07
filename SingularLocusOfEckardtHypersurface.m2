"SingularLocusOfEckardtHypersurface",
    	Version => "1.10", 
    	Date => "September, 2019",
    	Authors => {{Name => "Hanieh Keneshlou", 
		  Email => "hanieh.keneshlou@mis.mpg.de", 
		  HomePage => "https://keneshlou.wordpress.com"
		  Headline => "Cubic surfaces on singular locus of the Eckardt hypersurface",
                "This Macaulay2 file contains supporting code for the computational proofs 
		 in the paper 'Cubic surfaces on singular locus of the Eckardt hypersurface' "
		  }}

VerifyAssertion1:

kk=QQ--field of definition
R=kk[a_0..a_4];--coordinate ring of P4
--The elementary symmetric polynomials and invariants
s1=a_0+a_1+a_2+a_3+a_4;
s2=sum flatten apply(5,i->apply(5,j->(if i<j then a_i*a_j else 0)));
s3=sum flatten apply(5,i->flatten(apply(5,j->apply(5,k->(if i<j and j<k then a_i*a_j*a_k else 0)))));
s4=a_0*a_1*a_2*a_3+a_0*a_1*a_2*a_4+a_4*a_1*a_2*a_3+a_0*a_1*a_4*a_3+a_0*a_4*a_2*a_3;
s5=a_0*a_1*a_2*a_3*a_4;---
i8=s4^2-4*s3*s5;
i16=s5^3*s1;
i24=s5^4*s4;
i32=s5^6*s2;
i40=s5^8;----
L=apply(flatten entries vars R,i->{1,i,i^2,i^3,i^4})
M=matrix{L_0,L_1,L_2,L_3,L_4};
v=det M;
I=ideal(s5^18*v);--Eckardt hypersurface of degree 100
Delta=ideal jacobian gens I+I;--singular locus
dim Delta, degree Delta
dc=decompose Delta--decomposition to irreducible components
--The singular locus is not equidimensional.
--It has 5 component of dimension 3 and 25 components of dimension 2.


--We compute the map phi and the image of each component under phi
--Note:without lose of generality, as components differ by coordinate permutation
--in P4 we only compute the image of two random representative components of each group
T=kk[A,B,C,D,E,Degrees=>{1,2,3,4,5}]--Coordinate ring of P(1,2,3,4,5)
phi=map(R,T,{i8,i16,i24,i32,i40});
--The two surfaces as components of the image of singular locus
time Ic1=preimage_phi dc_10;--S1
dim Ic1,degree Ic1
time Ic2=preimage_phi dc_15;--S_2
dim Ic2,degree Ic2
J=Ic1+Ic2;--The intersection curve C1\cup C1
dim J, degree J


h=dc_6+dc_21;
mingens h--the line mapping to the rational curve C1
dim h, degree h, genus h
H=preimage_phi h;--the curve C1
dim H, degree H


p=dc_14+dc_15;
mingens p--the second line mapping to the rational curve C2
P=preimage_phi p;--the curve C2
dim P, degree P

Pt=P+H;--intersection of two curves C1 and C2
dim Pt, degree radical Pt
G=decompose Pt--the two points corresponding to the Clebsch surface and Fermat cubic surface



VerifyAssertion2:

D2=ideal jacobian gens Delta+Delta;
isTriple=apply(dc,l->isSubset(D2,l))--true means the component is a triple singularity
--we check the type of singularities determined by generic element of the two rational curves
isSubset(D2,h)
isSubset(D2,p)--true means the lines are triple singularity
