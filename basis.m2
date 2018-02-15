loadPackage "Polyhedra"

expList = (I) -> flatten exponents first flatten entries gens I

algGens = (I,J)->(
  B:=(expList(J))_(positions(expList(J),i->i!=0));
  A:=(expList(I))_(positions(expList(J),i->i!=0));
  L:=sort apply(A,B,(i,j)->i/j);
  C:=flatten {0,apply(L,i->numerator i),1};
  D:=flatten {1, apply(L,i->denominator i),0};
  M:=matrix{C,D};
  apply(unique flatten apply(#C-1,
                             i->hilbertBasis(posHull submatrix(M,{i,i+1}))),
        i->reverse flatten entries i)
)

S = QQ[x]
I = ideal(x)
algGens(I^3, I^2)
