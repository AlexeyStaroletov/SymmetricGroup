##### This file contains code for computing minimal polynomials of elements in irreducible representations of a finite group

n:=11;  ### degree of Sn or An

G := AlternatingGroup(n);     ### group
tbl := Irr(G);;              ### character table
cc := ConjugacyClasses(G);;  ### conjugacy classes
num := Size(cc);


############ this function calculates the polynomial by permutation whose roots are the products of t roots of cycle lengths

PolyByCycles := function(perm, k, n)
  local y, F, f, ord, i, j, eta, res, cycles, tmp, ord2, zet;
  ord := Order(perm);
  F := Field(E(ord));               ### use field Q[eta], where eta^ord=1
  y := Indeterminate(F, "y");
  eta := E(ord);                    ### primitive root
  res := Set([]);
  cycles := CycleLengths(perm, [1..n]);
  tmp := Combinations(cycles, k);
  for i in [1..Size(tmp)] do
    ord2 := Lcm(tmp[i]);
    zet := eta^(ord/ord2);
    for j in [1..ord2] do
      AddSet(res, zet^j);            ### add to set all the products of t roots of cycle lengths
    od;
  od;
  f:= y^0;
  for i in res do
    f := f*(y-i);           ### find the polynomial
  od;
  Print(f, "\n");           ### print the polynomial
end;

############ this function calculates the minimal polynomial of a permutation in i-th representation of tbl

PolyOfPerm := function(a, i)
  local root, y, F, f, ord, j, dim, k, xi;
  ord := Order(a);
  F := Field(E(ord));
  y := Indeterminate(F, "y");
  root := E(ord);
  f := y^0;
  for k in [0..ord-1] do
    xi := root^k;
    xi := xi^-1;
    dim := 0;
    for j in [1..ord] do
      dim := dim + xi^j*((a^j)^tbl[i]);              #### find the dimeinsion of eigenspace associated to xi
    od;
    if dim <> 0 then
      f := f * (y - xi^-1);
    fi;
  od;  
  Print("Poly = ", f);
end;

############ this function calculates the degree of the minimal polynomial of permutation in i-th representation of tbl
###########  we use the code from PolyOfPerm, but we only calculate the degree instead of the polynomial.

DegOfPoly := function(a, i)
  local root, y, F, f, ord, j, dim, deg, k, xi;
  ord := Order(a);
  F := Field(E(ord));
  y := Indeterminate(F, "y");
  root := E(ord);
  f := y^0;
  deg := 0;
  for k in [0..ord-1] do
    xi := root^k;
    xi := xi^-1;
    dim := 0;
    for j in [0..ord-1] do
      dim := dim + xi^j*((a^j)^tbl[i]);
    od;
    if dim <> 0 then
      deg := deg + 1;
    fi;
  od;
  return deg;
end;

############ this function calculates the polynomial by permutation whose roots are the products of t roots of cycle lengths multiplied by the sign
############ based on PolyByCycles

PolyByCyclesConj := function(perm, k, n)
  local  y, F, f, ord, i, j, eta, res, cycles, tmp, ord2, zet;
  ord := Order(perm);
  F := Field(E(ord));
  y := Indeterminate(F, "y");
  eta := E(ord);
  res := Set([]);
  cycles := CycleLengths(perm, [1..n]);
  tmp := Combinations(cycles, k);
  for i in [1..Size(tmp)] do
    ord2 := Lcm(tmp[i]);
    zet := eta^(ord/ord2);
    for j in [1..ord2] do
      AddSet(res, SignPerm(perm)*zet^j);
    od;
  od;
  f := y^0;
  for i in res do
    f := f*(y-i);
  od;
  Print(f, "\n");
end;

##### now try all premutations from G and all representations of tbl

for i in [1..num] do
  x := Representative(cc[i]);
  ord := Order(x);
  for k in [1..num] do
    if ()^tbl[k] < n then  #### do not check trivial, sign, standard and associated with standard representations
        continue;
    fi;
    deg := DegOfPoly(x, k);
    if deg < ord then
      Print("Permutation: ", x, ", degree of min poly = ", deg, " in representaion ", k, " with dimension ", ()^tbl[k],  "\n");
    fi;
  od;
od;
