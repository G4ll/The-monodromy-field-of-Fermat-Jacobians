SetLogFile("X20.out");
/*
We know that X_{20} is 4-dimensional and geometrically isogenous to E^4 for some elliptic
curve E. Two objectives: we want to compute the j-invariant of E (hence its minimal field
of definition) and the minimal field of definition of an isogeny X_{20} --> E^4
*/


K := Rationals();
R<x> := PolynomialRing(K);

/*
Step 1: we check that among the quotients of C there are elliptic curves
whose j-invariant generates Q(\sqrt{5})
*/

C1 := HyperellipticCurve(x^11+x);
G, map, map2 := AutomorphismGroup(C1);

for g in G do
	gp := AutomorphismGroup(C1, [map(g)]);
	X := CurveQuotient(gp);
	X;
	Genus(X);
end for;

C2 := HyperellipticCurve(x^5 - 5*x^3 + 5*x);

Qbar := AlgebraicClosure(Rationals());
C2bar := ChangeRing(C2, Qbar);

G, map, map2 := AutomorphismGroup(C2bar);
Absolutize(Qbar);
L := SplittingField(AbsolutePolynomial(Qbar));
L := OptimisedRepresentation(L);
L;

C2L := ChangeRing(C2, L);
G, map, map2 := AutomorphismGroup(C2L);

for g in G do
	gp := AutomorphismGroup(C2L, [map(g)]);
	X := CurveQuotient(gp);
	if Genus(X) eq 1 then	
		E := EllipticCurve(X);
		j := jInvariant(E);
		"Minimal polynomial of j-invariant", MinimalPolynomial(j);
		"Discriminant of Q(j):", Discriminant(MaximalOrder(NumberField(MinimalPolynomial(j))));
		"What's the CM order?", HasComplexMultiplication(E);
	end if;
end for;



/*
Step 2: we check the lattice of endomorphism rings of J_{20} over various
subfields of the field L over which all endomorphisms are defined. Only 
over L do we find a factor of multiplicity 4.
We note that, since we already know the endomorphism algebra of J_{20} and
the method correctly detects all the endomorphisms of J_{20}, the result
of this calculation is unconditional.
*/
prec := 100;
R<t> := PolynomialRing(RationalsExtra(prec));
f := t^20 + 1;
X := HyperellipticCurve(f);
HeuristicEndomorphismLatticeDescription(X);

exit;
