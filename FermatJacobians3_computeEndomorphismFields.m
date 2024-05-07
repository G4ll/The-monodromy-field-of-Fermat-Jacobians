SetLogFile("EndomorphismsExceptionalCases.out");
prec := 100;

m := 20;
load "FermatJacobians1_computeMTandKconn.m";
F20 := compute_Kconn();
print F20;

R<t> := PolynomialRing(RationalsExtra(prec));
f := t^20 + 1;
X := HyperellipticCurve(f);
time lat := HeuristicEndomorphismLattice(X);
print lat;
L20 := NumberField(R!lat[3][1][2]);
IsIsomorphic(F20, L20);

m := 24;
load "FermatJacobians1_computeMTandKconn.m";
F24 := compute_Kconn();

R<t> := PolynomialRing(RationalsExtra(prec));
f := t^24 + 1;
X := HyperellipticCurve(f);
time lat := HeuristicEndomorphismLattice(X);
print lat;
L24 := NumberField(R!lat[3][1][2]);
IsIsomorphic(F24, L24);


m := 60;
load "FermatJacobians1_computeMTandKconn.m";
time gens := compute_Kconn_generators();

K60 := CyclotomicField(60);
polys60 := [];

for p in gens do
	fact := Factorisation(ChangeRing(p, K60));
	for f in fact do
		if Degree(f[1]) gt 1 then
			Append(~polys60, f[1]);
		end if;
	end for;
end for;

assert {Degree(f) : f in polys60} eq {2};
discriminants := {Discriminant(f) : f in polys60};
discriminants := [discr : discr in discriminants] ;
KummerGenerators := [];

for d in discriminants do
	ss := Subsets( {k : k in KummerGenerators} );
	test := true;
	for s in ss do
		if test then
		if #s gt 0 then
			if IsSquare(&*s * d) then
				test := false;
			end if;
		end if;
		end if;
	end for;
	if test then
		Append(~KummerGenerators, d);
	end if;
end for;

KummerGenerators;

/*
We describe simplified generators
*/
z := K60.1;

g1 := 3*z^12 + 2*z^8 - 2*z^6 - z^4 - 2*z^2 -1;
g2 := 4*z^13 - 2*z^3;
g3 := 2*z^14 - 2*z^9 -z^6 - 2*z^4 + 2*z^3 + 2;

IsSquare(g1*KummerGenerators[1]);
IsSquare(g2*KummerGenerators[3]);
IsSquare(g3*KummerGenerators[2]);

exit;
