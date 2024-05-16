/* ------------------------------------------------------------ */
/*
We compute the endomorphism field of the excpetional factors
X20, X24, X60, of Corollary 3.5.5
*/
/* ------------------------------------------------------------ */

SetLogFile("Corollary355.out");
AttachSpec("./CHIMP/CHIMP.spec");
prec := 100;

m := 20;
load "FermatJacobians1_computeMTandKconn.m";
"** Computing the connected monodromy field of J20**";
time F20 := compute_Kconn();
"The connected monodromy field of J20 is the", F20;

R<t> := PolynomialRing(RationalsExtra(prec));
f := t^20 + 1;
X := HyperellipticCurve(f);
"** Computing the Endomorphism lattice of J20 via the methods of [CMSV19] **";
time lat := HeuristicEndomorphismLattice(X);
print lat;
L20 := NumberField(R!lat[3][1][2]);
"Do the two computation give the same Q(EndJm)?";
IsIsomorphic(F20, L20);

/* ------------------------------------------------------------ */

m := 24;
load "FermatJacobians1_computeMTandKconn.m";
"** Computing the connected monodromy field of J24**";
F24 := compute_Kconn();
"The connected monodromy field of J24 is the", F24;

R<t> := PolynomialRing(RationalsExtra(prec));
f := t^24 + 1;
X := HyperellipticCurve(f);
"** Computing the Endomorphism lattice of J24 via the methods of [CMSV19] **";
time lat := HeuristicEndomorphismLattice(X);
print lat;
L24 := NumberField(R!lat[3][1][2]);
"Do the two computation give the same Q(EndJm)?";
IsIsomorphic(F24, L24);

/* ------------------------------------------------------------ */

m := 60;
load "FermatJacobians1_computeMTandKconn.m";
"** Computing generators for the connected monodromy field of J60**";
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

"Are the numbers in the articles generators for the connected monodromy field of J60?";
IsSquare(g1*KummerGenerators[1]);
IsSquare(g2*KummerGenerators[3]);
IsSquare(g3*KummerGenerators[2]);

exit;
