/* This script computes the connected monodromy field
for m odd and composite in in the interval 3..111 */

SetLogFile("Example6410.out");

function ExtractKummerGenerators(m, gens)
	Km := CyclotomicField(m);
	polys := [];

	for p in gens do
		fact := Factorisation(ChangeRing(p, Km));
		for f in fact do
			if Degree(f[1]) gt 1 then
				Append(~polys, f[1]);
			end if;
		end for;
	end for;
	
	if #polys gt 0 then
		assert {Degree(f) : f in polys} eq {2};
	else
		return [];
	end if;
	discriminants := {Discriminant(f) : f in polys};
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

	return KummerGenerators;
end function;



m := 15;
"----------------------------";
load "FermatJacobians1_computeMTandKconn.m";
time gens := compute_Kconn_generators();
KG := ExtractKummerGenerators(m, gens);
"Number of generators:", #KG;
KG;

m := 21;
"----------------------------";
load "FermatJacobians1_computeMTandKconn.m";
time gens := compute_Kconn_generators();
KG := ExtractKummerGenerators(m, gens);
"Number of generators:", #KG;
KG;

m := 33;
"----------------------------";
load "FermatJacobians1_computeMTandKconn.m";
time gens := compute_Kconn_generators();
KG := ExtractKummerGenerators(m, gens);
"Number of generators:", #KG;
KG;

m := 35;
"----------------------------";
load "FermatJacobians1_computeMTandKconn.m";
time gens := compute_Kconn_generators();
KG := ExtractKummerGenerators(m, gens);
"Number of generators:", #KG;
KG;

m := 39;
"----------------------------";
load "FermatJacobians1_computeMTandKconn.m";
time gens := compute_Kconn_generators();
KG := ExtractKummerGenerators(m, gens);
"Number of generators:", #KG;
KG;

m := 45;
"----------------------------";
load "FermatJacobians1_computeMTandKconn.m";
time gens := compute_Kconn_generators();
KG := ExtractKummerGenerators(m, gens);
"Number of generators:", #KG;
KG;

m := 51;
"----------------------------";
load "FermatJacobians1_computeMTandKconn.m";
time gens := compute_Kconn_generators();
KG := ExtractKummerGenerators(m, gens);
"Number of generators:", #KG;
KG;


m := 55;
"----------------------------";
load "FermatJacobians1_computeMTandKconn.m";
time gens := compute_Kconn_generators();
KG := ExtractKummerGenerators(m, gens);
"Number of generators:", #KG;
KG;

m := 57;
"----------------------------";
load "FermatJacobians1_computeMTandKconn.m";
time gens := compute_Kconn_generators();
KG := ExtractKummerGenerators(m, gens);
"Number of generators:", #KG;
KG;

m := 63;
"----------------------------";
load "FermatJacobians1_computeMTandKconn.m";
time gens := compute_Kconn_generators();
KG := ExtractKummerGenerators(m, gens);
"Number of generators:", #KG;
KG;

m := 65;
"----------------------------";
load "FermatJacobians1_computeMTandKconn.m";
time gens := compute_Kconn_generators();
KG := ExtractKummerGenerators(m, gens);
"Number of generators:", #KG;
KG;

m := 69;
"----------------------------";
load "FermatJacobians1_computeMTandKconn.m";
time gens := compute_Kconn_generators();
KG := ExtractKummerGenerators(m, gens);
"Number of generators:", #KG;
KG;

m := 75;
"----------------------------";
load "FermatJacobians1_computeMTandKconn.m";
time gens := compute_Kconn_generators();
KG := ExtractKummerGenerators(m, gens);
"Number of generators:", #KG;
KG;

m := 77;
"----------------------------";
load "FermatJacobians1_computeMTandKconn.m";
time gens := compute_Kconn_generators();
KG := ExtractKummerGenerators(m, gens);
"Number of generators:", #KG;
KG;

m := 85;
"----------------------------";
load "FermatJacobians1_computeMTandKconn.m";
time gens := compute_Kconn_generators();
KG := ExtractKummerGenerators(m, gens);
"Number of generators:", #KG;
KG;

m := 87;
"----------------------------";
load "FermatJacobians1_computeMTandKconn.m";
time gens := compute_Kconn_generators();
KG := ExtractKummerGenerators(m, gens);
"Number of generators:", #KG;
KG;

m := 91;
"----------------------------";
load "FermatJacobians1_computeMTandKconn.m";
time gens := compute_Kconn_generators();
KG := ExtractKummerGenerators(m, gens);
"Number of generators:", #KG;
KG;

m := 93;
"----------------------------";
load "FermatJacobians1_computeMTandKconn.m";
time gens := compute_Kconn_generators();
KG := ExtractKummerGenerators(m, gens);
"Number of generators:", #KG;
KG;

m := 95;
"----------------------------";
load "FermatJacobians1_computeMTandKconn.m";
time gens := compute_Kconn_generators();
KG := ExtractKummerGenerators(m, gens);
"Number of generators:", #KG;
KG;

m := 99;
"----------------------------";
load "FermatJacobians1_computeMTandKconn.m";
time gens := compute_Kconn_generators();
KG := ExtractKummerGenerators(m, gens);
"Number of generators:", #KG;
KG;

m := 105;
"----------------------------";
load "FermatJacobians1_computeMTandKconn.m";
time gens := compute_Kconn_generators();
KG := ExtractKummerGenerators(m, gens);
"Number of generators:", #KG;
KG;

m := 111;
"----------------------------";
load "FermatJacobians1_computeMTandKconn.m";
time gens := compute_Kconn_generators();
KG := ExtractKummerGenerators(m, gens);
"Number of generators:", #KG;
KG;


exit;
