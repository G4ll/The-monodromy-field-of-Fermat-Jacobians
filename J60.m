SetLogFile("J60.out");

// Base field
K<z> := CyclotomicField(60);

// The three quadratic generators
delta1 := 3*z^12 + 2*z^8 - 2*z^6 - z^4 - 2*z^2 - 1;
delta2 := 4*z^13 - 2*z^3;
delta3 := 2*z^14 - 2*z^9 - z^6 - 2*z^4 + 2*z^3 + 2;

deltas := [delta1, delta2, delta3];

// Adjoin one square root to a number field L.
function AddSqrt(L, a)
    P<t> := PolynomialRing(L);
    return NumberField(t^2 - L!a);
end function;

// Trivial extension
K0 := K;

// Degree-2 extensions
K1<a1> := AddSqrt(K, delta1);
K2<a2> := AddSqrt(K, delta2);
K3<a3> := AddSqrt(K, delta3);

// Degree-4 extensions
K12<b12> := AddSqrt(K1, delta2);
K13<b13> := AddSqrt(K1, delta3);
K23<b23> := AddSqrt(K2, delta3);

// Degree-8 extension
K123<c123> := AddSqrt(K12, delta3);

split_primes := [ 1021, 1741, 1861 ];
non_split_primes := [ 61, 181, 241, 421, 541, 601, 661, 1201, 1321, 1381, 1621, 1801 ];

for L in [K1, K2, K3] do
	Labs := AbsoluteField(L);
	Labs := OptimisedRepresentation(Labs);
	O := MaximalOrder(Labs);
	"Number of primes in L above the split primes";
	for p in split_primes do
		I := ideal<O | p>;
		"The split prime", p, "factors into", #Factorisation(I), "primes in L";
	end for;
	"Number of primes in L above the non-split primes";
	for p in non_split_primes do
		I := ideal<O | p>;
		"The non-split prime", p, "factors into", #Factorisation(I), "primes in L";
	end for;
end for;