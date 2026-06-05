// Characteristic polynomial det(Frob_p - T) for
// C : y^2 = x^60 + 1,
// for primes p congruent to 1 modulo 60.
//
// The roots are
//
//     -chi^a(-1) * J(chi^a, phi)
//
// for 1 <= a <= 59, a != 30, where phi = chi^30.

SetLogFile("FrobCharPolysJ60.out");

function FrobCharpolyC60Jacobi(p : ReturnJacobiSums := false)

    K<z> := CyclotomicField(60);
    P<T> := PolynomialRing(K);

    // Let g be a primitive root mod p, and define chi(g)=z.
    g := Integers()!PrimitiveElement(GF(p));

    // Discrete logarithm table modulo p with respect to g.
    logs := [ -1 : i in [1..p] ];
    u := 1;
    for e in [0..p-2] do
        logs[u+1] := e;
        u := (u*g) mod p;
    end for;

    // counts[r+1][s+1] counts x in F_p with
    // log(x) = r mod 60 and log(1-x) = s mod 2.
    counts := [ [0, 0] : r in [0..59] ];

    for x in [1..p-1] do
        if x ne 1 then
            y := (1 - x) mod p;

            r := logs[x+1] mod 60;
            s := logs[y+1] mod 2;

            row := counts[r+1];
            row[s+1] +:= 1;
            counts[r+1] := row;
        end if;
    end for;

    zpow := [ z^i : i in [0..59] ];

    h := ((p-1) div 2) mod 60; // chi(-1)=z^h

    A := [ a : a in [1..59] | a ne 30 ];

    J := AssociativeArray(Integers());
    eigenvalues := AssociativeArray(Integers());

    fK := P!1;

    for a in A do
        Ja := K!0;

        for r in [0..59] do
            for s in [0..1] do
                c := counts[r+1][s+1];
                if c ne 0 then
                    e := (a*r + 30*s) mod 60;
                    Ja +:= c*zpow[e+1];
                end if;
            end for;
        end for;

        // Frobenius eigenvalue:
        // alpha_a = -chi^a(-1) J(chi^a, phi).
        chi_minus_one_a := zpow[((a*h) mod 60)+1];
        alpha := -chi_minus_one_a * Ja;

        J[a] := Ja;
        eigenvalues[a] := alpha;

        fK *:= T - alpha;
    end for;

    // Coerce the polynomial back to Z[T].
    Z := Integers();
    PZ<TZ> := PolynomialRing(Z);

    coeffs := [];
    for i in [0..Degree(fK)] do
        c := Coefficient(fK, i);
        ok, ci := IsCoercible(Z, c);
        Append(~coeffs, ci);
    end for;

    fZ := PZ!coeffs;

    if ReturnJacobiSums then
        return fZ, J, eigenvalues;
    else
        return fZ;
    end if;

end function;


// Given a Weil polynomial f(T) = det(T - Frob_q), compute:
//   p        : characteristic,
//   q        : size of the base field,
//   n        : minimal extension degree over F_q,
//   N        : minimal extension degree over F_p,
//   orders   : root-of-unity orders occurring among alpha_i/alpha_j.
//
// Thus all geometric endomorphisms are defined over F_{q^n}=F_{p^N}.
//
// Optional parameter:
//   SplitField := L
// If supplied, L should be a number field over which f splits.

function BaseFieldDataFromFrobCharpoly(f)
    P := Parent(f);

    d := Degree(f);
    g := d div 2;

    c := Coefficient(f, 0);
    c := Abs(Integers()!c);

    ok, q := IsPower(c, g);

    q := Integers()!q;

    fac := Factorization(q);

    p := fac[1][1];
    a := fac[1][2];

    return p, q, a;
end function;


function RootsWithMultiplicitiesOver(f, L)
    Qx<x> := PolynomialRing(Rationals());
    fQ := Qx!f;

    PL<T> := PolynomialRing(L);
    roots := Roots(PL!fQ);

    return roots;
end function;


function TateDimensionFromRoots(roots, n)
    vals := [];
    mults := [];

    for r in roots do
        beta := r[1]^n;
        m := r[2];

        found := false;
        for i in [1..#vals] do
            if beta eq vals[i] then
                mults[i] +:= m;
                found := true;
                break;
            end if;
        end for;

        if not found then
            Append(~vals, beta);
            Append(~mults, m);
        end if;
    end for;

    return &+[ m^2 : m in mults ];
end function;


function EndomorphismFieldFromFrobCharpoly(f : SplitField := false,
                                              ReturnField := false,
                                              Verbose := false)

    p, q, a := BaseFieldDataFromFrobCharpoly(f);

    Qx<x> := PolynomialRing(Rationals());
    fQ := Qx!f;

    if Type(SplitField) eq BoolElt then
        if Verbose then
            "Computing a splitting field...";
        end if;
        L := SplittingField(fQ);
    else
        L := SplitField;
    end if;

    roots := RootsWithMultiplicitiesOver(fQ, L);

    orders := { Integers() | 1 };

    for i in [1..#roots] do
        alpha := roots[i][1];

        for j in [1..#roots] do
            beta := roots[j][1];

            rho := alpha / beta;

            isrou, ord := IsRootOfUnity(rho);
            if isrou then
                Include(~orders, Integers()!ord);
            end if;
        end for;
    end for;

    n := LCM([ m : m in orders ]);

    // Base field is F_q = F_{p^a}; full endomorphism field is F_{q^n}.
    N := a*n;

    dim_base := TateDimensionFromRoots(roots, 1);
    dim_geom := TateDimensionFromRoots(roots, n);

    if Verbose then
        printf "Base field: F_%o = F_%o^%o\n", q, p, a;
        printf "Orders of root-of-unity ratios: %o\n", Sort(Setseq(orders));
        printf "Minimal degree over F_p: %o\n", N;
        printf "dim End^0 over base field: %o\n", dim_base;
        printf "dim geometric End^0: %o\n", dim_geom;
    end if;

    if ReturnField then
        kEnd := GF(p^N);
        return p, q, N, kEnd, orders, dim_base, dim_geom;
    else
        return p, q, N, orders, dim_base, dim_geom;
    end if;
end function;






/*
Expected dimension endomorphism algebra: 2*2^2 + 2*3^2 + 4*2^2 + 8*2^2 + 2*4^2 + 4*4^2;
170
*/







Ps := PrimesUpTo(2000);
Ps := [p : p in Ps | p mod 60 eq 1];
time chpolys := [FrobCharpolyC60Jacobi(p) : p in Ps];

split_primes := [];
non_split_primes := [];

for cp in chpolys do
	p, q, N, orders, dim_base, dim_geom := EndomorphismFieldFromFrobCharpoly(cp : Verbose := true);
	if dim_geom eq 170 then
		if dim_geom eq dim_base then
			Append(~split_primes, p);
		else
			Append(~non_split_primes, p);
		end if;
	end if;
end for;