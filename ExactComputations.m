CheckRelations := false;
m := 15;

"m =", m;
/*
Let m, d be positive integers with d | m.
This badly-written function computes the matrix
of the norm map from Q(\zeta_m)* to Q(\zeta_d)*,
where both are considered as the Q-rational points
of the algebraic tori Res_{Q(\zeta_m)/Q}(G_m) and 
Res_{Q(\zeta_d)/Q}(G_m)
*/
function MatrixOfNorm(m, d)
	r := EulerPhi(d);
	c := EulerPhi(m);
	M := ZeroMatrix(Integers(), r, c);

	ieff := 0;

	for i in [ u : u in [1..(d-1)] | GCD(u, d) eq 1 ] do
		ieff := ieff + 1;
		jeff := 0;
		for j in [ v : v in [1..(m-1)] | GCD(v, m) eq 1 ] do
			jeff := jeff +1;
			if (j mod d) eq i then
				M[ieff][jeff] := 1;
			end if;
		end for;
	end for;
	return Transpose(M);
end function;

/*
(Also badly written) helper function: given
an integer i with 0 < i < d and (i,d)=1,
return the integer k such that i is the k-th
integer in the interval [1, d-1] that is 
coprime to d
*/
function IndexToIndex(d, i)
	jeff := 0;
	for j in [ u : u in [1..(d-1)] | GCD(u, d) eq 1 and 2*u ne d ] do
		jeff := jeff + 1;
		if j eq i then
			return jeff;
		end if;
	end for;
end function;

/*
Compute the matrix of the reflex norm
corresponding to the CM type
{ i : i<d/2, (i,d)=1 }
*/
function MatrixOfReflexNorm(d)
	r := EulerPhi(d);
	M := ZeroMatrix(Integers(), r, r);
	ieff := 0;

	CMType := [ i : i in [1..d] | i lt d/2 and GCD(i, d) eq 1 ];
	ReflexCMType := [ InverseMod(i, d) : i in [1..d] | i lt d/2 and GCD(i, d) eq 1 ];

	for i in [ u : u in [1..(d-1)] | GCD(u, d) eq 1 ] do
		ieff := IndexToIndex(d, i);
		for k in ReflexCMType do
			M[IndexToIndex(d, (i*k) mod d)][ieff] := 1;
		end for;
	end for;
	return M;
end function;

/*
A "partial" Mumford-Tate matrix, obtained as the product
of the matrices of the norm and of the reflex norm
*/
function PartialMTMatrix(m, d)
	return MatrixOfNorm(m, d) * MatrixOfReflexNorm(d);
end function;

/*
Deligne's formula to attach a complex number to a collection
of numbers in the interval [1, m], seen as a character of
a product of groups (Z/mZ).
*/
function ComputeGammaProduct(exponents, m)
	C<i> := ComplexField(600);
	return &*[Gamma(C!j/m) : j in exponents] / (2*Pi(C)*i)^(&+exponents/m);
end function;

/*
Shioda's bijection (which to i associates the triple i, i, (-2i) mod m). I think this only applies to odd values of m.
*/
function TriplicateExponents(exponents, m)
	result := [];
	for i in exponents do
		result := result cat [ i, i, Integers()!( (-2*i) mod m ) ];
	end for;
	return result;
end function;

/*

*/
function Reorder(list, order)
  result := [];
  for i in [1..#order] do
  	for j in [1..#order] do
  		if order[j] eq i then
  			Append(~result, list[j]);
  		end if;
  	end for;
  end for;
  return result;
  end function;

/*
The full linear-algebraic information about 
the Mumford-Tate group of Jac(y^2 = x^m-1).
The kernel of this matrix, suitably interpreted,
gives equations for the Mumford-Tate group.
*/
function BuildMTMatrix(m)
	Ds := Divisors(m);
	Ds := [i : i in Ds | i ne 1 and i ne 2];
	r := EulerPhi(m);
	ReturnMatrix := Matrix(Integers(), r, 0, []);
	for d in Ds do
		PartialMTd := PartialMTMatrix(m, d);
		ReturnMatrix := HorizontalJoin(ReturnMatrix, PartialMTd);
	end for;
	return ReturnMatrix;
end function;

/*
The order in which characters appear in the matrix
of the Mumford-Tate group.
*/
function ListOrderCharacters(m)
	Ds := Divisors(m);
	Ds := [i : i in Ds | i ne 1 and i ne 2];
	Ordering := [];
	for d in Ds do
		for i in [1..d] do
			if GCD(i, d) eq 1 then
				Append(~Ordering, Integers()!((m/d)*i));
			end if;
		end for;
	end for;

	/*
	When m is even we need to shift the indices over m/2 by 1
	*/
	if m mod 2 eq 0 then
		for j in [1..#Ordering] do
			if Ordering[j] gt m/2 then
				Ordering[j] := Ordering[j]-1;
			end if;
		end for;
	end if;
	return Ordering;
end function;


"** Computing equations for the Mumford-Tate group **";
// 

MTMatrix := BuildMTMatrix(m);
B := Basis(Kernel(Transpose(MTMatrix)));
OrderCharacters := ListOrderCharacters(m);


"Matrix describing the Mumford-Tate group", MTMatrix;
"Basis of B";
OrderCharacters;


/* Here's a bit of code to print the MT equations */

function print_MT_equations()
	PP := PolynomialRing(Rationals(), m-1);
	RR := RationalFunctionField(PP);
	AssignNames(~PP, ["x_{" cat Sprint(j) cat "}": j in [1..(m-1)]]);
	
	all_equations := [];

	for v in B do
		eqn := RR!1;
		for i in [1..(m-2)] do
		    j := OrderCharacters[i];
		    if j ge (m/2) then
		        j := j+1;
		    end if;
		    if not(v[i] eq 0) then
		        eqn := eqn * (RR!(PP.j))^v[i];
		    end if;
		end for;
		Append(~all_equations, eqn);
	end for;
	
	return all_equations;
end function;





if CheckRelations then
"** Testing relations **";
// We now check that the eigenvalues of Frobenius,
// as computed using Kedlaya's algorithm,
// satisfy the equations found above.
// More precisely, they should satisfy said equations
// up to roots of unity, so the output should only
// include roots of unity (in my tests I only ever
// saw \pm 1)

R<x> := PolynomialRing(Rationals());
f := x^m-1;
C := HyperellipticCurve(f);

p := m+1;

for i in [1..5] do
	p := NextPrime(p);
	while not (p mod m eq 1) do
		p := NextPrime(p);
	end while;
	"Prime p =", p;
	Qp := pAdicField(p,30);
	M := ChangeRing(FrobeniusMatrix(C, p),Qp);

	if m mod 2 eq 0 then
		M := Matrix(Qp, #Rows(M)-1, #Rows(M)-1, [ M[i,j] : i in [ u : u in [1..#Rows(M)] | u ne Integers()!( (#Rows(M)+1)/2 ) ], j in [ u : u in [1..#Rows(M)] | u ne Integers()!( (#Rows(M)+1)/2 ) ] ]);
	end if;

	eigenvalues := [ M[i][i] : i in [1..#Rows(M)] ];

	for b in B do
		r2 := Reorder(b, OrderCharacters);
		&*[ eigenvalues[i]^r2[i] : i in [1..#r2] ];
	end for;
end for;

end if;





/*
Functions to compute the minimal polynomial of a product
of Gamma-factors algebraically.
*/


K<z> := CyclotomicField(4*m);
PrimeFactorsm := PrimeDivisors(m);
R<sqrttwopi> := PolynomialRing(K, 1+#PrimeFactorsm);	// the variables from the second to the last correspond to the prime divisors of m. The variable corresponding to a prime p represents p^(1/2m). The first variable is \sqrt{2\pi}
S<sqrttwopi> := FieldOfFractions(R);
MRelations := ZeroMatrix(Rationals(), m, Floor(m/2));
C := ComplexField(600);

EvaluationVector := [ Sqrt(2*Pi(C)) ] cat [ (C!k)^(1/(2*m)) : k in PrimeFactorsm ] ;

CorrectionFactors := [];
CorrectionFactorsC := [];

/*
Returns n mod m, with the convention that if 
n is 0 mod m, then m is returned instead of 0
*/
function SpecialMod(n, m)
	r := n mod m;
	if r ne 0 then
		return r;
	else
		return m;
	end if;
end function;


function ShiftCorrectionFactor(d, z, m)
	u := d*z/m;
	SCF := 1;
	while u gt 1 do
		SCF := SCF * (u-1);
		u := u-1;
	end while;
	return SCF;
end function;

/*
The sine of \pi x / m
*/
function ExactSine(x, m)
	return 1/2 * 1/z^m * ( z^(2*x) - z^(-2*x) );
end function;

/*
Given an integer d = p_1^{e_1} ... p_k^{e_k} Q,
where p_1, ..., p_k are the prime factors of d
that also divide m, return the monomial
x_{p_1}^{e_1} ... x_{p_k}^{e_k}
*/
function MonomialFromNumber(d)
	exps := [Valuation(d,p) : p in PrimeFactorsm];
	return &*[ S.(i+1)^exps[i] : i in [1..#exps] ];
end function;


/*
Compute and store the relations between the symbols
Gamma(i/m) generated by the reflection formula
Gamma(i/m)Gamma( (m-i)/m ) = \pi / \sin( \pi x / m  ).
The relation is stored in the following format:
- a vector in MRelations with 1 in positions i/m and (m-i)/m (or possibly, a 2 in position m/2)
- the corresponding "correction factor" \pi / \sin(\pi x/m), both as a complex number (stored in CorrectionFactorsC) 
  and as an exact number (stored in CorrectionFactors). Note that, since \pi is transcendental, we should (and do)
  store it as an indeterminate
*/
for x in [1..Floor(m/2)] do
	MRelations[x][x] := 1;
	MRelations[Integers()!((-x) mod m)][x] := MRelations[Integers()!((-x) mod m)][x] + 1; // might want to replace 1 with 2, not 0 with 1
	Append(~CorrectionFactorsC, Pi(C) / Sin(Pi(C) * x/m) );
	Append(~CorrectionFactors, 1/2 * (sqrttwopi)^2 / ExactSine(x, m) );
end for;


/*
Compute and store the relations among the symbols Gamma(i/m)
that arise from the multiplication theorem. For a prime factor
p of m, we represent p^(1/2m) as a variable in the polynomial ring R.
This allows us to work in a reasonably small number field,
instead of extending scalars all the way to Q(\zeta_{4m}, p_1^{1/2m},
..., p_k^{1/2m}).
*/
for d in Divisors(m) do
	if d ne 1 then
		codiv := Integers()!(m/d);
		for x in [1..codiv] do
			Mx := ZeroMatrix(Rationals(), m, 1);
			for a in [1..d] do // a corresponds to a/d = a(N/d) / N
				Mx[ SpecialMod(codiv*a+x, m) ][1] := 1;
			end for;
			Mx[SpecialMod(d*x, m)][1] := Mx[SpecialMod(d*x, m)][1] - 1;
			MRelations := HorizontalJoin(MRelations, Mx);
			exponent := Integers()!( m*(1-2*d*x/m) );
			Append(~CorrectionFactorsC, (2*Pi(C))^((d-1)/2) * d^( (1-2*d*x/m)/2 ) * ShiftCorrectionFactor(d, x, m) );
			// Append(~CorrectionFactors, (2*pi)^((d-1)*m) * d^( (m-2*d*x) ) * ShiftCorrectionFactor(d, x, m)^(2*m) );
			Append(~CorrectionFactors, (sqrttwopi)^(d-1) * MonomialFromNumber(d)^( (m-2*d*x) ) * ShiftCorrectionFactor(d, x, m) );
		end for;
	end if;
end for;

/*
Given an integer vector v, representing a Z-divisor supported on
(Z/mZ) \ {0}, the next function tries to express v as a Q-linear
combination of the distribution relations previously computed,
and (if successful) returns:
- a pair (u, denominator), where u = Gamma(v)^denominator
  and u is an exact algebraic integer (with fractional powers
  of the prime divisors of m represented as variables in the
  ring R)
- a complex number approximating Gamma(v), computed via the
  distribution relations;
- a complex number approximating Gamma(v), computed numerically
  from the definition
*/
function ComputeExactGammaProduct(v)
	/*
	We express the Gamma-product we have to compute in terms 
	of (ratios of) Gamma-products that we know how to express
	exactly in terms of powers of \pi and algebraic numbers
	*/
	vM := Matrix(Rationals(), m, 1, v cat [0]);
	M2 := HorizontalJoin(vM, MRelations);
	w1 := Basis(Kernel(Transpose(M2)))[1];
	w1 := [ -w1[i] : i in [2..(#CorrectionFactors+1)] ];
	// "Relation:", w1;
	denominator := LCM([Denominator(j) : j in w1]);
	w1Rescaled := [Integers()!(denominator*w1[i]) : i in [1..#w1]] ;
	Result1 := &*[CorrectionFactorsC[i]^(w1[i]) : i in [1..#w1]];	// computation via reflection formula + multiplication theorem
	Result2 := &*[Gamma(C!i/m)^v[i] : i in [1..(m-1)] ];			// direct computation in MAGMA
	u := &*[CorrectionFactors[i]^(w1Rescaled[i]) : i in [1..#w1]];	// this should be (value we want to compute)^(denominator)

	return u, denominator, Result1, Result2;
end function;


/*
Wrapper for ComputeExactGammaProduct. The only difference is the
input format: here the divisor is specified as a list where the
number of occurrences of an integer i is the coefficient of [i/m]
in the divisor. In other words, for m=3, the divisor [1/3] + 5*[2/3]
would be represented as [1, 5] in the case of the previous function,
and as [1, 2, 2, 2, 2, 2] in the case of the present one.
*/
function ComputeExactGammaProductFromList(list, m)
	v := [ 0 : i in [1..(m-1)] ];
	for i in list do
		v[i] := v[i]+1;
	end for;
	return ComputeExactGammaProduct(v);
end function;

/*
Given two lists listNumerator and listDenominator, expressed in the
same format as the input of ComputeExactGammaProductFromLists, this
function returns
ComputeExactGammaProductFromList(listNumerator) / ComputeExactGammaProductFromList(listDenominator)
*/
function ComputeExactGammaRatioFromLists(listNumerator, listDenominator, m)

	v := [ 0 : i in [1..(m-1)] ];
	for i in listNumerator do
		v[i] := v[i]+1;
	end for;

	for i in listDenominator do
		v[i] := v[i]-1;
	end for;

	u, denominator, Result1, Result2 := ComputeExactGammaProduct(v);
	return u, denominator, Result1, Result2;
end function;


function InterpretMonomialAsAlgebraicNumber(f, m)
	mon := Monomials(f);
	assert #mon eq 1;
	coeff := Coefficients(f)[1];
	
	exponents := Exponents(mon[1]);
	CommonDenominator := LCM([Denominator(e/(2*m)) : e in exponents]);
	
	radical := &*[ PrimeFactorsm[i]^( Integers()!( CommonDenominator/(2*m) * Exponents(mon[1])[i+1] ) ) : i in [1..#PrimeFactorsm] ];
	
	U<t> := PolynomialRing(K);
	f := t^CommonDenominator - radical;
	L := SplittingField(f);
	r := Roots(ChangeRing(f, L))[1][1];

	return coeff * r;
end function;


function InterpretRationalFunctionAsAlgebraicNumber(f, m)
	if f eq 0 then
		return 0;
	end if;
	n := Numerator(f);
	d := Denominator(f);
	return InterpretMonomialAsAlgebraicNumber(n,m) / InterpretMonomialAsAlgebraicNumber(d,m);
end function;















"** Generate Gamma-products **";
// Each equation corresponds to a character,
// and every character corresponds to a
// product of Gamma-factors, as in Deligne.
// We explicitly write down these (ratios of)
// products of Gamma factors, compute a minimal
// polynomial for each of them, and take the
// compositum of all the splitting fields.
// Hopefully, this should yield K^{conn}

K := QNF();
for b in B do
	bReorder := Reorder(b, OrderCharacters);
	// "Computing Gamma-product for relation", b;
	// "Reordered relation", bReorder;
	"Computing Gamma-product for relation", bReorder;
	NumeratorFactors := [];
	DenominatorFactors := [];
	e := 0;
	for i in [1..(2*Floor(((m-1)/2)))] do

		if m mod 2 eq 0 and i ge m/2 then
			ieff := i+1;
		else
			ieff := i;
		end if;

		if bReorder[i] gt 0 then
			NumeratorFactors := NumeratorFactors cat [ ieff : j in [1..bReorder[i]] ];
		else
			DenominatorFactors := DenominatorFactors cat [ ieff : j in [1..Abs(bReorder[i])] ];
		end if;
		e := e + ieff * bReorder[i];
	end for;

	elements12 := TriplicateExponents(NumeratorFactors, m);
	elements22 := TriplicateExponents(DenominatorFactors, m);


	u := ComputeGammaProduct(elements12, m) / ComputeGammaProduct(elements22, m);


	uExact, denExact, Result1, Result2 := ComputeExactGammaRatioFromLists(elements12, elements22, m);

	if m mod 2 eq 0 then
		uExact := uExact * (z^2 * 1/S.2^4)^(-e * denExact);	// S.2 is 2^(1/(2m)), so S.2^(-4) is 2^(-2/m) = 4^(-1/m)
	end if;

	// "Exact representation", uExact;
	// "Exponent of representation", denExact;

	// "Exact representation as algebraic number", InterpretRationalFunctionAsAlgebraicNumber(uExact, m);


	fExact := MinimalPolynomial(InterpretRationalFunctionAsAlgebraicNumber(uExact, m));
	if denExact gt 1 then
		gExact := Evaluate(fExact, Parent(fExact).1^denExact);
		fExact := Factorisation(gExact)[1][1];
	end if;

	"Minimal polynomial from exact representation as algebraic number", fExact;


	if m mod 2 eq 0 then
		// "Exponent e =", e;
		CC := Parent(u);
		u := u * ( Exp( 2*Pi(CC)*CC.1/(2*m) ) * (CC!4)^(-1/m) )^(-e);
	end if;
	
	f := MinimalPolynomial(u, 2*EulerPhi(m));
	"Minimal polynomial obtained numerically", f;

	F := OptimisedRepresentation(SplittingField(fExact));
	K := Compositum(K, F);
end for;

K := Compositum(K, CyclotomicField(m));

"K^{conn}:", K;


