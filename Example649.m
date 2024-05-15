load "FermatJacobians2_computeST.m";


pi := Pi(C);
x1 := P_gamma([9,12,7,2]);



function algebrify(algebraic_number)
	min_poly, _ := MinimalPolynomial(algebraic_number, 10, 1000);
    rr := Roots(min_poly, Kconn);
    residues := [Norm(iota(x[1]) - algebraic_number): x in rr];
    _, minimal_index := Min( residues );
    true_gamma := rr[minimal_index][1];
    
    return true_gamma;
end function;

function GammaNew(alpha1, alpha2)
	num := &*[ Gamma(C!(i/m))^2 * Gamma(C!((-2*i mod m)/m)) : i in alpha1];
    den := &*[ Gamma(C!(i/m))^2 * Gamma(C!((-2*i mod m)/m)) : i in alpha2];
    complex_gamma_value := num / den;
    algebraic_gamma_value := algebrify(complex_gamma_value);
    
    return algebraic_gamma_value;
end function;


/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */

x1 := GammaNew([9,12],[8,13]);
IsSubfield(SplittingField(MinimalPolynomial(x1)), K);

y1 := (1/2) * C!Sqrt(3/5) * C!(1/Sin(3*pi/15));
y1 := algebrify(y1); 
x1 eq y1;

/* ------------------------------------------------------------ */

x2 := GammaNew([11, 12], [9, 14]);
IsSubfield(SplittingField(MinimalPolynomial(x2)), K);

y2 := (1/8) * C!Sqrt(3) * C!(1/Sin(4*pi/15))^2 * C!(1/Sin(7*pi/15));
y2 := algebrify(y2);
x2 eq y2;

/* ------------------------------------------------------------ */

x3 := GammaNew([10, 12], [8, 14]);
IsSubfield(SplittingField(MinimalPolynomial(x3)), K);
IsSubfield(SplittingField(MinimalPolynomial(x3)), Kconn);

y3sq := (3/16) * Sqrt(5) * Sin(pi/15)^3 * Sin(6*pi/15) / ( Sin(2*pi/15) * Sin(3*pi/15)^4 * Sin(5*pi/15)^3 );
y3sq := algebrify(y3sq);
x3^2 eq y3sq;

