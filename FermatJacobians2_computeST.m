/* ------------------------------------------------------------ */
/* ------------------- 0. Global set up ----------------------- */
/* ------------------------------------------------------------ */


load "FermatJacobians1_computeMTandKconn.m";

K := SplittingField(CyclotomicPolynomial(m));

/* ------------------------------------------------------------ */

"** Converting equations for the Mumford-Tate group **";
all_f_list := usable_char_MT_equations();

/* ------------------------------------------------------------ */

"** Computing the connected monodromy field **";
Kconn := compute_Kconn();

_, KtoKconn := IsSubfield(K,Kconn);
Kconn_polynomial := DefiningPolynomial(Kconn);									
iota := hom< Kconn -> C | Roots(ChangeRing(Kconn_polynomial, C))[11][1] >;	
GaloisKconn, AutKconn, GaloisMap := AutomorphismGroup(Kconn);	

/* ------------------------------------------------------------ */

/* We define a big helper matrix whose entries are
algebraically indipendent variables */

ring_matrix_coeff := PolynomialRing(Kconn,(2*g)^2);
AssignNames(~ring_matrix_coeff, ["t" cat Sprint(i) cat "/" cat Sprint(j) : i in [1..(2*g)], j in [1..(2*g)]]);

matrix_wannabe := [[ring_matrix_coeff.(i+j*2*g): j in [0..(2*g-1)]]: i in [1..(2*g)]];
T := Matrix(ring_matrix_coeff, 2*g, 2*g, matrix_wannabe);



/* ------------------------------------------------------------ */
/* ---------- 1. Computing tensor representation -------------- */
/* ------------------------------------------------------------ */


/* Helper function:
if tau zeta^i = zeta^j,
tells me j. */
function ActionOnZeta(tau, i)
    for j in [0..(m-1)] do
        y := (GaloisMap(tau)(K.1^i)) / K.1^j;
        if y eq 1 then
            return j;
        end if;
    end for;
end function;


/* This function computes the algebraic
number P(gamma, omega_alpha) and stores it
as an element of the connected monodromy field Kconn.
The input is a *non-tripled* character (i_1, i_2, ... i_q),
as it takes the output of usable_char_MT_equations() */
function P_gamma(alpha)
	gamma_factor 	:= &*[ Gamma(C!(i/m))^2 / Gamma(C!(2*i/m)) : i in alpha];
    pi_factor 		:= ( 2*C.1*Pi(C) )^( Integers()!(#alpha/2) );
    algebraic_number := gamma_factor / pi_factor;
    
    min_poly, _ 	:= MinimalPolynomial(algebraic_number, 2*m);
    rr 				:= Roots(min_poly, Kconn);
    residues 		:= [Norm(iota(x[1]) - algebraic_number): x in rr];
    _, minimal_index := Min( residues );
    true_gamma 		:= rr[minimal_index][1];
    
    return true_gamma;
end function;


/* The rational coefficients coming
from Poincaré duality */
function mu(alpha)
	return &*[ (m-2*i)/m : i in alpha ];
end function;


/* The entry in position [alpha, u^(-1)alpha]
in the tensor representation matrix rho(tau) */
function STcoeff(tau, alpha)
	u := ActionOnZeta(tau, 1);
	uu := Modinv(u, m);
	
	a_start 	:= [ -i mod m : i in alpha ];
	a_target 	:= [ -uu*i mod m : i in alpha ];
	
	mu_start 	:= mu( a_start );
	mu_target 	:= mu( a_target );
	
	P_start 	:= P_gamma( a_start );
	P_target	:= P_gamma( a_target );
	
	return (( mu_target * GaloisMap(tau)(P_target) ) / ( mu_start * P_start ));
end function;


/* Given a character, build the tensor-representation
matrix rho(tau) on the [alpha]-generalized eigenspace */
function small_a_matrix(tau, alpha)
	M := ZeroMatrix(Kconn, m-1);
	u := ActionOnZeta(tau, 1);
	uu := Modinv(u, m);
		
	for b in coprimes_to_m do
		M[b, uu*b mod m] := STcoeff(tau, [b*i mod m : i in alpha]);
	end for;
	
	return M;
end function;



/* ------------------------------------------------------------ */
/* ---------------- 2. Computing relations -------------------- */
/* ------------------------------------------------------------ */


/* Given a character, compute the relation it imposes
on the matrix rho(tau) of the actual
(non-tensor) Galois representation */
function relation_from_character(tau, alpha)
	u := ActionOnZeta(tau, 1);
	uu := Modinv(u, m);
	
	monomial := &*[ T[i, uu*i mod m] : i in alpha ];
	constant := STcoeff(tau, alpha);
	
	return (monomial - constant);
end function;


/* Given a character, compute the relations that
ALL characters in the alpha-orbit impose
on the matrix rho(tau) of the actual
(non-tensor) Galois representation */
function relations_from_character_eigenspace(tau, alpha)
	relations := [];
	for x in coprimes_to_m do
		xa := [x*i mod m : i in alpha];
		x_rel := relation_from_character(tau, xa);
		Append(~relations, x_rel);
	end for;
	
	return relations;
end function;



/* ------------------------------------------------------------ */
/* ----------- 3. Computing connected components -------------- */
/* ------------------------------------------------------------ */


/* This function eats an element tau of GaloisKconn
and computes the ideal corresponding to tau-component in the
Sato-Tate group of Jm
in the polynomial ring where the variables are the entries
of a general 14x14 matrix */
function compute_component(tau)
	u := ActionOnZeta(tau, 1);
	uu := Modinv(u, m);
	TT := 0 * T;
	
	for i in [1..(2*g)] do
		j := i*uu mod m;
		TT[i, j] := T[i, j];
	end for;
	
	I := ideal< ring_matrix_coeff | 0 >;
	PermutationStructure := [ [i,j] : i in [1..(2*g)], j in [1..(2*g)] | j mod m ne (i*uu) mod m ];

  	Iperm := ideal< ring_matrix_coeff | [ T[sigma[1], sigma[2]] : sigma in PermutationStructure ] >;
  
  	I := I + Iperm;
	// I := ideal< ring_matrix_coeff | [&*[TT[i*uu^j mod m,i*uu^(j+1) mod m] : j in [1..8]]-1 : i in [1..m-1]] >;
	// I := ideal< ring_matrix_coeff | Determinant(TT)-1 >;
	// I := ideal< ring_matrix_coeff | Determinat(TT)-1, TT^8-1 >;
	Power8 := {(TT^8)[i][i]-1 : i in [1..(2*g)]};
	I := I + ideal< ring_matrix_coeff | [f : f in Power8] >;
	
	all_a_list := [ f[1] cat [ -x mod m : x in f[2]] : f in all_f_list ];
	for a in all_a_list do
		a_generators := relations_from_character_eigenspace(tau, a);
		J := Ideal(a_generators);
		I := I + J;
	end for;
	
	Groebner(I);
	return I;
end function;



/* ------------------------------------------------------------ */
/* ------------------- 4. The case m = 15 --------------------- */
/* ------------------------------------------------------------ */

function computeST15()
	if m ne 15 then
		return "m is not 15";
	end if;
	
	
	/* Check if the identity component
	is in the right isomorphism class for m=15 */
	function check_identity()
		tau := Identity(GaloisKconn);
		I := compute_component(tau);

		I_guess_gen := [T[5,5]-T[13,13]*T[3,3]*T[4,4], T[6,6]-T[14,14]*T[3,3]*T[4,4], T[7,7]-T[14,14]*T[13,13]*T[3,3]^2*T[4,4] ] cat [T[i,i]*T[m-i, m-i]-1 : i in [1..7]];
		I_guess := Ideal(I_guess_gen);

		return I_guess eq I;
	end function;


	"** Computing the identity component of ST **";
	check_identity();


	"** Computing the connected component of ST for the first generator **";
	tau := (GaloisKconn.1);
	I1 := compute_component(tau);
	J1 := I1;
	_, indipendent_variables := Dimension(J1);
	J1 := J1 + ideal< ring_matrix_coeff | [ring_matrix_coeff.t-1 : t in indipendent_variables] >;
	Groebner(J1);
	"Matrix in the component of GaloisKconn.1:", J1;
	non_trivial_relations := [];
	for b in Basis(J1) do
		if #Terms(b) ge 2 then
			Append(~non_trivial_relations, b);
		end if;
	end for;
	non_trivial_relations1 := non_trivial_relations;

	"** Computing the connected component of ST for the second generator **";
	tau := (GaloisKconn.2);
	I2 := compute_component(tau);
	J2 := I2;
	_, indipendent_variables := Dimension(J2);
	J2 := J2 + ideal< ring_matrix_coeff | [ring_matrix_coeff.t-1 : t in indipendent_variables] >;
	Groebner(J2);
	"Matrix in the component of GaloisKconn.2:", J2;
	non_trivial_relations := [];
	for b in Basis(J2) do
		if #Terms(b) ge 2 then
			Append(~non_trivial_relations, b);
		end if;
	end for;
	non_trivial_relations2 := non_trivial_relations;	
	
	return I1, J1, non_trivial_relations1, I2, J2, non_trivial_relations2;
end function;

/* ------------------------------------------------------------ */
/* --------------- 5. Useful out-put functions ---------------- */
/* ------------------------------------------------------------ */


"This file of code has 3 useful output functions:";
"	1. compute_component(tau);";
"	2. small_a_matrix(tau, alpha);";
"	3. computeST15();";



