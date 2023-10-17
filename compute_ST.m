/* -- Title -- */



/* ------------- set-up -------------- */


load "action_on_Halpha.m"; 									// include the previous code, computing the tau-action matrix
// SetLogFile("ST15-ideal.out");

m := 15;													// redundant - better safe then sorry
g := 7;								// genus of the hyperelliptic curve y^2=x^m+1	

equations2 := [ [[i], [i]] : i in coprimes_to_m ];			// the degree 2 equations defining MT
eq4_easy := [ [[i,15-i],[7,8]] : i in [1..6] ];				// the easy degree 4 equations
eq4_hard := [ [[9, 12],[8, 13]], [[11, 12],[9,14]], [[10, 12],[8, 14]] ]; // the hard ones
equations4 := eq4_easy cat eq4_hard;						// all degree 4 equations for MT	

all_f_list := equations2 cat equations4;					// all equations for MT (notice this works for m = 15 only!)

ring_matrix_coeff := PolynomialRing(Kconn,(2*g)^2);			// ring with matrix-entries as indeterminates
AssignNames(~ring_matrix_coeff, ["t" cat Sprint(i) cat "/" cat Sprint(j) : i in [1..(2*g)], j in [1..(2*g)]]);

matrix_wannabe := [[ring_matrix_coeff.(i+j*2*g): j in [0..(2*g-1)]]: i in [1..(2*g)]];
T := Matrix(ring_matrix_coeff, 2*g, 2*g, matrix_wannabe);	// now T[i,j] is the variable t_i,j

elements_of_GalKconn := [ x : x in Set(GaloisKconn)];		// list of elements in the Galois group of Kconn/Q


/* ------------------ functions ---------------------- */
/* The following function eats an element tau of the Galois group of Kconn
and an equation for MT, and gives back the relations that a general
14x14 matrix has to respect given the input.
The one below, computes the data for all equations f of MT. */

function relations_from_equation(f, tau)
	MM := Transpose(compute_matrix(f, tau));				// compute the tau action on H[alpha]
	non_tripled_alpha := f[1] cat [-x mod m : x in f[2]];	// the original alpha on Jm (not on Xm)
	
	molt := ActionOnZeta(tau, 1);							// set to zero the entries we know from theory shoul be zero, to ease computations
	TT := 0 * T;											// on the generic matrix T
	for i in [1..14] do
		j := i*molt mod m;
		TT[i, j] := T[i, j];
	end for;
	
	generators := [];										// here i compute how T^tensor4 acts on H[alpha]
	for i in [1..#coprimes_to_m] do							// and impose it equal to MM
		for j in [1..#coprimes_to_m] do
			ui_alpha := [a*coprimes_to_m[i] mod m : a in non_tripled_alpha];
			uj_alpha := [a*coprimes_to_m[j] mod m : a in non_tripled_alpha];
			monomial := &*[TT[ui_alpha[k], uj_alpha[k]]: k in [1..#non_tripled_alpha]];
			generators := Append(generators, monomial - MM[i, j]);
		end for;
	end for;
	
	return generators;										// return a list of equations, a subset of the generators of the ideal of the tau component
end function;

function all_relations(tau)									// just loop the previous function over all f to retrieve a complete list of generators
	all_generators := [];
	for f in all_f_list do
		f_generators := relations_from_equation(f, tau);
		all_generators := all_generators cat f_generators;
	end for;
	return all_generators;
end function;


/* --------------------------------------- */ 
/* The following function eats an element tau of GaloisKconn
and computes the ideal corresponding to tau-component in the Sato-Tate group of Jm
in the polynomial ring where the variables are the entries
of a general 14x14 matrix */

function compute_component(tau)		
	molt := ActionOnZeta(tau, 1);							// set to zero the entries we know from theory shoul be zero, to ease computations
	TT := 0 * T;
	for i in [1..14] do
		j := i*molt mod m;
		TT[i, j] := T[i, j];
	end for;
	
	I := ideal< ring_matrix_coeff | Determinant(TT)-1 >;

	for f in equations4 do
		f_generators := relations_from_equation(f, tau);
		J := Ideal(f_generators);
		I := I + J;
	end for;
	
	for f in equations2 do
		f_generators := relations_from_equation(f, tau);
		J := Ideal(f_generators);
		I := I + J;
	end for;
	
	Groebner(I);
	return I;
end function;

/* -------------- check id component is the right one --------------------- */

function check_identity()
	I := compute_component(Identity(GaloisKconn));

	I_guess_gen := [T[5,5]-T[13,13]*T[3,3]*T[4,4], T[6,6]-T[14,14]*T[3,3]*T[4,4], T[7,7]-T[14,14]*T[13,13]*T[3,3]^2*T[4,4] ] cat [T[i,i]*T[m-i, m-i]-1 : i in [1..7]];
	I_guess := Ideal(I_guess_gen);

	return I_guess eq I;
end function;

/* --------------------------------------------------------------------- */
/* the following function build a matrix in the tau component
of the Sato-Tate group */

function build_matrix_in_component(tau)
	I := compute_component(tau);
	J:= I + ideal< ring_matrix_coeff | [T[i, ActionOnZeta(tau, i)] -1 : i in [1,2,3,5] ] >;
	Groebner(J);
	pols := Basis(J);

	molt := ActionOnZeta(tau, 1);
	T_tau := ZeroMatrix(Kconn, 14, 14);
	for i in [1..14] do
		j := i*molt mod m;
		for pol in pols do
			if Terms(pol)[1] eq T[i, j] then
				x := - Kconn!Terms(pol)[2];
				//x := (KtoKconn^(-1))(x);
				T_tau[i, j] := x ;
			end if;
		end for;
	end for;
	return T_tau;
end function;


/* -------------------------------------------------------------------- */

