/* This routine computes the matrix rho(tau)^vee acting on the
rank 1 Q([alpha])-module of Hodge cycles in H(X_m^n)[alpha] */

/* --------------------- global set-up ------------------------- */

m := 15;																	// <--
C := ComplexField(400);														// <--

coprimes_to_m := [ u : u in [1..(m-1)] | GCD(u, m) eq 1 ];

P<x> := PolynomialRing(Rationals());
z_m := Roots(ChangeRing(CyclotomicPolynomial(m), C))[1][1];

K := SplittingField(CyclotomicPolynomial(m));
V, K_coordinates := VectorSpace(K, Rationals());

PolinomioDiKconn := x^16 - 5*x^15 + 14*x^14 - 30*x^13 + 57*x^12 - 100*x^11 + 157*x^10- 215*x^9 + 250*x^8 - 240*x^7 + 183*x^6- 110*x^5 + 57*x^4 - 30*x^3 + 16*x^2 - 5*x + 1;											// <--
Kconn := SplittingField(PolinomioDiKconn);
_, KtoKconn := IsSubfield(K,Kconn);											// Fix an embedding of K in Kconn
iota := hom< Kconn -> C | Roots(ChangeRing(PolinomioDiKconn, C))[11][1] >;	// Fix an embedding of Kconn in C
GaloisKconn, AutKconn, GaloisMap := AutomorphismGroup(Kconn);				// Compute the Galois group of Kconn/Q



/* ------------------------------------------------------------ */


function compute_alpha(Num, Den)											// Given an equation f, computes the corresp. char.
	alpha_dualized := Num cat [ -b mod m : b in Den];						// Numerator's pedices cat -Denominator's ones
	alpha_tripled := &cat[[i, i, -2*i mod m] : i in alpha_dualized];    	// triplicate entries - according to Shioda
    return alpha_tripled;
end function; 


/* ------------------------------------------------------------ */

function P_sigma(alpha)														// Computes P(sigma, omega_alpha)
	gamma_factor := &*[Gamma(C!1-a/m): a in alpha];							// we add one to make the computation happens (we rescale back later)
    pi_factor := ( 2*C.1*Pi(C) )^(&+alpha/m);								// usual 2*pi*i facor
    
    algebraic_number := (gamma_factor / pi_factor);							// this number is algebraic and in Kconn by theory
    min_poly, _ := MinimalPolynomial(algebraic_number, 30);					// numerically compute its minimal polynomial
    rr := Roots(min_poly, Kconn);											// take its roots in Kconn
    residues := [Norm(iota(x[1]) - algebraic_number): x in rr];				// let's search for the root closest to algebraic_number, under the fixed embedding iota of Kconn into C
    _, minimal_index := Min( residues );									// the index of the closest root in the array rr
    true_gamma := rr[minimal_index][1];										// this is the representative of algebraic_number in Kconn
    
    true_xi := &*[KtoKconn(1-(K.1)^a): a in alpha];							// the rescaling factor xi(alpha)
    rescaling_factor := &*[Kconn!(-a/m): a in alpha];						// the rescaling factor we need to actually compute Gamma(-a/m)
    
    return (true_gamma * true_xi / rescaling_factor);						// return the number as element of Kconn
end function;


function compute_gamma_coordinates(alpha)									// Computes the coordinate of gamma in the omega_u*alpha basis
    gamma := [];
    for u in [ u : u in [1..(m-1)] | GCD(u, m) eq 1 ] do					// notice gamma[u] will have some undefined entries
    	u_alpha := [(a*u mod m): a in alpha];								// computes the other characters in the orbit
    	gamma[u] := P_sigma(u_alpha);										// this is a number in Kconn
    end for;
    return gamma;
end function;


function compute_zetagamma_coordinates(alpha)								// Computes the coordinates of all z^i*gamma in the basis omega_u*alpha
	/* in order to compute the omega_alpha coordinates of K.1^i*gamma, we pick the charachter
	zeta_underscore = (1 everywhere but K.1^j in an index a sucht that (a,15)=1 and a*j=1 mod 15) */
	gamma := compute_gamma_coordinates(alpha);								// the coordinates of gamma in the basis omega_alpha^\vee
	zetagamma_coord := [ ];													// the coordinates of z^i * gamma in the omega_alpha basis, for i from 0 to 7.
	for j in [1..EulerPhi(m)] do
    	zetagamma_coord[j] := [ K.1^(-(j-1)*u) * gamma[u] : u in coprimes_to_m ];        // the minus is there because of the dual (??)
	end for;
	zetagamma_coord := Matrix(Kconn, EulerPhi(m), EulerPhi(m), zetagamma_coord);		// put everything in a base-change matrix
	
	return zetagamma_coord;
end function;



/* ------------------------------------------------------------ */


function lambda(tau, alpha)        											// we compute directly lambda^-1(tau) !!
    numerello := P_sigma(alpha);											// compute P(sigma, omega_alpha)
    rapporto := GaloisMap(tau)(numerello) / numerello;						// compute Galois action, we should get a number in K
	return (KtoKconn^(-1))(rapporto);										// ... this works!
end function;


function ActionOnZeta(tau, i)														// There might have been a better way to find out where tau send K.1^i
    for j in [0..(m-1)] do
        y := (GaloisMap(tau)(K.1^i)) / K.1^j;
        if y eq 1 then
            return j;
        end if;
    end for;
end function;


function tau_action_on_zetagamma(alpha, tau)								// Computes the matrix of tau acting on z^i*gamma
	tau_acts_on_zetagamma := [];
	lambda_tau := lambda(tau, alpha);										// the lambda factor as an element of K
	for i in [1..EulerPhi(m)] do
		eigenvector_in_Qalpha := K.1^ActionOnZeta(tau, i-1) * lambda_tau;			// multiply by the corresponding tau(z^i)
		eigen_coord := K_coordinates(eigenvector_in_Qalpha);					// coordinates in terms of the z^i basis of K
		tau_acts_on_zetagamma := Append(tau_acts_on_zetagamma, eigen_coord);// get a matrix...
	end for;
	tau_acts_on_zetagamma := Matrix(Kconn, EulerPhi(m), EulerPhi(m), tau_acts_on_zetagamma);
	return tau_acts_on_zetagamma;
end function;


/* ------------------------------------------------------------ */


function compute_matrix(f, tau)												// Computes the matrix rho(tau)^vee
	alpha_f := compute_alpha(f[1], f[2]);									// this is the character alpha corresp. to the equation f
	
	zetagamma_coord := compute_zetagamma_coordinates(alpha_f);				// the base-change matrix in bitween omegas and z^i*gammas
	tau_acts_on_zetagamma := tau_action_on_zetagamma(alpha_f, tau);			// the matrix rho(tau)^vee in z^i basis
	
	tau_per_omega := zetagamma_coord^-1 * tau_acts_on_zetagamma * zetagamma_coord; 	// base-changing
	
	return tau_per_omega;
end function;


/* ------------------------------------------------------------ */
