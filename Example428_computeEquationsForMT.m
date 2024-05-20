/* The script computes equations for MTJ15 */

SetLogFile("Examples428and4210.out");

m := 15;
load "FermatJacobians1_computeMTandKconn.m";

MTmatrix := Matrix(Rationals(), [[b[i] : i in [1..Ncols(b)] ]: b in B ]);


CambioBase := Matrix(Rationals(),9,9,[
 [ 1, 1, 0, 0,  0, -1, 0, 0, 0 ],
 [ 0, 1, 0, 0,  0,  0, 0, 0, 0 ],
 [ 0, 0, 1, 0,  0,  0, 0, 0, 0 ],
 [ 0, 0, 0, 1,  1, -1, 0, 0, 0 ],
 [ 0, 0, 0, 0,  1,  0, 0, 0, 0 ],
 [ 0, 0, 0, 0, -1,  1, 0, 0, 0 ],
 [ 0, 0, 0, 0,  0,  0, 1, 0, 0 ],
 [ 0, 0, 0, 0,  0,  0, 0, 1, 0 ],
 [ 0, 0, 0, 0,  0,  0, 0, 0, 1 ]]);

"Is the base change invertible?";
IsInvertible(CambioBase);
 
BB := CambioBase*MTmatrix;

PP := PolynomialRing(Rationals(), m-1);
RR := RationalFunctionField(PP);
AssignNames(~PP, ["x_{" cat Sprint(j) cat "}": j in [1..(m-1)]]);
all_equations := [];
for x in [1..Nrows(BB)] do
	v := BB[x];
	eqn := RR!1;
		for i in [1..(m-1)] do
		    j := OrderCharacters[i];
		    if not(v[i] eq 0) then
		        eqn := eqn * (RR!(PP.j))^(Integers()!v[i]);
		    end if;
		end for;
	Append(~all_equations, eqn);
end for;

"Equations for the MT group of J15:";
all_equations;

/* ------------------------------------------------------------ */

m := 10;
load "FermatJacobians1_computeMTandKconn.m";
"Equations for the MT group of J10:";
print_MT_equations();

