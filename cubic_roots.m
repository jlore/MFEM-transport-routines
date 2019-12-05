function z = cubic_roots(A,B,C,D)
% y = A*x^3 + B*x^2 + C*x + D


d0 = B^2 - 3*A*C;
d1 = 2*B^3 - 9*A*B*C + 27*A^2*D;

c = ((d1 + sqrt(d1^2 - 4*d0^3))/2)^(1/3);


xi = (-1 + sqrt(-3))/2;

z1 = -1/(3*A)*(B + xi  *c + d0/(xi  *c));
z2 = -1/(3*A)*(B + xi^2*c + d0/(xi^2*c));
z3 = -1/(3*A)*(B + xi^3*c + d0/(xi^3*c));

z = [z1; z2; z3];
