clearvars;

pc = phys_const;
e = pc.e;
mp = pc.mp;

n = 3e19;
a = 0.4;
b = 0.64;
chi = 1;

% Tb = 11.0623;
P = 2e6;
R0 = 3.9;
qup = P/(4*pi^2*R0*a);
beta = 1e-3;
gamma = 5.5;
M = 1*mp;

Erc = 31;

A = beta*n*gamma*sqrt(2*e/M)*e;
B = 0;
C = beta*n*sqrt(2*e/M)*Erc*e;
D = -a*qup/b;
% p = [3 -2 -4 1];
p = [A B C D];

Tb_test = linspace(1,50,1000);
ftest = A*Tb_test.^(1.5) + C*sqrt(Tb_test) + D;
% 
% figure; hold on; grid on;
% plot(Tb_test,ftest)
% 
% 
% figure; hold on; grid on;
% plot(sqrt(Tb_test),ftest)

% asdfasd

r = roots(p)
% r = r(imag(r)==0)
% Tb = r.^2


z = cubic_roots(A,B,C,D)


asdfdsf
Tb= 11.0623
% polyval(p,Tb)



b/(a*qup)*beta*n*sqrt(2*e*Tb/M)*(gamma*Tb+Erc)*e

beta*n*sqrt(2*e*Tb/M)*(gamma*Tb+Erc)*e - (a*qup)/b



