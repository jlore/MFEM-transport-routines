clearvars;

%% Domain and geometry
a = 0.4;
b = 0.64;
nr = 100;
r = linspace(a,b,nr);

R0 = 3.9;
beta = 1e-3;

Asurf = 4*pi^2*R0*a;

geo.beta = beta;

%% Plasma quantities
Psol = 2e6;
density = 3e19;
chi = 3;
gamma = 5.5;
mi_amu = 1;
Erc = 31;

qup = Psol/Asurf;

plasma.qup = qup;
plasma.density = density;
plasma.chi = chi;
plasma.gamma = gamma;
plasma.mi_amu = mi_amu;
plasma.Erc = Erc;

%% Solve
transport2d_analytic_solution(r,geo,plasma)

