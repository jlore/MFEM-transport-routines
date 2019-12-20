clearvars;

%% Domain and geometry
geo.a = 0.4;
geo.b = 0.64;
nr = 50;
r = linspace(geo.a,geo.b,nr);

R0 = 3.9;
beta = 1e-3;

Asurf = 4*pi^2*R0*geo.a;

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

%% 
ntime = 40;
t0 = 0;
t1 = 1;
time = linspace(t0,t1,ntime);

%% get analytic solution
T_analytic = transport2d_analytic_solution(r,geo,plasma);
T0 = icfun(r,geo,plasma,T_analytic);

figure; hold on; box on; grid on;
plot(r,T_analytic,'k','linew',3)
plot(r,T0,'linew',4)
xlabel('r (m)')
ylabel('T (eV)')
legend('analytic','IC')


m = 1;
options = optimset;
sol = pdepe(m,@pdefun,@icfun,@bcfun,r,time,options,geo,plasma,T_analytic);

Tsol = sol(:,:,1);

% plot(r,Tsol(:,:))
plot(r,Tsol(end,:),'o')
legend('analytic','IC','Sol at t1')

figure; hold on; box on; grid on;
ieval = floor(nr/2);
plot(time,Tsol(:,ieval))
title('Evolution at grid midpoint')
xlabel('t (s)')
ylabel('T (eV)')

figure; hold on; box on; grid on;
plot(r,abs(Tsol(end,:)-T_analytic))
set(gca,'yscale','log')
title('T(t1) - T_{analytic}')
xlabel('r (m)')


figure;
surf(r,time,Tsol)
xlabel('r')
ylabel('t')
zlabel('T(r,t)')
view([150 25])

function [c,f,s] = pdefun(r,t,T,dTdr,geo,plasma,T_analytic)
    c = 1.5*plasma.density;
    f = plasma.density*plasma.chi*dTdr;
    s = 0;
end

function T0 = icfun(r,geo,plasma,T_analytic)
Tmax = 500;
Tmin = 5;
Te_exp = 0;

rnorm = (r-geo.a)/(geo.b-geo.a);

% Tmin = T_analytic(1);
% Tmax = T_analytic(end);
% T0 = Tmax - (Tmax-Tmin)*rnorm;

% T0 = 100*sin(r*pi);

rs = (r - 0.5*(geo.a + geo.b)).^2;
T0 = Tmin + (Tmax - Tmin)*(0.5 + 0.5*cos(pi*rnorm)) + 0.5*Te_exp*exp(-400*rs);
end


function [pl,ql,pr,qr] = bcfun(rl,Tl,rr,Tr,t,geo,plasma,T_analytic)
e = 1.602176487000000e-19;
amu = 1.660539040000000e-27;

pl = plasma.qup;
ql = e; 
M = amu*plasma.mi_amu;

pr = plasma.density*geo.beta*sqrt(2*e*Tr/M)*(plasma.gamma*Tr + plasma.Erc)*e;
% pr = plasma.density*geo.beta*sqrt(e*2/M)*plasma.gamma*Tr^1.5*e;
qr = e;     

end

% % Flux IB + const T outer
% function [pl,ql,pr,qr] = bcfun(rl,Tl,rr,Tr,t,geo,plasma,T_analytic)
% Tb = T_analytic(end);
% 
% e = 1.602176487000000e-19;
% 
% pl = plasma.qup;
% ql = e; 
% 
% pr = Tr - Tb;
% qr = 0;
% 
% end

% Dirichlet
% function [pl,ql,pr,qr] = bcfun(rl,Tl,rr,Tr,t,geo,plasma,T_analytic)
% Ta = T_analytic(1);
% Tb = T_analytic(end);
% pl = Tl - Ta;
% ql = 0; 
% pr = Tr - Tb;
% qr = 0;
% end