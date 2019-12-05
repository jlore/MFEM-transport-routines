function T = transport2d_analytic_solution(r,geo,plasma)

plotit = 0;

mi_amu = plasma.mi_amu;
qup = plasma.qup;
beta = geo.beta;
n_i = plasma.density;
gamma = plasma.gamma;
Erc = plasma.Erc;
chi = plasma.chi;

pc = phys_const;

a = r(1);
b = r(end);

M = mi_amu*pc.mp;

A = beta*n_i*gamma*sqrt(2*pc.e/M)*pc.e;
B = 0;
C = beta*n_i*sqrt(2*pc.e/M)*Erc*pc.e;
D = -a*qup/b;
rt = roots([A B C D]);
rt = rt(imag(rt)==0);
Tb = rt.^2;

T = Tb - a*qup/(n_i*chi*pc.e)*log(r/b);
% Ta = T(1);
% Tr  = Tb + (Ta - Tb)*log(r/b)/log(a/b);
% Tr2 = Tb + (Ta - Tb)*log(b./r)./log(b/a);

if plotit
[rk,Tk] = kob_data;

figure; hold on; box on; grid on;
plot(r,T,'k')
plot(rk,Tk,'o')
% plot(r,Tr,'x-k')
% plot(r,Tr2,'<-')
xlabel('r (m)')
ylabel('T (eV)')
set(gcf,'color','w')
end


end


