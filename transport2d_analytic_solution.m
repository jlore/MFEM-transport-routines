clearvars;

pc = phys_const;
e = pc.e;
mp = pc.mp;

n = 3e19;
a = 0.4;
b = 0.64;
chi = 3;

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
rt = roots([A B C D]);
rt = rt(imag(rt)==0);
Tb = rt.^2

r = linspace(a,b,100);


T = -a*qup/(n*chi*e)*log(r/b) + Tb;
Ta = T(1);
Tr  = Tb + (Ta - Tb)*log(r/b)/log(a/b);
Tr2 = Tb + (Ta - Tb)*log(b./r)./log(b/a);

[rk,Tk] = kob_data ;

figure; hold on; box on; grid on;
plot(r,T,'k')
plot(rk,Tk,'o')
% plot(r,Tr,'x-k')
% plot(r,Tr2,'<-')
xlabel('r (m)')
ylabel('T (eV)')
set(gcf,'color','w')

function [r,T] = kob_data 
d=100*[    0.004008739271464   4.347547148147084
   0.004034522203404   4.283419927089014
   0.004122936820048   4.085225478657748
   0.004211355742458   3.892865344317823
   0.004299774664868   3.700505209977897
   0.004377153601056   3.548963745443063
   0.004465572523466   3.356603611103137
   0.004565074489767   3.181767947871515
   0.004638763383759   3.030219307058586
   0.004708775152855   2.896166432241582
   0.004808268507621   2.709662140827281
   0.004951969170709   2.424060625197347
   0.005125185865603   2.132682205700835
   0.005209936274651   1.969486465539512
   0.005460557740334   1.561572466056204
   0.005700117690961   1.165305565921290
   0.005788579671039   1.031288572494760
   0.005873355914688   0.903098716881475
   0.006031855498464   0.670034733185979
   0.006201403679995   0.407820707868070
   0.006367283348164   0.174771076728764
   0.006407830754657   0.116506874874413];
r = d(:,1);
T = d(:,2);
end