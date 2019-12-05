clearvars;

%% Domain and geometry
geo.a = 0;
geo.b = 1;
nr = 105;
r = linspace(geo.a,geo.b,nr);

%% 
ntime = 10;
t0 = 0;
t1 = 10;
time = linspace(t0,t1,ntime);



T0 = icfun(r);

figure; hold on; box on; grid on;
plot(r,T0,'k','linew',2)
xlabel('r (m)')
ylabel('T (eV)')


m = 0;
options = optimset;
sol = pdepe(m,@pdefun,@icfun,@bcfun,r,time,options);

Tsol = sol(:,:,1);

plot(r,Tsol(:,:))
plot(r,Tsol(end,:),'linew',3)

figure; hold on; box on; grid on;
ieval = floor(nr/2);
plot(time,Tsol(:,ieval))
xlabel('t (s)')
ylabel('T (eV)')



function [c,f,s] = pdefun(r,t,T,dTdr)
c = 1;
f = dTdr;
s = 0;
end

function T0 = icfun(r)
T0 = 100*sin(r*pi);
end

function [pl,ql,pr,qr] = bcfun(rl,Tl,rr,Tr,t)

pl = 100;
ql = 1; 

pr = Tr;
qr = 1; 

end