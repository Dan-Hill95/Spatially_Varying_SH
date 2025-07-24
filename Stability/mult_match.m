%Implicit Plots of j'th root of k-index matching condition, for all 0 < k < kmax
%Example: run   >> mult_match(20,0);

%Note: currently doesn't really work for j>0, maybe need to include
%Jacobian in matching.m
function fig=mult_match(kmax)
close all
j=0
k=0:1:kmax;  %consider all indices up to kmax
rmax=20;     %fix a maximum radius
rmin = fsolve(@(r) matching(0,1,r),1); %fix the minimum radius to be when the zero'th mode solves eps=1

r=linspace(rmin,rmax,500);
eps0 = (2.*j+1).*pi./(2.*r);   %initial guess close to eps = (2*j+1)*pi/(2*R)
for q=1:kmax+1
F = @(eps) matching(k(q),eps,r);
options = optimset('Jacobian','off','Display','iter','TolFun',1e-5,'DerivativeCheck','off');
Eps=fsolve(F, eps0,options);
fig=figure(1);
plot(Eps,r)
hold on
end

axis([0 1 0 rmax]);
pbaspect([1 1 1]);
end
