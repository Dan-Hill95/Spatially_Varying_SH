function [fig1,fig2] = Eigenfunction_plots(k,j,R)

%% Initialisation
% clear all; close all; clc;

% k = 0;                      % Bessel index

% R = 5;                     % Point of discontinuity

r = linspace(0,20,2000);   % Radial coordinate

% j = 0;                      % Number of root for k'th matching equation

%% Finding critical epsilon

eps0=fsolve(@(eps) matching(k,eps,R), (2*j+1)*pi/(2*R));

sprintf('eps_%d,%d = %d',k,j,eps0)

%% Finding linear eigenfunction

z = linear_eig(k,eps0,r,R);   % z is linear eigenfunction
if max(abs(z)) == max(z)
 
z = z/(max(abs(z)));          % normalising z, so that max(|z|) = 1

else

z = -z/(max(abs(z)));          % normalising z, so that max(|z|) = 1

end

%% Plots

if nargout>1

% Radial plot

fig1=figure;
plot(r,z,r,0.*z,'k');xline(R);xlabel('r');ylabel(['v_',num2str(k),'_,_',num2str(j)]);
axis([r(1) r(end) -1.1  1.1]);
pbaspect([1 1 1]);
set(gca,'Color','none','XColor','none','YColor','none','ZColor','none');
  grid off

end
% Planar plot

t = linspace(0,2*pi,4000);

[rr,tt] = meshgrid(r,t);

[zz,~]=meshgrid(z,t);
zz= zz.*cos(k.*tt);

fig2 = figure;
p = surf(rr.*cos(tt), rr.*sin(tt),zz);
view(0,90); 
p.FaceColor='flat';
p.FaceAlpha=1;
p.EdgeColor='none';
load("Cool_Warm.mat","mat")
newmap = mat(round(linspace(1, 256, 11)), :);
colormap(fig2,newmap);
clim([-1.1 1.1]);
set(gca,'Color','none','XColor','none','YColor','none','ZColor','none');
  grid off
axis off
pbaspect([1 1 1]);
% axis([-50*sqrt(2) 50*sqrt(2) -50*sqrt(2) 50*sqrt(2) -2 3]);
end