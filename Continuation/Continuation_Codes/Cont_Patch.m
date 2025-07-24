% Written by Dan Hill (2021) - University of Surrey, 
% adapted from codes by David Lloyd - University of Surrey and
% Daniele Avitabile - Vrije Universiteit Amsterdam;
% see Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on
% Mathematical Neuroscience, Antibes Juan-les-Pins, 2016.
%
%% Inputs
%           p=[mu, nu, kappa, k, r0, m, N+1]   - Initial parameters of system; 
%                                         mu is the bifurcation parameter, 
%                                         nu is the quadratic coefficient,
%                                         kappa is the cubic coefficient,
%                                         m is the dihedral symmetry,
%                                         N+1 is the dimension of the ODE system
%           Dir                         - Direction of parameter continuation;
%                                         must be either 'pl' (plus), 'mn' (minus),
%                                         'sp' (small plus), or 'sm' (small minus).
%
%% Purpose
% Solves a N'th order Galerkin system for the stationary 2D 2-3 Swift-Hohenberg
% equation via finite-difference methods and continues solutions in mu-parameter
% space. For U(r,theta) such that
%
% U(r,theta) = u[0](r) + 2*sum_{i=1}^{N} u[i](r)*cos(m*i*theta), 
% then each u[i](r) satisfies 
%
% 0= F[i](u[i]):= -(1+d^2_r + 1/r*d_r - (2*m*i/r)^2)^2 u[i] - mu*u[i] + f[i](U).
%
% After solving an algebraic matching condition solutions are then 
% continued in mu-parameter space.
%
%% Outputs
%           branch- [Step, -1, mu, EucNorm(u), L2Norm(u[0]), max(|u[0]|), ..., max(|u[N]|)]
%
% All data is stored in a folder named as "D[m]_Patch_[Dir]" 
%
function branch = Cont_Patch(p,Dir)

close all, clc;

k = p(4);
r0 = p(5);
m = round(p(6));                   % Dihedral lattice
n = p(7);                          % ODE dimension

% Setting up mesh parameters, collected as mesh_params, including
% finite-difference matrices
SetupDiffMats_Patch;

r0= mesh_params.r0;

F = @(eps) matching(k,eps,r0);

eps0 = sqrt(p(1));

options = optimset('Jacobian','off','Display','iter','TolFun',1e-5,'DerivativeCheck','off');

Eps=fsolve(F, eps0,options);

[fig1,fig2] = Eigenfunction_plots(k,0,r0);
pause;close(fig1);close(fig2);

% if p(2)==0
% p(1) = Eps^2 - sign(p(3))*1e-04;
% else
% p(1) = Eps^2 - sign(19*p(2)^2/18 - 3*p(3)/4)*1e-04;
% end
% p(1) = Eps^2-1e-03;

%% Initial guess
r=mesh_params.r;
u0 = linear_eig(k,Eps,r,r0);
p0 = p;

uu0 = -1*[zeros((k/m)*N,1);u0;zeros((n-1-(k/m))*N,1)];

myproblemHandle = @(u) Equation_Patch(u,p0,mesh_params);
    
% Fsolve options
options = optimset('Jacobian','on','Display','iter','MaxIter',50,'TolFun',1e-7,'DerivativeCheck','off');

%Solving the "k"'th-order Galerkin truncation of the Swift-Hohenberg
%equation close to the initial guess
iter=1;
[u_out,fval,exitflag,output,jacobian] = fsolve(myproblemHandle,uu0,options);
while norm(u_out,2)<1e-06 && iter<20
[u_out,fval,exitflag,output,jacobian] = fsolve(myproblemHandle,(-1)^(iter)*(1+iter*100)*uu0,options);
uu0=u_out;
iter=iter+1;
end

u1 = u_out;
if (k/m)*N>2
u1(i<(k/m)*N+1)=0;
end

%% Plotting the surface of the found solution

hfig1=figure;
t = 0:0.01:2*pi;                % Angular variable
t = t';
[R,T]=meshgrid(r,t);
[U,T0]=meshgrid(u1,t);
UU(:,:,1)=U(:,1:N)/2;
  for i=1:n-1
      UU(:,:,i+1)=U(:,1+i*N:(i+1)*N);
      UU(:,:,i+1)=UU(:,:,i+1).*cos(m*i.*T);
  end
z=surf(R.*cos(T), R.*sin(T),sum(UU,3));
view(0,90);
z.FaceColor='flat';
z.FaceAlpha=1;
z.EdgeColor='none';
  load("Cool_Warm.mat","mat")
newmap = mat(round(linspace(1, 256, 11)), :);
colormap(hfig1,newmap);
 clim([-max(abs(u1)) max(abs(u1))]);
axis([-0.5*L*sqrt(2) 0.5*L*sqrt(2) -0.5*L*sqrt(2) 0.5*L*sqrt(2) -1 2]);
  
pause;
close(hfig1);

%% Define handle to right-hand side and time output function
prob     = @(u,p) Equation_Patch(u,p,mesh_params);
plotSol  = @(u,p,parent) PlotSurface_Patch(u,p,parent,mesh_params);
solMeas  = @(step,u,p) SolutionMeasures_Patch(step,u,p,mesh_params);
compSpec = [];%@(u,p) ComputeEigenvalues_Patch(u,p,prob);
plotSpec = [];%@(d,p,parentHandle) PlotSpectrum_SH(d,p,parentHandle);

%% Assign problem 
stepperPars.iContPar      = 1;
% Define the direction and step distance for continuation
if Dir == 'pl'
if n == 1
stepperPars.s0            = 0.05;
stepperPars.sMin          = 1e-4;
stepperPars.sMax          = 0.1;
else
stepperPars.s0            = 0.01;
stepperPars.sMin          = 1e-9;
stepperPars.sMax          = 0.05;
end
elseif Dir == 'mn'
if n == 1
stepperPars.s0            = -0.05;
stepperPars.sMin          = 1e-4;
stepperPars.sMax          = 0.1;
else
stepperPars.s0            = -0.01;
stepperPars.sMin          = 1e-9;
stepperPars.sMax          = 0.05;
end
elseif Dir== 'sp' 
stepperPars.s0            = 1e-4;
stepperPars.sMin          = 1e-9;
stepperPars.sMax          = 1e-3;
elseif Dir== 'sm' 
stepperPars.s0            = -1e-4;
stepperPars.sMin          = 1e-9;
stepperPars.sMax          = 1e-3;
else
    error('Final argument must be "pl" (forwards), "mn" (backwards), "sp" (small values forwards), or "sm" (small values backwards).');
end

stepperPars.pMin          = 0;
stepperPars.pMax          = 5;
stepperPars.maxSteps      = 2000;
stepperPars.nPrint        = 1;
stepperPars.nSaveSol      = 1;
stepperPars.finDiffEps    = 1e-7;
stepperPars.fsolveOptions = optimset('Display','off',...
                                     'DerivativeCheck','off',...
                                     'TolFun',1e-7,...
                                     'Jacobian','on',...
                                     'MaxIter',10);
stepperPars.optNonlinIter = 10;
stepperPars.dataFolder    = ['D',num2str(k),'_Bif_',num2str(r0),'_', Dir];
stepperPars.PlotSolution  = plotSol;
stepperPars.BranchVariables = solMeas;
% stepperPars.PlotBranchVariableId = [];%2;
stepperPars.ComputeEigenvalues = compSpec;
stepperPars.PlotSpectrum = plotSpec;      
stepperPars.PlotBranchVariableId = 1;


branch = SecantContinuation(prob,u1,p,stepperPars,'Branch');
end