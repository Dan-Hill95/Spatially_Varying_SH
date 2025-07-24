% Copyright 2016 by Daniele Avitabile
% See https://www.maths.nottingham.ac.uk/personal/pmzda/
%
% If you use this code, please cite
% Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on
% Mathematical Neuroscience, Antibes Juan-les-Pins, 2016

function F = SolutionMeasures_Patch(step,uu,pp,mesh_params)

  %% Rename parameters
     r = mesh_params.r;
     N = mesh_params.N;
     n = floor(length(uu)/N);
     u = uu(1:n*N);
  %% Compute quadrature weights and l2norm
      Lr = max(r);
    hr = abs(r(2)-r(1));   
  w = ones(n*N,1); W = ones(N,1); 
l2Norm = sqrt(sum( hr *w.*repmat(r,n,1).* u.^2)/(2*Lr));
v = zeros(n,1);
% q = u(15);
for j=1:n
% v(j) = max(abs(u(1+(j-1)*N:j*N)));
% v(j) = sign(max(u(1+(j-1)*N:j*N)) + min(u(1+(j-1)*N:j*N)))*max(abs(u(1+(j-1)*N:j*N)));
v(j) = sign(max(u(1+(j-1)*N:j*N)) + min(u(1+(j-1)*N:j*N)))*sqrt(sum( hr * W.*r.*u(1+(j-1)*N:j*N).^2)/(2*Lr));
end

  %% Allocate
     F = zeros(1,n+2);

  %% Assign branch variables
     F = [l2Norm v'];
%      F = [l2Norm v' q];
  
end
