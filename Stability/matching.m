%k-index matching function F(eps,r)
function F = matching(k,eps,r)
%k is bessel index
%eps is linear pumping coefficient
%r is the discontiuity point
%Assume that k>-1, eps>=0, and r>0

ap=sqrt(1+eps);
am=sqrt(1-eps);
b=sqrt(1+1i*eps);

f1 = b.*besselh(k-1,b.*r).*besselj(k,ap.*r)-ap.*besselj(k-1,ap.*r).*besselh(k,b.*r);
f2 = conj(b.*besselh(k-1,b.*r)).*besselj(k,am.*r)-am.*besselj(k-1,am.*r).*conj(besselh(k,b.*r));

F = real(f1.*f2);
end