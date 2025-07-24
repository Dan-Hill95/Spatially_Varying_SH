function z = linear_eig(k,eps,r,R)
%k is Bessel index                      k natural number
%eps is linear pumping coefficient      0<eps<1
%R is point of discontinuity            0<R<L
%r is radius                            0<r<L

ap=sqrt(1+eps);
am=sqrt(1-eps);
b=sqrt(1+1i*eps);

z = 0.*r;

r(r==R)=R+1e-04;
r1=r(r<R); i1=length(r1);
r2=r(r>R); i2=length(r2);

u1 = besselj(k,ap.*r1)./besselj(k,ap.*R);
u2 = besselj(k,am.*r1)./besselj(k,am.*R);
u3 = besselh(k,b.*r2)./besselh(k,b.*R);

ephi=(((conj(b.*besselh(k-1,b.*R)).*besselj(k,am.*R) - am.*besselj(k-1,am.*R).*conj(besselh(k,b.*R))).*besselh(k,b.*R))./conj((conj(b.*besselh(k-1,b.*R)).*besselj(k,am.*R) - am.*besselj(k-1,am.*R).*conj(besselh(k,b.*R))).*besselh(k,b.*R))).^(1/2);

z1=1/sqrt(2).*(real(ephi).*u1 + imag(ephi).*u2);
z2= real(((1-1i)/sqrt(2)).*ephi.*u3);


z(1:i1) = z1; 
z(i1+1:i1+i2) = z2;
end