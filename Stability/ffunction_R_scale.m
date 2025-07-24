function F = ffunction_R_scale(R,epsilon,k,mu)

%    psi=atan2(sqrt(epsilon^2+mu)/sqrt(epsilon^2-mu),1);
%sqrt(epsilon^2-mu)
%sqrt(epsilon^2+mu)
    psi=angle(sqrt(1-mu)+1i*sqrt(1+mu));
    alphap=sqrt(1+epsilon*sqrt(1-mu));
    alpham=sqrt(1-epsilon*sqrt(1-mu));
    beta=sqrt(1+1i*epsilon*sqrt(1+mu));

    u1=besselj(k,alphap*R);
    u2=besselj(k,alpham*R);
    u3=besselh(k,beta*R);

    du1=-alphap*besselj(k+1,alphap*R)+2.0*k/R*besselj(k,alphap*R);
    du2=-alpham*besselj(k+1,alpham*R)+2.0*k/R*besselj(k,alpham*R);
    du3=-beta*besselh(k+1,beta*R)+2.0*k/R*besselh(k,beta*R);

    Wu1u3=u1*du3-u3*du1;
    Wu1u3b=u1*conj(du3)-conj(u3)*du1;
    Wu2u3=u2*du3-u3*du2;
    Wu2u3b=u2*conj(du3)-conj(u3)*du2;
    
    exp_phi=1i*(Wu1u3b*u3)/(Wu1u3*conj(u3));
%    angle(exp_phi)
%    abs(exp_phi)

    F=(exp(4*1i*psi)*Wu2u3b/Wu2u3-Wu1u3b/Wu1u3);
%    G=imag(exp(4*1i*psi)*Wu2u3b/Wu2u3-Wu1u3b/Wu1u3);
%    F=real(exp(4*1i*psi)*Wu2u3b*Wu1u3-Wu1u3b*Wu2u3);

end