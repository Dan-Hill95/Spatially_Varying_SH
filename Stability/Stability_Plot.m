function fig1 = Stability_Plot(R,k)

% clear
clc

%options = optimoptions('fsolve','Display','iter','FunctionTolerance',1E-12,'StepTolerance',1E-12);

mu=0;
% R=5.0;
%epsilon=1.12;
% k=6;

num=100;

epsilon=linspace(0,1,num);

for m=num:-1:1

        m
        if (m==num)
%         mu2=(num-1)/num*epsilon(m)^2/2;
%         F2=ffunction_R(R,epsilon(m),k,mu2);
        mu2=(num-1)/num/2;
        F2=ffunction_R_scale(R,epsilon(m),k,mu2);
%         mu1=(num-2)/num*epsilon(m)^2/2;
%        mu1=-epsilon(m)^2;
%         F1=ffunction_R(R,epsilon(m),k,mu1);
        mu1=(num-2)/num/2;
        F1=ffunction_R_scale(R,epsilon(m),k,mu1);
        else
        mu2=mu_val+2*(epsilon(2)-epsilon(1));
%         F2=ffunction_R(R,epsilon(m),k,mu2);
        F2=ffunction_R_scale(R,epsilon(m),k,mu2);
        mu1=mu_val+1*(epsilon(2)-epsilon(1));
%        mu1=-epsilon(m)^2;
%         F1=ffunction_R(R,epsilon(m),k,mu1);
F1=ffunction_R_scale(R,epsilon(m),k,mu1);
        end
        while (real(F1)*real(F2)>0 || imag(F1)*imag(F2)>0 )
            mu2=mu1;
            F2=F1;
%             mu1=mu1-1/num*epsilon(m)^2/1000;
            mu1=mu1-1/num/1000;
%             if (mu1<=-epsilon(m)^2)
            if (mu1<=-1)
%            mu1=-epsilon(m)^2;
                break
            end
%             F1=ffunction_R(R,epsilon(m),k,mu1);
            F1=ffunction_R_scale(R,epsilon(m),k,mu1);
            end

%            F1
%            F2
%     if (mu1<=-epsilon(m)^2)
    if (mu1<=-1)
        mu_val=NaN;
    else
        mu_val=real(mu1-F1*(mu2-mu1)/(F2-F1));        
    end

    mu_val=real(mu1-F1*(mu2-mu1)/(F2-F1));

%    if (mu_val>0)
%        muplot(m)=1.0;
%    else
%        muplot(m)=-1.0;
%    end
    muplot(m)=mu_val;

%    figure(3);plot(epsilon(1:m),muplot(1:m),'b.-',epsilon(1:m),epsilon(1:m)*0,'k-')
end

for m=num:-1:1

        m
        if (m==num)
%         mu2=(num-1)/num*epsilon(m)^2/1;
%         F2=ffunction_R(R,epsilon(m),k,mu2);
mu2=(num-1)/num/1;        
F2=ffunction_R_scale(R,epsilon(m),k,mu2);
%         mu1=(num-2)/num*epsilon(m)^2/1;
%        mu1=-epsilon(m)^2;
%         F1=ffunction_R(R,epsilon(m),k,mu1);
mu1=(num-2)/num/1;
        F1=ffunction_R_scale(R,epsilon(m),k,mu1);
        else
        mu2=mu_val+0.5*(epsilon(2)-epsilon(1));
%         F2=ffunction_R(R,epsilon(m),k,mu2);
        F2=ffunction_R_scale(R,epsilon(m),k,mu2);
        mu1=mu_val+0.1*(epsilon(2)-epsilon(1));
%        mu1=-epsilon(m)^2;
%         F1=ffunction_R(R,epsilon(m),k,mu1);
        F1=ffunction_R_scale(R,epsilon(m),k,mu1);
        end
        while (real(F1)*real(F2)>0 || imag(F1)*imag(F2)>0 )
            mu2=mu1;
            F2=F1;
%             mu1=mu1-1/num*epsilon(m)^2/10;
mu1=mu1-1/num/10;
%             if (mu1<=-epsilon(m)^2)
if (mu1<=-1)
%            mu1=-epsilon(m)^2;
                break
            end
%             F1=ffunction_R(R,epsilon(m),k,mu1);
            F1=ffunction_R_scale(R,epsilon(m),k,mu1);
            end

%            F1
%            F2
%     if (mu1<=-epsilon(m)^2)
if (mu1<=-1)
        mu_val=NaN;
    else
        mu_val=real(mu1-F1*(mu2-mu1)/(F2-F1));        
    end

    mu_val=real(mu1-F1*(mu2-mu1)/(F2-F1));

%    if (mu_val>0)
%        muplot(m)=1.0;
%    else
%        muplot(m)=-1.0;
%    end
    muplot2(m)=mu_val;

%    figure(3);plot(epsilon(1:m),muplot(1:m),'b.-',epsilon(1:m),epsilon(1:m)*0,'k-')
end


    fig1=figure(3);plot(epsilon,epsilon*0,'k-',epsilon,epsilon.^2.*muplot,'b-',epsilon, epsilon.^2.*muplot2,'b-',epsilon,-epsilon.^2,'r.','Linewidth',2);xline(0);
axis([0 1 -1 1])

end