function [e,rate,flag] = Generate_beamforming_e(N, M, G_U, G_E, ...
                    f_ini, e_ini,  noise_maxpower, trans_maxpower,miu_bs,miu_u )
     
cvx_solver mosek
cvx_save_prefs

cvx_begin quiet
    variable E(M+1,M+1) hermitian
    variable r(4,1)    
    variable p(4,1) 
               
    T_u = G_U*f_ini*f_ini'*G_U';   
    D_u = (1+miu_u)*miu_bs*G_U*diag(diag( f_ini*f_ini' ))*G_U';
    temp=trace(((1+miu_u)*T_u+D_u)*E)+noise_maxpower-r(1);
    constraint(1) = -temp;
    constraint(2)=trace((miu_u*T_u+D_u)*E)+noise_maxpower-r(2);
    r_ini(2)=trace(e_ini'*(miu_u*T_u+D_u)*e_ini)+noise_maxpower;
    
    T_e = G_E*f_ini*f_ini'*G_E';  
    D_e = miu_bs*G_E*diag(diag( f_ini*f_ini' ))*G_E';
    constraint(3)=trace((T_e+D_e)*E)+noise_maxpower-r(3);
    r_ini(3)=trace(e_ini'*(T_e+D_e)*e_ini)+noise_maxpower;
    temp=trace(D_e*E)+noise_maxpower-r(4);
    constraint(4) = -temp;
    
    
    constraint(5) = p(1)-log(r(1))/log(2);
    constraint(6) = log2(r_ini(2)) + (r(2)-r_ini(2)) / (r_ini(2)*log(2)) -p(2);
    constraint(7) = log2(r_ini(3)) + (r(3)-r_ini(3)) / (r_ini(3)*log(2)) -p(3);
    constraint(8) = p(4)-log(r(4))/log(2);
      
    maximize p(1) - p(2) - p(3) + p(4)
  
    subject to
    
         real(constraint)<=0;
         diag(E) == 1;
         E       == hermitian_semidefinite(M+1);

cvx_end             
                
                
if cvx_status(1)=='S'  ||  cvx_status(3)=='a'
    flag=1;

[t1,t2]=eig(E);
location=find( abs(diag(t2))>10^(-6));
if size(location,1)==1
    e_hat=t1(:,location)*t2(location,location)^(1/2);
        %%%%%  Obj value  %%%%%
    numi_u=trace(e_hat'*((1+miu_u)*T_u+D_u)*e_hat)+noise_maxpower;
    deno_u=trace(e_hat'*(miu_u*T_u+D_u)*e_hat)+noise_maxpower;
    rate_u=log2( numi_u/deno_u );

    numi_e=trace(e_hat'*(T_e+D_e)*e_hat)+noise_maxpower;
    deno_e=trace(e_hat'*D_e*e_hat)+noise_maxpower;
    rate_e=log2( numi_e/deno_e );
    rate=max(rate_u-rate_e,0);  
else
    for i=1:500
        flag_2=1;
        b1(:,i)=t1*t2^(1/2)*sqrt(1/2)*(randn(M+1,1) + sqrt(-1)*  randn(M+1,1));
        b2(:,i)=exp(1j*angle(b1(:,i)/b1(M+1,i)));
        %%%%%  Obj value  %%%%%
        numi_u=trace(b2(:,i)'*((1+miu_u)*T_u+D_u)*b2(:,i))+noise_maxpower;
        deno_u=trace(b2(:,i)'*(miu_u*T_u+D_u)*b2(:,i))+noise_maxpower;
        rate_u=log2( numi_u/deno_u );
        
        numi_e=trace(b2(:,i)'*(T_e+D_e)*b2(:,i))+noise_maxpower;
        deno_e=trace(b2(:,i)'*D_e*b2(:,i))+noise_maxpower;
        rate_e=log2( numi_e/deno_e );
        Obj(i)=max(rate_u-rate_e,0);   
    end
    [X,locat]=max(Obj);
    e_hat=b2(:,locat);
    rate=Obj(locat);
end
e=exp(1j*angle(e_hat/e_hat(M+1)));

else
    flag=0;
    e=ones(M+1,1);
    rate=0;
end

end