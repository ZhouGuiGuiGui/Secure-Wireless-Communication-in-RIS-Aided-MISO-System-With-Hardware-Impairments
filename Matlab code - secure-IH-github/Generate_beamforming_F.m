function [f,rate,rate_p,flag] = Generate_beamforming_F(N, M, G_U, G_E, ...
                    f_ini, e_ini,  noise_maxpower, trans_maxpower,miu_bs,miu_u )

cvx_solver mosek
cvx_save_prefs

cvx_begin quiet
    variable f(N,1) 
    variable r(4,1)    
    variable p(4,1) 
      
    temp_u = G_U'*e_ini*e_ini'*G_U;     
    A_u    = (1+miu_u)*temp_u+(1+miu_u)*miu_bs*diag(diag( temp_u ));
    temp_e = G_E'*e_ini*e_ini'*G_E; 
    A_e    = miu_bs*diag(diag( temp_e ));
     
    constraint(1)=norm(f,2)-sqrt(trans_maxpower);
    temp=(2*real(f_ini'*A_u*f)-f_ini'*A_u*f_ini)+noise_maxpower-r(1);
    constraint(2) = -temp;   
    B2=miu_u*temp_u+(1+miu_u)*miu_bs*diag(diag( temp_u ));
    [t1,t2]=eig(B2);
    constraint(3)=pow_pos(norm(f'*t1*sqrt(t2),2),2)+noise_maxpower -r(2);
    r_ini(2)=f_ini'*B2*f_ini+noise_maxpower;
    
    
    B3=temp_e+miu_bs*diag(diag( temp_e ));
    [t1,t2]=eig(B3);
    constraint(4)=pow_pos(norm(f'*t1*sqrt(t2),2),2)+noise_maxpower -r(3);
    r_ini(3)=f_ini'*B3*f_ini+noise_maxpower;
    temp=(2*real(f_ini'*A_e*f)-f_ini'*A_e*f_ini)+noise_maxpower-r(4);
    constraint(5) = -temp;
    
 
    constraint(6) = p(1)-log(r(1))/log(2);
    constraint(7) = log2(r_ini(2)) + (r(2)-r_ini(2)) / (r_ini(2)*log(2)) -p(2);
    constraint(8) = log2(r_ini(3)) + (r(3)-r_ini(3)) / (r_ini(3)*log(2)) -p(3);
    constraint(9) = p(4)-log(r(4))/log(2);
    
    
    maximize p(1) - p(2) - p(3) + p(4)
  
    subject to
    
         real(constraint)<=0;

cvx_end

if cvx_status(1)=='S' || cvx_status(3)=='a' 
    flag=1;
    rate_p=max(p(1) - p(2) - p(3) + p(4),0);  
    
    %%%%% TEST RATE
    rate_u_ori=log2( (trace(f_ini'*A_u*f_ini)+noise_maxpower)/(trace(f_ini'*B2*f_ini)+noise_maxpower));
    rate_e_ori=log2( (trace(f_ini'*B3*f_ini)+noise_maxpower)/(trace(f_ini'*A_e*f_ini)+noise_maxpower));
    rate_ori=rate_u_ori-rate_e_ori;
    
    rate_u_new=log2( (trace(f'*A_u*f)+noise_maxpower)/(trace(f'*B2*f)+noise_maxpower));
    rate_e_new=log2( (trace(f'*B3*f)+noise_maxpower)/(trace(f'*A_e*f)+noise_maxpower));
    rate_new=rate_u_new-rate_e_new;
    rate=max(rate_new,0); 
    %%%%%  END
else
    flag=0;
    f=ones(N,1);
    rate=0;
end

end

