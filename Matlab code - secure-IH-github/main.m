%By Zhou Gui
%From 2019-10-19 to 
close all;
clear all;clc;
warning('off');
rand('twister',mod(floor(now*8640000),2^31-1));
%% Parameters Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% '1' stands for source-relays; '2' stands for relays-destination
N          = 4;            % array number of BS
M          = 64;            % array number of IRS

SNR_dB     = 10;     % dBW
%%%%% noise
N0=10^((-174-30) / 10); %-174dBm  
B=10^7; %10MHz
 noise_maxpower_original   = N0*B;            % % W
noise_maxpower_original   = 10^((-80-30) / 10);            % % W
%%%%% ends
temp=10:5:40;
trans_maxpower_all =10.^((temp-30) ./ 10); % trans_power=1 
miu_bs=0.01;
miu_u=0.01;
%% Simulation loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_loop = 500; 
load('G_U_all');
load('G_E_all');

for loop =1 : num_loop
    outerflag=1;   
%     T1=cputime;          
    G_U=[G_U_all(1:M,1:N,loop); G_U_all(size(G_U_all,1),1:N,loop)]/sqrt(noise_maxpower_original);
    G_E=[G_E_all(1:M,1:N,loop); G_E_all(size(G_U_all,1),1:N,loop)]/sqrt(noise_maxpower_original);
    noise_maxpower=noise_maxpower_original/noise_maxpower_original;

%     G_U=[G_U_all(1:M,1:N,loop); G_U_all(size(G_U_all,1),1:N,loop)];
%     G_E=[G_E_all(1:M,1:N,loop); G_E_all(size(G_U_all,1),1:N,loop)];
%     noise_maxpower=noise_maxpower_original;
%%  For different SNR  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('  loop |  num_J  |  SNR  |  i  |  trans_SNR | relay_SNR \n');
for i_p     = 1 : length(trans_maxpower_all)
 
    t0=cputime;
    trans_maxpower=trans_maxpower_all(i_p);
     
    %%%%%  Initialize F and e beamforming  %%%%%
    e_temp=randn(M+1,1) + sqrt(-1)*  randn(M+1,1);
    e_ini=exp(1j*angle(e_temp));
    e_ini=ones(M+1,1);

    f_ini=ones(N,1)*sqrt(trans_maxpower/(N));
%     f_ini=randn(N,1)*sqrt(trans_maxpower/(N));
    f(:,1)=f_ini;
    e(:,1)=e_ini;

    num_iterative = 15000;
    for n  = 1 : num_iterative
                 %%%%%  Optimize F  %%%%%
        [f_1,rate_sub,rate_p,flag_f] = Generate_beamforming_F(N, M, G_U, G_E, ...
                    f(:,n), e(:,n),  noise_maxpower, trans_maxpower,miu_bs,miu_u );
        rate_f(n+1)=real(rate_sub);
        rate_obj(n+1)=real(rate_p);
        f(:,n+1)=f_1;
        if flag_f==0
            outerflag=0;
            break;
        end    

                %%%%%  Optimize e  %%%%%
        [e_1,rate_sub,flag_e] = Generate_beamforming_e(N, M, G_U, G_E, ...
                    f(:,n+1), e(:,n),  noise_maxpower, trans_maxpower,miu_bs,miu_u );
        rate_e(n+1)=real(rate_sub);
        e(:,n+1)=e_1;
        if flag_e==0
            outerflag=0;
            break;
        end     
      
        %%%%%  stop criterion  %%%%%
        
        fprintf('   %g  |  %g  |  %g  \n',loop, trans_maxpower_all(i_p), n);
        if abs(rate_f(n+1)-rate_f(n))<10^(-4) 
            break;
        end
        
    end
    if outerflag==0
        break;
    end
 
    Rate(loop,i_p)=real(rate_f(n+1));

    t2=cputime;
    CPU_Time(loop,i_p)=t2-t0;
    iteration(loop,i_p)=n;

end
    save('Rate_RIS_001','Rate');
    T2=cputime;
  
end
a=1;
    
