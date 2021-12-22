
close all;
clear all;clc;
N          =10;            % array number of BS
M          =64;            % array number of IRS
K          = 1;            % number of users in each group
Rican_BR   =10;
Rican_RU   =10; Rican_RE   =10;
Rican_BU   =0;  Rican_BE   =0;
%%%%% Large scale path loss
PL_0=10^(-30/10); %dB the channel gain at the reference distance
x_bs=0;
y_bs=0;
x_irs=50;
y_irs=0;

d_BI=sqrt((x_irs-x_bs)^2+(y_irs-y_bs)^2); %m distance from the BS to IRS
pathloss_BR=sqrt(PL_0*(d_BI)^(-2.2));    % Large-scale pass loss from the BS to the IRS


x_user=50; y_user=2;
d_IU=sqrt((x_irs-x_user)^2+(y_irs-y_user)^2);  %m distance from the IRS to the users
pathloss_RU=sqrt(PL_0*(d_IU)^(-2.2));  % Large-scale pass loss from the IRS to the users
d_BU=sqrt((x_bs-x_user)^2+(y_bs-y_user)^2);%sqrt(d^2+d_v^2);  %m distance from the BS to the users
pathloss_BU=sqrt(PL_0*(d_BU)^(-3.6));  % Large-scale pass loss from the BS to the users

x_eve=45; y_eve=2;
d_IU=sqrt((x_irs-x_eve)^2+(y_irs-y_eve)^2);  %m distance from the IRS to the users
pathloss_RE=sqrt(PL_0*(d_IU)^(-2.2));  % Large-scale pass loss from the IRS to the users
d_BE=sqrt((x_bs-x_eve)^2+(y_bs-y_eve)^2);%sqrt(d^2+d_v^2);  %m distance from the BS to the users
pathloss_BE=sqrt(PL_0*(d_BE)^(-3.6));  % Large-scale pass loss from the BS to the users
%% Simulation loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_loop = 500; 

for loop = 1 : num_loop
    %% LOS of BS-RIS  %%%%%
    N_t = linspace(0,N-1,N).';
    N_M = linspace(0,M-1,M).';
    ang_AOD_BU=-pi/2+rand(1,1)*pi;
    theta_AOA_BU=0;
    steering_AOD_BR   = exp(-1j*pi*N_t*sin(ang_AOD_BU).');
    steering_AOA_BR   = exp(-1j*pi*N_M*sin(theta_AOA_BU));
    H_LOS=steering_AOA_BR*steering_AOD_BR';
    %%%%%  channel of BS-RIS   %%%%%
    H_NLOS=sqrt(1/2)*(randn(M,N) + sqrt(-1)*  randn(M,N));
    H=(sqrt(Rican_BR/(1+Rican_BR))*H_LOS+...
                         sqrt(1/(1+Rican_BR))*H_NLOS)*diag(pathloss_BR);
    %% LOS of user-IRS  %%%%%
    N_M = linspace(0,M-1,M).';
    sin_ang_IU=-pi/2+rand(K,1)*pi;
    pha_AOA_IU=pi/6;
    for m=1:M
        steering_AOA_IU(m,:)   = exp(-1j*pi* (floor(m/4)*sin_ang_IU.'*sin(pha_AOA_IU)...
                                 +(m-floor(m/4)*4)*sin_ang_IU.'*cos(pha_AOA_IU)) );
    end
    h_LOS=steering_AOA_IU;
    %%%%%  channel of RIS-user   %%%%%
    h_NLOS=sqrt(1/2)*(randn(M,1) + sqrt(-1)*  randn(M,1));
    h_IU=(sqrt(Rican_RU/(1+Rican_RU))*h_LOS+...
                         sqrt(1/(1+Rican_RU))*h_NLOS)*diag(pathloss_RU);
    
   
    %% LOS of user-IRS  %%%%%
    N_M = linspace(0,M-1,M).';
    sin_ang_IE=-pi/2+rand(K,1)*pi;
    pha_AOA_IE=pi/6;
    for m=1:M
        steering_AOA_IE(m,:)   = exp(-1j*pi* (floor(m/4)*sin_ang_IE.'*sin(pha_AOA_IE)...
                                 +(m-floor(m/4)*4)*sin_ang_IE.'*cos(pha_AOA_IE)) );
    end
    h_LOS=steering_AOA_IE;
    %%%%%  channel of RIS-user   %%%%%
    h_NLOS=sqrt(1/2)*(randn(M,1) + sqrt(-1)*  randn(M,1));
    h_IE=(sqrt(Rican_RE/(1+Rican_RE))*h_LOS+...
                         sqrt(1/(1+Rican_RE))*h_NLOS)*diag(pathloss_RE);
                     
                     
    %%     channel of BS-user   %%%%%  
    H_d_temp=sqrt(1/2)*(randn(N,1) + sqrt(-1)*  randn(N,1)); % small scale pass loss from the BS to the user
    h_dU=pathloss_BU*H_d_temp;
    %%     channel of BS-eve   %%%%%  
    H_d_temp=sqrt(1/2)*(randn(N,1) + sqrt(-1)*  randn(N,1)); % small scale pass loss from the BS to the user
    h_dE=pathloss_BE*H_d_temp;
              
    
    %%   Final channel
    G_U_all(:,:,loop)=[diag(h_IU)'*H; h_dU'];
    G_E_all(:,:,loop)=[diag(h_IE)'*H; h_dE'];
    
end

save('G_U_all','G_U_all');
save('G_E_all','G_E_all');