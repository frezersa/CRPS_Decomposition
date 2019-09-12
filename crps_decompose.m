function [meanCRPS_3,Reli,Reso,UNCR] = crps_decompose(fcst,obs)

% Decomposition of CRPS based on Hersbach 2000
% Written by Frezer Awol, Sept 2018
% awolf@mcmaster.ca

% Note that there are three ways to compute mean crps in the paper which might have different values - refer meanCRPS_3, meanCRPS_2 & meanCRPS_1

%% references
% Hersbach, H. (2000). Decomposition of the Continuous Ranked Probability Score for Ensemble Prediction Systems. 
% Weather and Forecasting, 15(5), 559–570. https://doi.org/10.1175/1520-0434(2000)015<0559:DOTCRP>2.0.CO;2


%% inputs
ens_mat = fcst; %ensembles: K*N, K is the time seires(or cases in Hersbach's term), N is the number of ensembles
obs = obs; % verifing analysis: K*1

%%
ens_sort = sort(ens_mat,2); 

%Commulative dist function forecasted by an ensemble system = p_x = pi
N=size(ens_sort,2);
K=size(ens_sort,1);
for k=1:K
    for i=1:N-1 % 0 < i < N
        P_x(k,i) = i/N; %P_x = pi eqn 23
    end
end
s=size(P_x);
%left outlier, i = 0
P_0 = zeros(s(1,1),1);
%right outlier, i = N
P_N = ones(s(1,1),1);
p_i=P_x;
P_x_all = [P_0 P_x P_N]; 

% decomposition
% first compute ci ~ CRPS
cntbn=0;cntou=0;
idx_0=[];idx_i=[];idx_N=[];
for k=1:K
%     for i=1:N-1 % 0 < i < N
    switch(obs(k) < ens_sort(k,1) || obs(k) > ens_sort(k,N))
        case 0
            cntbn=cntbn+1;
            for i=2:N % 0 < i < N, eqn 26    
                if obs(k) > ens_sort(k,i)
                    a(k,i) = ens_sort(k,i)-ens_sort(k,i-1);
                    b(k,i) = 0;
                    %compute ci, eqn 25
                    c_i_t(k,i) = a(k,i) * (p_i(k,i-1))^2;                    
                elseif ens_sort(k,i) > obs(k) && obs(k) > ens_sort(k,i-1)
                    a(k,i) = obs(k)-ens_sort(k,i-1);
                    b(k,i) = ens_sort(k,i)-obs(k);
                    %compute ci, eqn 25
                    c_i_t(k,i) = a(k,i) * (p_i(k,i-1))^2 + b(k,i) * (1 - p_i(k,i-1))^2;                   
                elseif obs(k) < ens_sort(k,i-1)
                    a(k,i) = 0;
                    b(k,i) = ens_sort(k,i)-ens_sort(k,i-1);
                    %compute ci, eqn 25
                    c_i_t(k,i) = b(k,i) * (1 - p_i(k,i-1))^2;                  
                end
            end
            idx_i=[idx_i; k];            
         
        case 1
            cntou=cntou+1;
            %outliers eqn 27
            if obs(k) < ens_sort(k,1)
                a(k,1) = 0; 
                b(k,i) = ens_sort(k,1) - obs(k);
                 %compute ci, eqn 25
                idx_0=[idx_0; k];                 
                c_i_t(k,1) = b(k,i) * (1 - P_0(k,1))^2;
                             
                
            elseif obs(k) > ens_sort(k,N)
                a(k,N+1) = obs(k) - ens_sort(k,N); 
                b(k,N+1) = 0; 
                 %compute ci, eqn 25  
                idx_N=[idx_N; k];
                c_i_t(k,N+1) = a(k,N+1) * (P_N(k,1))^2; 
                             
            end
    end
    
    %CRPS, eqn 24
    CRPS(k)=sum(c_i_t(k,:),2); 
end
%CRPS, eqn 24
meanCRPS_1 = mean(CRPS); 

a_avg=zeros(1,N+1);b_avg=zeros(1,N+1);
o_avg=zeros(1,N+1);g_avg=zeros(1,N+1);
cc_i=zeros(1,N+1);c_i_t=zeros(K,N+1);
Reli_i=zeros(1,N+1);C_pot_i=zeros(1,N+1);

% %Uncertainity (adapted from EVS):
inc = 1/N;
for i=1:K
    p_sam = 0;unc = 0;
    for j = 1:N-1
   start = round((i/K) * (length(obs_sort) - 1));
   stop = round((i+1)/K * (length(obs_sort) - 1));
   gap = obs_sort(stop) - obs_sort(start);
%    gap = ens_sort(i,j+1) - ens_sort(i,j);
   p_sam = p_sam + inc;
   p_sam_m(i,j) = p_sam;
   unc = unc + p_sam * (1 - p_sam) * gap;
   unc_m(i,j) = unc;
    end
end
%  UNCR = unc;
UNCR = mean(unc_m(:,N-1));

% average over cases
a_avg = mean(a); %give equal weight to each cases (1 up to k)
b_avg = mean(b);
switch(obs(k) < ens_sort(k,1) || obs(k) > ens_sort(k,N))
        case 0
            for i=2:N % 0 < i < N
                %eqn 28
                cc_i(i) = a_avg(i) * (p_i(1,i-1))^2 + b_avg(i) * (1 - p_i(k,i-1))^2; 

                %observed frequency 1, eqn 31
                o_avg(i) = b_avg(i)/(a_avg(i)+b_avg(i));
                %avg bin width, eqn 32
                g_avg(i) = a_avg(i) + b_avg(i); % mean(ens_sort(i+1) - ens_sort(i));
%             end
            %eqn 36-37
%             for i=2:N
                Reli_i(i) = g_avg(i) * (o_avg(i) - p_i(k,i-1))^2;
                C_pot_i(i) = g_avg(i) * o_avg(i) * (1 - o_avg(i));
            end
        case 1
            %outliers
            if obs(k) < ens_sort(k,1)
                %left outlier, %eqn 28
                cc_i(1) = b_avg(1) * (1 - P_0(1))^2; 
                %observed frequency for outliers, eqn 33
                for k=1:K
                    o_tmp_0(k) = ((1/K) * heaviside((ens_sort(k,1)-obs(k))));
                end
                o_avg(1) = sum(o_tmp_0,2);
                g_avg(1) = b_avg(1)/o_avg(1);   
                %relia and cpot
                Reli_i(1) = g_avg(1) * (o_avg(1) - P_0(1))^2;
                C_pot_i(1) = g_avg(1) * o_avg(1) * (1 - o_avg(1));
            elseif obs(k) > ens_sort(k,N)
                %right outlier, %eqn 28
                cc_i(N+1) = a_avg(N+1) * (P_N(1))^2;
                %observed frequency for outliers, eqn 33
                for k=1:K
                    o_tmp_N(k) = ((1/K) * heaviside((ens_sort(k,N)-obs(k))));
                end
                o_avg(N+1) = sum(o_tmp_N,2);
                g_avg(N+1) = a_avg(N+1)/(1-o_avg(N+1));  
                %relia and cpot
                Reli_i(N+1) = g_avg(N+1) * (o_avg(N+1) - P_N(1))^2;
                C_pot_i(N+1) = g_avg(N+1) * o_avg(N+1) * (1 - o_avg(N+1));                
            end
end
%eqn 28
% meanCRPS_2 = sum(cc_i,2); %%%%% mean CRPS alternative estimation 

%Check
% eqn34: input here

% get the decompositions, eqn 35
Reli = sum(Reli_i,2); % Reliability component
C_pot = sum(C_pot_i,2); % Potential CRPS component

meanCRPS_3 = Reli + C_pot; %%%%% mean CRPS from summation of the main decomposed components

Reso = UNCR - C_pot; % Resolution is dependent on Uncertainity component

% Draw reliability plot: fig 4,5,6 top
% Reliability plots can be drawn using P_x_all and S_o_avg

end
