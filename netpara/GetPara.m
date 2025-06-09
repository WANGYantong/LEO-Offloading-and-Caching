function para = GetPara(set)
% Getpara function generates simulation parameters of Table I.
% input set setting (the first row group in Tab.I), including:
%       set.N_satellite: # of satellites
%       set.N_service: # of service categories
%       set.N_index: # of NFs in each service chain
%       set.N_terminal: # of terminals
%       set.N_function: # of NFs 
%      
% output para (the second and third row group in Tab. I), including:
%       para.R_uj: 0-1 indicator for terminal u request of service j
%       para.A_us: 0-1 indicator for terminal u access via LEO s
%       para.V_ijk:0-1 indicator for ith function in service j is NF k
%       para.c_k: required cycles to execute 1 bit in NF k
%       para.f_ks:allocated computation resources of NF k in LEO s
%       para.C_s: total computation resources in LEO s
%       para.f_g: computation resources for each terminal in data center
%       para.l_k: required storage space of NF k
%       para.L_s: storage capacity in each LEO s
%       para.l_uij: input data size of ith function in service j from
%                   terminal u
%       para.r_{u,s,g}: data rate of uplink, ISL, and downlink
%       para.p_{u,s,g}: transmission power of uplink, ISL, and downlink
%       para.d_{u,s,g}: distance from terminal-to-satellite, 
%                       satellite-to-satellite, satellite-to-data center
%       para.kappa: energy consumption coefficient

I=set.N_index;
J=set.N_service;
K=set.N_function;
U=set.N_terminal;
S=set.N_satellite;

%% Indicators Generation (Second Row Group)
R_uj=zeros(U,J);
A_us=zeros(U,S);
V_ijk=zeros(I,J,K);

for ii=1:U
    R_uj(ii,randi([1,J]))=1;
    A_us(ii,randi([1,2]))=1; % only satellite 1 & 2 provide access service 
end

for jj=1:J
    order=sort(randperm(K,I));
    for ii=1:I
        V_ijk(ii,jj,order(ii))=1;
    end
end

para.R_uj=R_uj;
para.A_us=A_us;
para.V_ijk=V_ijk;

%% Other Parameters Generation (Third Row Group)
para.c_k=100;  % 100cycles/bit
para.f_ks=2e9; % default 2Gcycles/s  [0.5:0.5:3]GC/s
para.C_s=10e9; % default 10Gcycles/s [2:2:12]GC/s
para.f_g=2e9; % 2Gcycles/s
para.l_k=randi([8,40],K,1)*1e7; %80-400Mbit
para.L_s=8e8;  % default 800Mb, [400:200:1400]Mb

gamma=0.25; 
l_uij=ones(I,1)*1e8; % default 100Mb, [50:50:300]Mb
for ii=1:I-1
    l_uij(ii+1)=l_uij(ii)*gamma; 
end
para.l_uij=l_uij;

para.r_u=2e8; % 200Mbps
para.p_u=2; % 2Watt
para.r_s=1e10; % 10Gbps
para.p_s=1e3; % 30dbW
para.r_g=1e8; % 100Mbps
para.p_g=1e2; % 20dbW

para.d_u=1e6; % 1000km
para.d_g=3e6; % 2000km
para.d_s=8e5; % 800km

para.c=3e8; % light speed

para.kappa=1e-28; % 

end

