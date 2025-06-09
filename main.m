clear;
clc;

addpath(genpath(pwd));

%% Network Paramters
topo=GetTopo(); % topo.Hss-# of hops

set.N_satellite=length(topo.Hop);
set.N_terminal=10; %[5,10,15,20,25,30] default 10
set.N_service=10;
set.N_function=10;
set.N_index=4; % [1,2,3,4,5] default 4

Iter=100;
para=cell(Iter,1);
base_cost=cell(Iter,1);
sol_ILP=cell(Iter,1);
sol_Greedy=cell(Iter,1);
sol_NFC=cell(Iter,1);
parfor ii=1:Iter
    para{ii}=GetPara(set);
    para{ii}.alpha=0.5; % [0.1:0.1:0.9], default 0.5
    para{ii}.M=1e6; % sufficiently large number of linearization
    para{ii}.r_u=2e8*5/set.N_terminal; % 200Mbps /1,2,3,4,5,6 default 100

    % baseline
    base_cost{ii}=NoCachingSolver(topo,para{ii},set);
    %% Optimization Model
    % ILP solver
    sol_ILP{ii}=ILPSolver(topo,para{ii},set,base_cost{ii});

    % Greedy solver
    sol_Greedy{ii}=GreedySolver(topo,para{ii},set,base_cost{ii});

    % No-NFC solver
    sol_NFC{ii}=NFCSolver(topo,para{ii},set,base_cost{ii});
  
end

%% Performance Comparison
Iter=1e3;
cost_ILP=zeros(Iter,1);
cost_Greedy=zeros(Iter,1);
cost_NFC=zeros(Iter,1);
ILP_energy_ratio=zeros(Iter,1);
ILP_delay_ratio=zeros(Iter,1);
Greedy_energy_ratio=zeros(Iter,1);
Greedy_delay_ratio=zeros(Iter,1);
NFC_energy_ratio=zeros(Iter,1);
NFC_delay_ratio=zeros(Iter,1);
ILP_time=zeros(Iter,1);
Greedy_time=zeros(Iter,1);
NFC_time=zeros(Iter,1);
for ii=1:Iter
   cost_ILP(ii)=sol_ILP{ii}.fval2;
   ILP_energy_ratio(ii)=sol_ILP{ii}.energy_ratio;
   ILP_delay_ratio(ii)=sol_ILP{ii}.delay_ratio;
   Greedy_energy_ratio(ii)=sol_Greedy{ii}.energy_ratio;
   Greedy_delay_ratio(ii)=sol_Greedy{ii}.delay_ratio;
   NFC_energy_ratio(ii)=sol_NFC{ii}.energy_ratio;
   NFC_delay_ratio(ii)=sol_NFC{ii}.delay_ratio;
   cost_Greedy(ii)=sol_Greedy{ii}.fval2;
   cost_NFC(ii)=sol_NFC{ii}.fval2;
   ILP_time(ii)=sol_ILP{ii}.time;
   Greedy_time(ii)=sol_Greedy{ii}.time;
   NFC_time(ii)=sol_NFC{ii}.time;
end
