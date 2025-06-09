function base_cost = NoCachingSolver(topo,para,set)

K=set.N_function;
U=set.N_terminal;
S=set.N_satellite;

x=zeros(K,S);
y=zeros(U,K,S);

tic;
base_cost=CostCalculator(x,y,topo,para,set);
base_cost.time=toc;

end

