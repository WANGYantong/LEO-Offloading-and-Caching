function solution = GreedySolver(topo,para,set,base_cost)

I=set.N_index;
J=set.N_service;
K=set.N_function;
U=set.N_terminal;
S=set.N_satellite;

L_s=para.L_s*ones(S,1);
C_s=para.C_s*ones(S,1);

x=zeros(K,S);
y=zeros(U,K,S);
index_y=zeros(U,K,S);

tic;
for ii=1:I
    p_ks=transpose(transpose(para.A_us)*para.R_uj*...
        squeeze(para.V_ijk(ii,:,:)));
    for ss=1:S
        [value,index]=sort(p_ks(:,ss),'descend');
        for kk=1:K
            for uu=1:U
                if (value(kk)>0) && (para.l_k(index(kk))<=L_s(ss)) &&...
                        (para.R_uj(uu,:)*transpose(para.V_ijk(ii,:,index(kk)))>0)
                    if x(index(kk),ss)==0
                        x(index(kk),ss)=1;
                        L_s(ss)=L_s(ss)-para.l_k(index(kk));
                    end
                    if (para.f_ks<=C_s(ss)) &&...
                            index_y(uu,index(kk),ss)==0
                        y(uu,index(kk),ss)=1;
                        index_y(uu,index(kk),:)=1;
                        C_s(ss)=C_s(ss)-para.f_ks;
                    end
                elseif value(kk)>0 &&...
                        (para.R_uj(uu,:)*transpose(para.V_ijk(ii,:,index(kk)))>0)
                    [~,index2]=sort(topo.Hop(ss,:));
                    for ss2=2:S
                        for uu2=1:U
                            if para.l_k(index(kk))<=L_s(index2(ss2)) &&...
                                    (para.R_uj(uu2,:)*transpose(para.V_ijk(ii,:,index(kk)))>0)
                                if x(index(kk),index2(ss2))==0
                                    x(index(kk),index2(ss2))=1;
                                    L_s(index2(ss2))=L_s(index2(ss2))-para.l_k(index(kk));
                                end

                                if para.f_ks<=C_s(index2(ss2)) &&...
                                        index_y(uu2,index(kk),ss)==0
                                    y(uu2,index(kk),index2(ss2))=1;
                                    index_y(uu2,index(kk),:)=1;
                                    C_s(index2(ss2))=C_s(index2(ss2))-para.f_ks;
                                end
                            end
                        end
                    end
                elseif value(kk)==0
                    break;
                end
            end
        end
    end
end
time=toc;

solution=CostCalculator(x,y,topo,para,set);
solution.time=time;
solution.energy_ratio=(1-para.alpha)*solution.energy/base_cost.energy;
solution.delay_ratio=para.alpha*solution.delay/base_cost.delay;
solution.fval2=solution.energy_ratio+solution.delay_ratio;

end

