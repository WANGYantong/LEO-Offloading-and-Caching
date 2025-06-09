function cost = CostCalculator(x,y,topo,para,set)

I=set.N_index;
J=set.N_service;
K=set.N_function;
U=set.N_terminal;
S=set.N_satellite;

%% Constraints Evaluation
constr1=para.f_ks*sum(sum(y))<=para.C_s;

if all(constr1,'all')
    cost.constr1=0;
else
    cost.constr1=sum(1-constr1,'all');
end


% l_k-->l_ks
l_ks=repmat(para.l_k,1,S);
constr2=squeeze(sum(l_ks.*x,1))<=para.L_s;

if all(constr2,'all')
    cost.constr2=0;
else
    cost.constr2=sum(1-constr2,'all');
end

% x_ks-->x_uks
x_ex=reshape(x,1,K*S);
x_ex=repmat(x_ex,[U,1]);
x_ex=reshape(x_ex,U,K,S);

constr3=y<=x_ex;
if all(constr3,'all')
    cost.constr3=0;
else
    cost.constr3=sum(1-constr3,'all');
end

constr4=sum(y,3)<=1;
if all(constr4,'all')
    cost.constr4=0;
else
    cost.constr4=sum(1-constr4,'all');
end

%% Auxiliary Building
% z_uk1s1k2s2
y_ex1=reshape(y,1,U*K*S);
y_ex1=repmat(y_ex1,[1,K,S]);
y_ex1=reshape(y_ex1,U,K,S,K,S);

y_ex2=repmat(y,[S,1,1,1]);
y_ex2=reshape(y_ex2,U,S,K,S);
y_ex2=repmat(y_ex2,[K,1,1,1,1]);
y_ex2=reshape(y_ex2,U,K,S,K,S);

z=y_ex1.*y_ex2;

% theta_iu
buffer=zeros(I-1,U,K);
buffer2=zeros(I-1,U,K);
for ii=1:I-1
    for uu=1:U
        for kk=1:K
            buffer(ii,uu,kk)=para.R_uj(uu,:)*transpose(para.V_ijk(ii,:,kk));
            buffer2(ii,uu,kk)=para.R_uj(uu,:)*transpose(para.V_ijk(ii+1,:,kk));
        end
    end
end
buffer_ex=reshape(buffer,1,(I-1)*U*K);
buffer_ex=repmat(buffer_ex,[1,S]);
buffer_ex=reshape(buffer_ex,(I-1),U,K,S);
buffer2_ex=reshape(buffer2,1,(I-1)*U*K);
buffer2_ex=repmat(buffer2_ex,[1,S]);
buffer2_ex=reshape(buffer2_ex,(I-1),U,K,S);
%y_uks-->y_iuks
y_ex=reshape(y,1,U*K*S);
y_ex=repmat(y_ex,[(I-1),1]);
y_ex=reshape(y_ex,(I-1),U,K,S);
g1_auxi=sum(buffer_ex.*y_ex,[3,4]);
g2_auxi=sum(buffer2_ex.*y_ex,[3,4]);
theta=xor(1-g1_auxi,1-g2_auxi);

% pi_u
buffer=zeros(U,K);
for uu=1:U
    for kk=1:K
        buffer(uu,kk)=sum(para.R_uj(uu,:)*transpose(para.V_ijk(:,:,kk)));
    end
end
pi=zeros(U,1);
for uu=1:U
    if sum(buffer.*(1-sum(y,3)),2)>=1
        pi(uu)=1;
    end
end

%% Cost Computing

%%%% latency model %%%%
%%%% U-->S
%%% transmission delay
delay_buffer=sum(para.R_uj*para.l_uij(1),2);
delay_buffer=repmat(delay_buffer,1,S);
delay_trans_US=sum(para.A_us.*delay_buffer,'all')/para.r_u;

delay_US=delay_trans_US+2*(para.d_u/para.c)*U;

%%%% S-->S'
%%% transmission delay
%R_uj*l_1*V_1jk--->uks
delay_buffer=para.R_uj*para.l_uij(1)*squeeze(para.V_ijk(1,:,:));
delay_buffer=reshape(delay_buffer,1,U*K);
delay_buffer=repmat(delay_buffer,[1,S]);
delay_buffer=reshape(delay_buffer,U,K,S);
%A_us--->A_uks
delay_buffer2=repmat(para.A_us,[K,1,1]);
delay_buffer2=reshape(delay_buffer2,U,K,S);
%u,k,s1,s2
delay_buffer3=delay_buffer.*delay_buffer2;
delay_buffer3=reshape(delay_buffer3,1,U*K*S);
delay_buffer3=repmat(delay_buffer3,[1,S]);
delay_buffer3=reshape(delay_buffer3,U,K,S,S);
%y_uks2--->y_uks1s2
y_ex=repmat(y,[1,S]);
y_ex=reshape(y_ex,U,K,S,S);

delay_buffer=squeeze(sum(sum(delay_buffer3.*y_ex)));
delay_trans_SS1=sum(delay_buffer.*topo.Hop,"all")/para.r_s;

%l_(i+1)*V_ijk1*V_(i+1)jk2--->ujk1k2
delay_buffer=zeros(J,K,K);
for jj=1:J
    for kk1=1:K
        for kk2=1:K
            delay_buffer(jj,kk1,kk2)=sum(para.l_uij(2:I).*para.V_ijk(1:I-1,jj,kk1).*...
                para.V_ijk(2:I,jj,kk2));
        end
    end
end
delay_buffer=reshape(delay_buffer,1,J*K*K);
delay_buffer=repmat(delay_buffer,[U,1]);
delay_buffer=reshape(delay_buffer,U,J,K,K);
%R_uj--->R_ujk1k2
delay_buffer2=reshape(para.R_uj,1,U*J);
delay_buffer2=repmat(delay_buffer2,[1,K,K]);
delay_buffer2=reshape(delay_buffer2,U,J,K,K);
%ujk1k2--->uk1k2--->uk1s1k2s2
delay_buffer3=squeeze(sum(delay_buffer.*delay_buffer2,2));
delay_buffer3=repmat(delay_buffer3,[1,S]);
delay_buffer3=reshape(delay_buffer3,U,K,S,K);
delay_buffer3=reshape(delay_buffer3,1,U*K*S*K);
delay_buffer3=repmat(delay_buffer3,[1,S]);
delay_buffer3=reshape(delay_buffer3,U,K,S,K,S);

delay_buffer=squeeze(sum(delay_buffer3.*z,[1,2,4]));
delay_trans_SS2=sum(delay_buffer.*topo.Hop,"all")/para.r_s;

delay_trans_SS=delay_trans_SS1+delay_trans_SS2;

%%% propagation delay
%R_uj*V_1jk--->uks
delay_buffer=para.R_uj*squeeze(para.V_ijk(1,:,:));
delay_buffer=reshape(delay_buffer,1,U*K);
delay_buffer=repmat(delay_buffer,[1,S]);
delay_buffer=reshape(delay_buffer,U,K,S);
%A_us--->A_uks
delay_buffer2=repmat(para.A_us,[K,1,1]);
delay_buffer2=reshape(delay_buffer2,U,K,S);
%u,k,s1,s2
delay_buffer3=delay_buffer.*delay_buffer2;
delay_buffer3=reshape(delay_buffer3,1,U*K*S);
delay_buffer3=repmat(delay_buffer3,[1,S]);
delay_buffer3=reshape(delay_buffer3,U,K,S,S);
%y_uks2--->y_uks1s2
y_ex=repmat(y,[1,S]);
y_ex=reshape(y_ex,U,K,S,S);

delay_buffer=squeeze(sum(sum(delay_buffer3.*y_ex)));
delay_propa_SS1=sum(delay_buffer.*topo.Hop,"all")*para.d_s/para.r_s;

%R_uj*V_1jk--->uks
delay_buffer=para.R_uj*squeeze(para.V_ijk(I,:,:));
delay_buffer=reshape(delay_buffer,1,U*K);
delay_buffer=repmat(delay_buffer,[1,S]);
delay_buffer=reshape(delay_buffer,U,K,S);
%A_us--->A_uks
delay_buffer2=repmat(para.A_us,[K,1,1]);
delay_buffer2=reshape(delay_buffer2,U,K,S);
%u,k,s1,s2
delay_buffer3=delay_buffer.*delay_buffer2;
delay_buffer3=reshape(delay_buffer3,1,U*K*S);
delay_buffer3=repmat(delay_buffer3,[1,S]);
delay_buffer3=reshape(delay_buffer3,U,K,S,S);
%y_uks2--->y_uks1s2
y_ex=repmat(y,[1,S]);
y_ex=reshape(y_ex,U,K,S,S);

delay_buffer=squeeze(sum(sum(delay_buffer3.*y_ex)));
delay_propa_SS2=sum(delay_buffer.*topo.Hop,"all")*para.d_s/para.r_s;

%V_ijk1*V_(i+1)jk2--->ujk1k2
delay_buffer=zeros(J,K,K);
for jj=1:J
    for kk1=1:K
        for kk2=1:K
            delay_buffer(jj,kk1,kk2)=sum(para.V_ijk(1:I-1,jj,kk1).*...
                para.V_ijk(2:I,jj,kk2));
        end
    end
end
delay_buffer=reshape(delay_buffer,1,J*K*K);
delay_buffer=repmat(delay_buffer,[U,1]);
delay_buffer=reshape(delay_buffer,U,J,K,K);
%R_uj--->R_ujk1k2
delay_buffer2=reshape(para.R_uj,1,U*J);
delay_buffer2=repmat(delay_buffer2,[1,K,K]);
delay_buffer2=reshape(delay_buffer2,U,J,K,K);
%ujk1k2--->uk1k2--->uk1s1k2s2
delay_buffer3=squeeze(sum(delay_buffer.*delay_buffer2,2));
delay_buffer3=repmat(delay_buffer3,[1,S]);
delay_buffer3=reshape(delay_buffer3,U,K,S,K);
delay_buffer3=reshape(delay_buffer3,1,U*K*S*K);
delay_buffer3=repmat(delay_buffer3,[1,S]);
delay_buffer3=reshape(delay_buffer3,U,K,S,K,S);

delay_buffer=squeeze(sum(delay_buffer3.*z,[1,2,4]));
delay_propa_SS3=sum(delay_buffer.*topo.Hop,"all")*para.d_s/para.r_s;

delay_propa_SS=delay_propa_SS1+delay_propa_SS2+delay_propa_SS3;


%%% computation delay
delay_buffer=zeros(U,K);
for uu=1:U
    for kk=1:K
        delay_buffer(uu,kk)=transpose(para.l_uij)*para.V_ijk(:,:,kk)...
            *transpose(para.R_uj(uu,:));
    end
end
delay_buffer=reshape(delay_buffer,1,U*K);
delay_buffer=repmat(delay_buffer,[1,S]);
delay_buffer=reshape(delay_buffer,U,K,S);

delay_comput_SS=sum(delay_buffer.*y,'all')*para.c_k/para.f_ks;

delay_SS=delay_trans_SS+delay_propa_SS+delay_comput_SS;

%%%% S-->G
%%% transmission delay
delay_buffer=zeros(U,K);
for uu=1:U
    for kk=1:K
        delay_buffer(uu,kk)=para.l_uij(1)*para.R_uj(uu,:)*...
            transpose(para.V_ijk(1,:,kk));
    end
end
delay_buffer2=1-sum(y,3);
delay_trans_SG1=sum(delay_buffer.*delay_buffer2,'all')/para.r_g;

delay_buffer3=zeros(I-1,U);
for ii=1:I-1
    for uu=1:U
        delay_buffer3(ii,uu)=sum(para.l_uij(ii+1)*para.R_uj(uu,:));
    end
end
delay_trans_SG2=sum(delay_buffer3.*theta,'all')/para.r_g;
delay_trans_SG=delay_trans_SG1+delay_trans_SG2;

%%% computation delay
delay_buffer=zeros(I,U,K);
for ii=1:I
    for uu=1:U
        for kk=1:K
            delay_buffer(ii,uu,kk)=para.l_uij(ii)*para.R_uj(uu,:)*...
                transpose(para.V_ijk(ii,:,kk));
        end
    end
end
delay_buffer=squeeze(sum(delay_buffer));
delay_buffer2=1-sum(y,3);
delay_comput_SG=sum(delay_buffer.*delay_buffer2*para.c_k,'all')/para.f_g;

%%% propagation delay
% delay_buffer=zeros(U,K);
% for uu=1:U
%     for kk=1:K
%         delay_buffer(uu,kk)=sum(para.R_uj(uu,:)*para.V_ijk(:,:,kk));
%     end
% end
delay_propa_SG=sum(pi)*2*para.d_g/para.c;

delay_SG=delay_trans_SG+delay_comput_SG+delay_propa_SG;

delay=delay_US+delay_SS+delay_SG;

%%%% energy model %%%%
%%%% U-->S
energy_trans_US=para.p_u*delay_trans_US;

%%%% S-->S'
energy_trans_SS=para.p_s*delay_trans_SS;

energy_buffer=zeros(U,K);
for uu=1:U
    for kk=1:K
        energy_buffer(uu,kk)=para.R_uj(uu,:)*transpose(para.V_ijk(:,:,kk))...
            *para.l_uij;
    end
end
energy_buffer2=sum(y,3);
energy_comput_SS=sum(energy_buffer.*energy_buffer2,'all')*para.c_k*...
    (para.f_ks)^2*para.kappa;

%%%% S-->G
energy_trans_SG=para.p_g*delay_trans_SG;

energy=energy_trans_US+energy_trans_SS+energy_trans_SG+energy_comput_SS;


cost.energy=energy;
cost.energy_trans_US=energy_trans_US;
cost.energy_trans_SS=energy_trans_SS;
cost.energy_trans_SG=energy_trans_SG;
cost.energy_comput_SS=energy_comput_SS;

cost.delay=delay;
cost.delay_US=delay_US;
cost.delay_SS=[delay_trans_SS,delay_propa_SS,delay_comput_SS];
cost.delay_SG=[delay_trans_SG,delay_comput_SG,delay_propa_SG];

cost.fval=para.alpha*delay+(1-para.alpha)*energy;

end

