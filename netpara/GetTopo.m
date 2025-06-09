function topo = GetTopo()
% GetTopo Function generates the LEO network topology
% Input: Null
% Output: topo.Hop: # of hops between two satellites

debug=0; % 1 enable drawing graph

%% graph generation
s=[1,1,1,1,2,2,2,5,3];
t=[2,3,4,5,6,7,8,6,8];
X_axis=[2,3,2,1,2,3,4,3];
Y_axis=[2,2,1,2,3,3,2,1];
weights=ones(size(s));

G=graph(s,t,weights);

if debug==1
    plot(G,'EdgeLabel',G.Edges.Weight,'LineWidth',1.5,'Marker','o',...
        'MarkerSize',10,'XData',X_axis,'YData',Y_axis);
end

%% topology parameter
source=1:8;
destination=1:8;

Hop=zeros(length(source),length(destination));

for ii=1:length(source)
    for jj=1:length(destination)
        [~,Hop(ii,jj)]=shortestpath(G,source(ii),destination(jj));
    end
end

topo.Hop=Hop;

end

