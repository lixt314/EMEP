function [fitness,label]=objectivefunction(PredictClass,para,K,truelabels,X1)
dcASRS = 0.8;
R = 5;
dcSRS = 0.8;
dcCTS = 0.8;
%%%%1 for cts, 2 forASRS, 3 for SRS

if para(1)==1   
S = cts(PredictClass, dcCTS);
elseif para(1)==2
S = asrs(PredictClass, dcASRS);    
else
S = srs(PredictClass, dcSRS, R);   
end

%%%%1 for kmeans, 2 forAL, 3 for CL, 4 for 
dis = stod(S); %convert similarity matrix to distance vector


if para(2)==1   
label = kmeans(S,K);
while length(unique(label)) ~= K;
        label = kmeans(S,K);
end
elseif para(2)==2
fig_flag=0;
label1 = spectral_clustering(S, K, fig_flag);
[~,label]=max(label1');
label=label';
elseif para(2)==3
dist = pdist2(S,S);
[cluster_lables] = cluster_dp(dist,K);
label=cluster_lables';
end

V = cleval(X1, label, truelabels);

fitness=V(1:2);
%fitness(1)= V(2);
%fitness(2)= V(3);

