   function Value=Devation(X, labels) 
%==========================================================================
% FUNCTION: cp = valid_compactness(X, labels) 
% DESCRIPTION: A function for computing Compactness validity index 
%              for a clustering result
%
% INPUTS:  X = a dataset, rows of X correspond to observations; columns
%              correspond to variables (exclude class labels!!)
%     labels = cluster labels from a clustering result (N-by-1 vector)
%
% OUTPUT: cp = compactness score
%==========================================================================
% copyright (c) 2010 Iam-on & Garrett
%==========================================================================

n = length(labels); % number of data points
C = unique(labels); %all clusters
k = length(C); %number of clusters
cp = 0; %initialize compactness
for i=1:k %for each cluster
    ind = find(labels==C(i)); %find data point members for the i-th cluster
    center=sum(X(ind,:))./size(ind,1);
    sum_d = 0;

    for ij=1:1:size(ind,1)
        ix = X(ind(ij),:);
        xx=norm(center-ix);
     
        sum_d= sum_d+ xx;
    end

    cp = cp + sum_d;
end
Value=cp;