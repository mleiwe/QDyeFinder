function [clusterIDs,centroids] = mnl_MeanShift(X,bandwidth)
%Perform Mean Shift Clustering on the input data X with the given bandwidth
%Params:
%X : Input matrix
%bandwidth: Bandwidth controlling kernel size
%Returns:
%clusterIDs, centroids

[n,p]=size(X);

%Initialize clusters to zero
clusterIDs=zeros(n,1);

%Init centroids
centroids = NaN(1,p);

for i=1:n
    x_i = X(i,:);
    shift = Inf;
    while norm(shift) > eps
        dist = pdist2(X,x_i);
        kernel = exp(-(dist.^2)/(2*bandwidth^2));
        shift = sum(bsxfun(@times, kernel, X), 1) / sum(kernel) - x_i;
        x_i = x_i+shift;
    end

    if isnan(centroids(1,:))
        %First cluster
        centroids(1,:) = x_i;
        clusterIDs(i) = 1;
    else
        dist = pdist2(x_i,centroids);
        [minDist,idx] = min(dist);
        if minDist < bandwidth
            clusterIDs(i) = idx;
            centroids(idx,:) = (centroids(idx,:)+x_i)/2;
        else
            centroids(end+1,:)=x_i;
            clusterIDs(i) = size(centroids,1);
        end
    end
end
end


