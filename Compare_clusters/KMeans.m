%Non-weighted K-means clustering
function [cluster_idx,centroids] = KMeansClustering(X,K)
%Inputs:
%X: data matrix
%K: number of clusters
%Outputs:
%cluster_idx: cluster index for each data point
%centroids: cluster centroids

%Initialization
[N,D] = size(X);

%Random init centroids
rand_idx = randperm(N,K);
centroids = X(rand_idx,:);

%Setup params
start = val;
max_iter = val;
threshold = val;

while start < max_iter
    %Assign data to nearest centroids
    cluster_idx = zeros(N,1);
    for i = 1:N
        dist = zeros(K,1);
        for j = 1:K
            dist(j) = norm(X(i,:)-centroids(j,:));
        end
        [~,cluster_idx(i)] = min(dist);
    end

    %Update centroids
    centroids1 = zeros(K,D);
    c = zeros(K,1);
    for i = 1:N
        centroids1(cluster_idx(i),:) = centroids1(cluster_idx(i),:) + X(i,:);
        c(cluster_idx(i)) = c(cluster_idx(i)) + 1;
    end
    for i=1:K
        centroids1(i,:) = centroids1(i,:)/c(i);
    end

    %Convergence check
    start = start + 1;
    if sum(sqrt(sum((centroids1-centroids).^2,2))) < threshold
        break;
    end
    centroids = centroids1;
end
