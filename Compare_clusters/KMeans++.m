function [cluster_labels, centroids] = KMeansPlusPlus(data, K, true_labels)
% Inputs:
%   data: NxD matrix of N data points in D dimensions
%   K: number of clusters
%   true_labels: Nx1 vector of true labels for the data points (optional)
% Outputs:
%   cluster_labels: Nx1 vector indicating the cluster index for each data point
%   centroids: KxD matrix of K centroids in D dimensions

% Initialize variables
[N, D] = size(data);
centroids = zeros(K, D);

% Choose first centroid randomly
centroids(1,:) = data(randi(N),:);

% Repeat for each remaining centroid
for i = 2:K
    % Compute the squared distance of each data point to the nearest centroid
    d = zeros(N, 1);
    for j = 1:N
        d(j) = min(vecnorm(data(j,:) - centroids, 2, 2).^2);
    end
    
    % Normalize the distances
    d = d / sum(d);
    
    % Choose the next centroid using a weighted random sample
    r = rand();
    cumsum_d = cumsum(d);
    centroids(i,:) = data(find(cumsum_d >= r, 1), :);
end

% Run K-means on the initial centroids
[cluster_labels, centroids] = KMeans(data, K, centroids);
