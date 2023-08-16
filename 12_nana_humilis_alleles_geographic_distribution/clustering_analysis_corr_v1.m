function [ totalSumsofDistsToClusterCentroids ] = clustering_analysis_v2( my_data, N_runs, k_max )

   % my_data: data array (columns correspond to variables)
   % N_runs: number of runs for each k value
   % k_max: maximum k value considered

totalSumsofDistsToClusterCentroids = zeros(k_max,1);

for k=1:k_max
    
    for j = 1:N_runs
        [idx,ClusterCentroids,sumsofDistsToClusterCentroids] = kmeans(my_data,k,'dist','corr');

    % Note: When using kmeans(X,k), rows of X correspond to points and columns correspond to variables.
        
        totalSumsofDistsToClusterCentroids(k,j)=sum(sumsofDistsToClusterCentroids);     % total sum of Distances to the Cluster Centroids
    
    end
    
end

%totalSumsofDistsToClusterCentroids(1,:)=[];


end

