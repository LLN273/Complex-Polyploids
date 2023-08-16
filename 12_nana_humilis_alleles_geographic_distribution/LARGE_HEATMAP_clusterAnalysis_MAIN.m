
% k-mean cluster analysis of dataset that shows frequency with which an introgressed allele is observed, by population

% Luis Leal
% 2022

clear
close all hidden
format shortG


%______________________________________________________________________________________________________________________________  

% Load data

    MPoly_GE = readtable('db_nana_pop_ALL.txt');               % nana data; all genes  
    %MPoly_GE = readtable('db_hum_pop_ALL.txt');               % humilis data; all genes
    
    
%___________________________________________________________________________________________________________________________  

%%%%%% k-mean clustering analysis: conditions
N_runs = 50;                                                            % number of runs for each k-value   
Kmax = 10;                                                              % maximum k-value considered

% clustering analysis
GE_matrix = [MPoly_GE{:,2:8}];
[  totalSumsofDistsToClusterCentroids_1, delta_totalSumsofDistsToClusterCentroids  ] = clustering_control_corr(GE_matrix, N_runs, Kmax, 1, 'k-means clustering analysis' );    

%%%%%% Data for best k (k*) 
[clustering_new,ClusterCentroids_new,sumsofDistsToClusterCentroids_new] = kmeans(GE_matrix,6);    % clustering for best k

% number of elememts per cluster
tag_count(1)= sum(clustering_new(:) == 1)
tag_count(2)= sum(clustering_new(:) == 2)
tag_count(3)= sum(clustering_new(:) == 3)
tag_count(4)= sum(clustering_new(:) == 4)
tag_count(5)= sum(clustering_new(:) == 5)
tag_count(6)= sum(clustering_new(:) == 6)
tag_count(7)= sum(clustering_new(:) == 7)
%tag_count(8)= sum(clustering_new(:) == 8)
%tag_count(9)= sum(clustering_new(:) == 9)

% reorganize expression matrix according to clustering label
matrix_aux = [MPoly_GE table(clustering_new)];
matrix_by_label = sortrows(matrix_aux,9);

% save matrix to file
writetable(matrix_by_label,'presence_matrix_after_kmeans_nana.txt')
%writetable(matrix_by_label,'presence_matrix_after_kmeans_hum.txt')





