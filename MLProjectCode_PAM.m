%% Load PAM50 data and PAM50 grouping

load data_PAM50_77patients
load PAM50_groups

%% Check number of missing data per column

nrows = size(M_PAM,1);
plot(sum(isnan(M_PAM),1)/nrows*100)
title('Percentage of missing values')

%% Remove features with large amount of missing values
% Using 10% and 20% as threshold for maximum percentage missing per protein
% Input: M_PAM

M_PAM_10 = M_PAM(:,sum(isnan(M_PAM),1)/nrows*100<10); 

%% Impute missing values using k-Nearest Neighbors
% Using k=5 and k=8
% Default distance measure: Euclidean

M_PAM_10_KNN8 = knnimpute(M_PAM_10',8)';

%% Mean Centering
% Input: M_PAM_10_KNN8 from previous section

mc_M_PAM_10_KNN8 = ................................

%% Hierarchical Clustering

% Squared Euclidean
Y = pdist(mc_M_PAM_10_KNN8,'squaredeuclidean');
Z = linkage(Y,'weighted');
figure;
[H, T] = dendrogram(Z,0,'Colorthreshold',400);
xlabel('Patients')
ylabel('Distance (squared Euclidean)')
C = cluster(Z,'maxclust',4);
figure;
[s,h] = silhouette(mc_M_PAM_10_KNN8,C)
ARI_HsE = adjrand(C,PAM_groups)

% Cosine
Y = pdist(mc_M_PAM_10_KNN8,'cosine');
Z = linkage(Y,'weighted');
figure;
[H, T] = dendrogram(Z,0,'Colorthreshold',0.98);
xlabel('Patients')
ylabel('Distance (cosine)')
C = cluster(Z,'maxclust',4);
figure;
[s,h] = silhouette(mc_M_PAM_10_KNN8,C)
ARI_HC = adjrand(C,PAM_groups)

%% K-means Clustering

[IDX, C, SUMD, D] = kmeans(M_PAM_20_KNN8, 5, 'Distance', 'sqeuclidean', 'Replicates',20)
[IDX PAM_groups]
Table = crosstab(IDX,PAM_groups)
[AR,RI]=RandIndex(IDX,PAM_groups)

%% Fuzzy clustering 

       [center,U,obj_fcn] = fcm(data,3);
        maxU = max(U);
        % Find the data points with highest grade of membership in cluster 1
%          index1 = find(U(1,:) == maxU);
%          % Find the data points with highest grade of membership in cluster 2
%          index2 = find(U(2,:) == maxU);
%          index3 = find(U(3,:) == maxU);
%          index4 = find(U(4,:) == maxU);
         
         classify=zeros(77,1);
         for i=1:77;
             [max_value,row]=max(U(:,i));
             classify(i,1)=row;
         end
         
        [AR,RI]=RandIndex(classify,PAM_groups)


%% ANOVA


