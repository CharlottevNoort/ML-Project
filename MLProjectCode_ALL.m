%% Load data

load data_breastcancer

%% Average values for samples that belong to the same patient
% Sets of samples averaged here: 1 and 10; 2 and 68; 3 and 74

%In order to be able to combine two samples into a matrix (N), take them as columns
M = M';

%Replace column 1 by average of itself and column 10
N = [M(:,1) M(:,10)];
O = mean(N','omitnan');
M(:,1) = O';

%Replace column 2 by average of itself and column 68
N = [M(:,2) M(:,68)];
O = mean(N','omitnan');
M(:,2) = O';

%Replace column 3 by average of itself and column 74
N = [M(:,3) M(:,74)];
O = mean(N','omitnan');
M(:,3) = O';

clear N O

%Return matrix to 77x12553 orientation for further use
M = M';

%Remove remaining columns to keep one per patient
M(74,:) = [];
M(68,:) = [];
M(10,:) = [];

%% Assess amount of missing data

%To count number of NaNs per protein (column) in M:
MVP = sum(isnan(M));

%To calculate for each protein the percentage of patients with no data:
MVP_percent = MVP/77*100;
  
%To count number of proteins for each patient that have NaN:
MVS = sum(isnan(M'))';

%To calculate for each patient the percentage of proteins with no data:
MVS_percent = MVS/12553*100;

%To calculate total number of empty (NaN) cells:
sum_MV = sum(MVP)
percent_MV = sum_MV/966581*100

%% Remove features with large amount of missing values
% Using 10% and 20% as threshold for maximum percentage missing per protein
% Input: M from first section

nrows = size(M,1);

M10 = M(:,sum(isnan(M),1)/nrows*100<10)
clear nrows;

%% Impute missing values using k-Nearest Neighbors
% Using k=5 and k=8
% Default distance measure: Euclidean
% Input: M10 and M20 from previous section

M10_knn8 = knnimpute(M10',8)';

clear M10

%% Principal Component Analysis
% For reduction of dimensionality
% Input: Imputed matrices from previous section
%How can we make this a loop that saves output under different names for
%each input matrix?? To get it to do something like this, but in a loop:

display('Running Principal Component Analysis ...')

[loadings,score,latent,tsquared,explained,mu] = pca(M10_knn8, 'Economy', false);

%% Whitening

display('Whitening data ...')

X = score;
sigma = X * X' / size(X,2);
[Uwhite,Swhite,~] = svd(sigma);
clear sigma;
epsilon = 0.000001; %small constant
xPCAwhite = diag(1./sqrt(diag(Swhite(:,:)) + epsilon)) * Uwhite(:,:)' * X;

%% Plotting PCA vs Whitened PCA
% Shows that there are differences, for ecplanation of why whitening is needed

diagS = diag(Swhite);
ax3 = subplot(2,2,3);
plot(ax3,1:100,100*cumsum(PCAVar(1:100))/sum(PCAVar(1:100)),'r-^',...
    1:100,100*cumsum(diagS(1:100))/sum(diagS(1:100)),'b-o');
xlabel('Number of Principal Components');
ylabel('Percent Variance Explained in X');
legend({'PCA' 'PCA Whitening'},'location','SE');
grid('on');
clear diagS S PCAVar;

%% K-fold cross validation
% To determine optimal number of Principal Components

%for k = [5 8]
[RMSECV,CI] = pcaKFold(M10_knn8,10); %CI gives the confidence interval

errorbar(RMSECV,CI','DisplayName','10-fold Cross-validation',...
    'MarkerFaceColor',[1 0 0],...
    'MarkerEdgeColor',[1 0 0],...
    'Marker','*',...
    'Color',[0 0 0]);

grid('on');

ylim([0 4]);
xlabel('Number of Principal Components');
ylabel('Cross-validation error (RMSECV(c))');

%% Automatic evaluation of number of clusters

disp('Evaluation of number of clusters...')

for i=2:size(RMSECV,1)
    prev = RMSECV(i-1);
    value = RMSECV(i);
    
    %get the index, considering error below 0.03 and a growing tendency
    if ((prev < value) && (value < 0.03))
        maxCV = i-1;
        break;
    end
end

ResultCosine = [];
ResultEuclidean = [];

NR = 10;

%we have 4 subtypes of cancer on the database. 
for num_of_cluster = 1:4
    
    sCosine = 0;
    for i=1:NR %small number of observations
        idxCosinex = kmeans(xPCAwhite',num_of_cluster,'MaxIter',10000,...
            'distance','cosine','replicate',100);
        sCosine = sCosine + mean(silhouette(xPCAwhite',idxCosinex,'cosine'));
    end
    
    ResultCosine(:,num_of_cluster) = sCosine/NR;
    
    sEuclidean = 0;
    for i=1:NR %small number of observations
        idxEuclideanx = kmeans(xPCAwhite',num_of_cluster,'MaxIter',10000,...
            'distance','sqeuclidean','replicate',100);
        sEuclidean = sEuclidean + mean(silhouette(xPCAwhite',idxEuclideanx,'euclidean'));
     end
    
    ResultEuclidean(:,num_of_cluster) = sEuclidean/NR;
    
end

figure;

plot1 = plot([ResultCosine' ResultEuclidean'],'MarkerSize',6,'LineStyle','--');
set(plot1(1),'DisplayName','K-means (distance: cosine)',...
    'MarkerFaceColor',[1 0 0],...
    'MarkerEdgeColor',[1 0 0],...
    'Marker','o',...
    'Color',[0 0 0]);
set(plot1(2),'DisplayName','K-means (distance: euclidean)',...
    'MarkerFaceColor',[0 0 0],...
    'MarkerEdgeColor',[0 0 0],...
    'Marker','square',...
    'Color',[0 0 1]);

% Create xlabel
xlabel('Number of Clusters');
% Create ylabel
ylabel('Silhouette coefficient');
grid('on');
% Create legend
legend('show');

%% K-means Clustering

[IDX,C,sumd,d] = kmeans(M_PAM_20_KNN8,4,'distance','cosine','Replicates',20);

[B,ix] = sort(IDX);

imagesc(M_PAM_20_KNN8(ix,:))
lines = find(diff(B));
for l = lines'
    hl = refline(0,l);
    set(hl,'LineWidth',3,'color','k');
end
title('Expression of PAM50 proteins for patiens ordered by cluster')

D = squareform(pdist(M_PAM_20_KNN8(ix,:),'cosine'));
imagesc(D)
lines = find(diff(B));
for l = lines'
    hl = refline(0,l);
    set(hl,'LineWidth',3,'color','k');
end
title('Distance matrix ordered with kmeans 4 clusters')

%% Hierarchical clustering

X = score(:,1:54)

% Squared Euclidean distance
Y = pdist(X,'squaredeuclidean');
Z = linkage(Y,'weighted');

R = [];
for cluster_nr = 1:10
  C = cluster(Z,'maxclust',cluster_nr);
  [s,h] = silhouette(X,C);
  R = [ R; cluster_nr mean(s)];
end
figure;
plot( R(:,1),R(:,2),'r*-.');
xlabel('Number of clusters')
ylabel('Average Silhouette score')

figure;
[H, T] = dendrogram(Z,0,'Colorthreshold',44900);
xlabel('Patients')
ylabel('Distance (squared Euclidean)')
C = cluster(Z,'maxclust',4);

figure;
[s,h] = silhouette(X,C)
adjrand(C,PAM_groups)

% Cosine
Y = pdist(X,'squaredeuclidean');
Z = linkage(Y,'weighted');

R = [];
for cluster_nr = 1:10
  C = cluster(Z,'maxclust',cluster_nr);
  [s,h] = silhouette(X,C);
  R = [ R; cluster_nr mean(s)];
end
figure;
plot( R(:,1),R(:,2),'r*-.');
xlabel('Number of clusters')
ylabel('Average Silhouette score')
Y = pdist(X,'cosine');
Z = linkage(Y,'weighted');
figure;
[H, T] = dendrogram(Z,0,'Colorthreshold',0.98);
xlabel('Patients')
ylabel('Distance (cosine)')
C = cluster(Z,'maxclust',4);

figure;
[s,h] = silhouette(X,C)
adjrand(C,PAM_groups)

%% Fuzzy clustering 
%fcm(matrix,number_of_clusters)
       options=[NaN 3000 0 NaN]
       [center,U,obj_fcn] = fcm(data,3,options);
        maxU = max(U);
% U contains values for each patients association with each cluster
% maxU is the highest value for each patient
% gives each patient a label, 1-n, n=#of clusters        
         classify=zeros(77,1);
         for i=1:77;
             [max_value,row]=max(U(:,i));
             classify(i,1)=row;
         end
         
 %calc RAND, PAM_groups=patient clustered by cancer type        
        AR=adjrand(classify,PAM_groups)
                 

%% ANOVA
