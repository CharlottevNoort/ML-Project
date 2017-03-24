%Implementation of the Cross-validation by Eigenvector
%Source: Bro et al. (2008) Cross-validation of component models:
%a critical look at current methods

function [sumErrorPCA, CI] = pcaKFold(X, K)
indices = crossvalind('Kfold',size(X,1),K);

[U,W,V] = svd(X,'econ');
maxCV = size(diag(W),1);

X = U(:,1:maxCV)*W(1:maxCV,1:maxCV)*V(1:maxCV,1:maxCV)';

disp(strcat('Processing Cross-validation... '));

for k = 1:K %% K folds
    
    disp(strcat('Fold: ',num2str(k)));
    
    test = (indices == k); train = ~test;
    
    Xtrain = X(train,:);
    Xtest = X(test,:);
    
    avg = mean(Xtrain);
    Xtrain = bsxfun(@minus, Xtrain, avg);
    
    [~,~,VPCA] = svd(Xtrain);
    
    Xtest = bsxfun(@minus, Xtest, avg);
    
    
    %loop over the number of the left-out samples (get by the k-fold
    %method)
    for f=1:maxCV
        for j=1:size(Xtest,2)
            scorePCA = Xtest(:,[1:j-1 j+1:end])*pinv(VPCA([1:j-1 j+1:end],1:f))'*VPCA(:,1:f)';
            predErrPCA(j) = (Xtest(j) - scorePCA(j));
        end
        
        PRESSPCA(f) = sum(predErrPCA(:).^2);
    end
        
    %RMSECV error
    sumErrorPCA(:,k) = sqrt(PRESSPCA' / size(Xtest,1));
    
end

for i=1:size(sumErrorPCA,1)
    SEM = std(sumErrorPCA(i,:))/sqrt(size(sumErrorPCA,2));  % Standard Error
    ts = tinv([0.025  0.975],size(sumErrorPCA,2)-1);      % T-Score 95% confidence
    CI(i) = ts(1,2)*SEM;    % Confidence Interval
end

disp(strcat('Complete! '));

%Mean error each fold
sumErrorPCA = sum(sumErrorPCA,2) / K;
