function [mcx,mx] = mncn(x)
%MNCN Mean center scales matrix to mean zero.
%  Mean centers matrix (x), returning a matrix with
%  mean zero columns (mcx) and the vector of means
%  (mx) used in the scaling.
%
%From Biosystems Data Analysis course (UvA)
%
%I/O: [mcx,mx] = mncn(x);


[m,n] = size(x);
mx    = mean(x);
mcx   = (x-mx(ones(m,1),:));
