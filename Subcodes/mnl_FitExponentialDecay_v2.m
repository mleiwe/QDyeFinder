function [beta,modelfunc]=mnl_FitExponentialDecay_v2(X,Y,beta0)
%This is a much more streamlined version of mnl_FitExponentialDecay. Where
%it doesn't work out the RMSE etc.. it simply provides the co-ordinates.
%Also it doesn't faff around making a table. The addition is that I prevent
%negative numbers
%
%Marcus Leiwe Spetember 2022
%% Make sure the it is all in a single column
if size(X,1)==1
    X=X';
elseif size(X,1)>1 && size(X,2)>1
    error("Your X values are not a single row/column")
end

if size(Y,1)==1
    Y=Y';
elseif size(Y,1)>1 && size(Y,2)>1
    error("Your Y values are not a single row/column")
end

%% Decay Model
modelfunc=  @(b,x)  b(1) * exp(-b(2)*x(:, 1));

%check if initial beta values exist
if isempty(beta0)==1 || exist('beta0','var')
    beta0= [1,0.01];
end
%% Fit
beta = nlinfit(X, Y, modelfunc, beta0);
end