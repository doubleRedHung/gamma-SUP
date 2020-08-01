% gamma-SUP (Hung, 09/30/2013)
%
% ----- INPUT -----
% X: n by p data matrix 
% s: minor tuning parameter (suggest to be 0.025)
% tau: critical tuning parameter, needs to be tuned by phase-transition
% blurring: 1 = blurrning (default)
%
% ----- OUTPUT -----
% c_id: cluster membership id
% Y_new: updated data points at convergence

function [c_id, Y_new] = gamma_sup(X, s, tau, blurring)

if nargin < 4
    blurring = 1;
end
err_th = max(max(abs(X)))/10000; % error threshold for stopping
err = 999;
max_iter = 150;
[n,p] = size(X); % n: sample size; p: dim
Y_new = X/tau;

Y_0 = X/tau;
Y3_0 = diag(Y_0*Y_0')*ones(1,n);

iter_count=0;
while err>err_th && iter_count < max_iter
    
    Y=Y_new; % n x p; updated
    Y2=Y*Y'; % n x n
    Y3=diag(Y2)*ones(1,n);
    if blurring ==1
        YD = Y3 + Y3'- 2*Y2; % n x n; (i,j)th entry: || y_i-y_j ||^2
    else  
        YD = Y3 + Y3_0'- 2*(Y*Y_0');   % for non-blurring 
    end
    YD=max(YD,0); % dis-similarity matrix
    W=exp_q(-YD, 1-s);
    W=W./(sum(W,2)*ones(1,n));  
    Y_new=W*Y;  % update \mu_i; weight on model
    
    err=max(max(abs(Y-Y_new)));
    iter_count = iter_count+1;
end

% assign cluster membership
c_id=zeros(n,1);
C=ceil(log10(err_th));
Z=round(Y_new/10^C) * 10^C;
c_n=1;
while sum(c_id==0)>0
    k=find(c_id==0);
    k=k(1);
    Zd=Z-ones(n,1)*Z(k,:);
    Zd=sum(abs(Zd),2);
    tmp= Zd<=p*10^C;
    c_id(tmp)=c_n; % assign data indexed by tmp to cluster c_n
    c_n=c_n+1;
end


% q-exponential: e_q(u)
function f = exp_q(u, q)
    f = max(0, 1+(1-q)*u);
f = f.^(1/(1-q));







