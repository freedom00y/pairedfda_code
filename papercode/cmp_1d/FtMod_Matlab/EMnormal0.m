function [theta,sigma2,AIC,BIC] = EMnormal0(id,t,x,a,b,nk,order,lsm)
%function [theta,sigma2,AIC,BIC] = EMnormal0(id,t,x,a,b,nk,order,lsm)
%Semiparametric Normal-model estimators of the mean (via EM algorithm)
%
%NOTE: since the time grids may be sparse and irregular, the data is input 
% as a concatenated vector rather than a matrix; a vector of labels is used
% to identify the individuals
%
%INPUT:
%   id:     Individual labels           (N x 1, char, cell or numeric)
%   t:      Time grid                   (N x 1)
%   x:      Observations                (N x 1)
%   a:      Time range, lower bound     (scalar)
%   b:      Time range, upper bound     (scalar)
%   nk:     Number of knots (excluding endpoints)  (scalar)
%             (Equispaced knots in [a,b] will be used)
%   order:  Spline order                (scalar; 4 is cubic)
%   lsm:    Smoothing parameter         (non-negative scalar)
%
%OUTPUT:
%   theta:  Estimated spline coefficients ((nk+order) x 1)
%   sigma2: Estimated variance            (scalar)
%   AIC:    AIC criterion               (scalar)
%   BIC:    BIC criterion               (scalar)
%
%
%External programs called: BSPL
%
% Version: April 2016

if nargin==7
    lsm = 0;
end

indiv = unique(id);
n = length(indiv);

knots = linspace(a,b,nk+2);
p = nk+order;

B = cell(n,1);
for i = 1:n
    if isnumeric(id)
        iab = id==indiv(i) & t>=a & t<=b;
    else
        iab = strcmp(id,indiv(i)) & t>=a & t<=b;
    end
    if any(iab)
        B{i} = bspl(t(iab),order,knots,0);
    end
end

if order>=3
    tg = linspace(a,b,300);
    B2 = bspl(tg,order,knots,2);
    Omega = (B2'*B2)*(tg(2)-tg(1));
else
    Omega = zeros(p,p);
end

A1 = zeros(p,p);
A2 = zeros(p,1);
N = 0;
for i = 1:n
    if isnumeric(id)
        iab = id==indiv(i) & t>=a & t<=b;
    else
        iab = strcmp(id,indiv(i)) & t>=a & t<=b;
    end
    if any(iab)
        xx = x(iab);
        m = length(xx);
        A1 = A1 + B{i}'*B{i};
        A2 = A2 + B{i}'*xx;
        N = N + m;
    end
end
theta = (A1+n*lsm*Omega)\A2;
ss = 0;
for i = 1:n
    if isnumeric(id)
        iab = id==indiv(i) & t>=a & t<=b;
    else
        iab = strcmp(id,indiv(i)) & t>=a & t<=b;
    end
    if any(iab)
        xx = x(iab);
        ss = ss + norm(xx-B{i}*theta)^2;
    end
end
sigma2 = ss/N;

loglik = 0;
for i = 1:n
    if isnumeric(id)
        iab = id==indiv(i) & t>=a & t<=b;
    else
        iab = strcmp(id,indiv(i)) & t>=a & t<=b;
    end
    if any(iab)
        xx = x(iab);
        m = length(xx);
        S0 = sigma2*eye(m);
        loglik = loglik - 0.5*m*log(2*pi) - 0.5*log(det(S0)) ...
            - 0.5*(xx-B{i}*theta)'*(S0\(xx-B{i}*theta));
    end
end

p_gen = trace((A1+n*lsm*Omega)\A1);  % This is just p if lsm=0
df = p_gen + 1;
AIC = -loglik + df;
BIC = -loglik + (log(n)/2)*df;
