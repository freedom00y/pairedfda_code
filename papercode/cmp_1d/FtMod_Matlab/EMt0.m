function [theta,sigma2,w,AIC,BIC] = ...
    EMt0(id,t,x,a,b,theta,sigma2,nk,order,nu,maxit,dsp)
%function [theta,sigma2,w,AIC,BIC] = ...
%   EMt0(id,t,x,a,b,theta,sigma2,nk,order,nu,maxit,dsp)
%Semiparametric t-estimators of the mean (via EM algorithm)
%
%NOTE: since the time grids may be sparse and irregular, the data is input 
% as a concatenated vector rather than a matrix; a vector of labels is used
% to identify the individuals
%
%INPUT:
%   id:     Individual labels           (N x 1)
%   t:      Time grid                   (N x 1)
%   x:      Observations                (N x 1)
%   a:      Time range, lower bound     (scalar)
%   b:      Time range, upper bound     (scalar)
%   theta:  Initial spline coefficients ((nk+order) x 1)
%             (Enter [] for default)
%   sigma2: Initial variance estimator  (scalar)
%             (Enter [] for default)
%   nk:     Number of knots (excluding endpoints) (scalar)
%             (Equispaced knots in [a,b] will be used)
%   order:  Spline order                (scalar)
%   nu:     t model degrees of freedom  (scalar)
%             (Use nu = 1 for Cauchy estimators)
%   maxit:  Max. number of iterations   (scalar)
%   dsp:    Display on each iteration   ('on'/'off')
%
%OUTPUT:
%   theta:  Estimated spline coefficients ((nk+order) x 1)
%   sigma2: Estimated variance            (scalar)
%   w:      Data weights (nu+m)/(nu+s)  (n x 1)
%   AIC:    AIC criterion               (scalar)
%   BIC:    BIC criterion               (scalar)
%
%External programs called: BSPL
%
% Version: May 2010

indiv = unique(id);
n = length(indiv);

knots = linspace(a,b,nk+2);
p = nk+order;

if isempty(theta)
    theta = zeros(p,1);
end
if isempty(sigma2)
    sigma2 = var(x);
end

B = cell(n,1);
for i = 1:n
    iab = id==indiv(i) & t>=a & t<=b;
    if any(iab)
        B{i} = bspl(t(iab),order,knots,0);
    end
end

w = zeros(n,1);
err = 1;
iter = 0;
while err>1e-3 && iter<maxit
    sigma20 = sigma2;
    iter = iter + 1;
    A1 = zeros(p,p);
    A2 = zeros(p,1);
    ss = 0;
    N = 0;
    for i = 1:n
        iab = id==indiv(i) & t>=a & t<=b;
        if any(iab)
            xx = x(iab);
            m = length(xx);
            d2 = (xx-B{i}*theta)'*(xx-B{i}*theta)/sigma2;
            Eu = (nu+m)/(nu+d2);
            w(i) = Eu;
            A1 = A1 + Eu*(B{i}'*B{i});
            A2 = A2 + Eu*B{i}'*xx;
            ss = ss + Eu*norm(xx-B{i}*theta)^2;
            N = N + m;
        end
    end
    theta = A1\A2;
    sigma2 = ss/N;
    err = abs(sigma2/sigma20-1);
    if strcmp(dsp,'on')
        disp(['Iter: ' num2str(iter) ...
            ', S^2 = ' num2str(sigma2) ', Error = ' num2str(err)])
    end
end

loglik = 0;
for i = 1:n
    iab = id==indiv(i) & t>=a & t<=b;
    if any(iab)
        xx = x(iab);
        m = length(xx);
        S0 = sigma2*eye(m);
        loglik = loglik - 0.5*m*log(nu*pi) + gammaln(0.5*(nu+m)) ...
            - gammaln(0.5*nu) - 0.5*log(det(S0)) ...
            - 0.5*(nu+m)*log(1+(xx-B{i}*theta)'*(S0\(xx-B{i}*theta))/nu);
    end
end
df = p + 1;
AIC = -loglik + df;
BIC = -loglik + (log(n)/2)*df;
