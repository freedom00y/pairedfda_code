function [theta,etas,lambdas,sigma2,y,xhat,w,AIC,BIC] = ...
    EMt(id,t,x,a,b,theta,etas,lambdas,sigma2,nk,order,nu,maxiter,dsp)
%function [theta,etas,lambdas,sigma2,y,xhat,w,AIC,BIC] = ...
%   EMt(id,t,x,a,b,theta,etas,lambdas,sigma2,nk,order,nu,maxiter,dsp)
%Semiparametric t-estimators of the mean and principal components
%
%NOTE 1: since the time grids may be sparse and irregular, the data is 
% input as a concatenated vector rather than a matrix; a vector of labels 
% is used to identify the individuals.
%
%NOTE 2: although the program accepts arbitrary initial estimators, it is 
% best used in a stepwise manner, i.e.: fit a mean + error model first 
% (using program EMt0), then fit a one-component model using the output 
% THETA from EMt0 as initial THETA, and so on.
%
%INPUT:
%   id:     Individual labels           (N x 1)
%   t:      Time grid                   (N x 1)
%   x:      Observations                (N x 1)
%   a:      Time range, lower bound     (scalar)
%   b:      Time range, upper bound     (scalar)
%   theta:  Initial mean coefficients   ((nk+order) x 1)
%            (See Note 2 above)
%   etas:   Initial PC coefficients     ((nk+order) x q)
%            (See Note 2; we recommend etas = [etas0, ones((nk+order),1)]
%            where etas0 is the output from previous (q-1)-component model)
%   lambdas: Initial eigenvalues        (q x 1)
%            (See Note 2; we recommend lambdas = [lambdas0; lambdas0(end)/2]
%            where lambdas0 is output from previous (q-1)-component model)
%   sigma2: Initial variance estimator  (scalar)
%             (See Note 2 above)
%   nk:     Number of knots (excluding endpoints)  (scalar)
%             (Equispaced knots in [a,b] will be used)
%   order:  Spline order                (scalar)
%   nu:     t model degrees of freedom  (scalar)
%             (Use nu = 1 for Cauchy estimators)
%   maxiter: Max. number of iterations  (scalar)
%   dsp:    Display on each iteration   ('on'/'off')
%
%OUTPUT:
%   theta:  Estimated mean coefficients ((nk+order) x 1)
%   etas:   Estimated PC coefficients   ((nk+order) x q)
%   lambdas: Estimated eigenvalues      (q x 1)
%   sigma2: Estimated variance          (scalar)
%   y:      Predicted component scores  (q x n)
%   xhat:   Predicted values of X       (N x 1)
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
q = length(lambdas);
if q==0
    error('For a mean + error model use EMt0')
end
if size(etas,1)~=p || size(etas,2)~=q
    error('Wrong dimension of ETAS')
end
if length(theta)~=p
    error('Wrong dimension of THETA')
end

H = etas*diag(sqrt(lambdas));
grN = max(200,10*nk);
B0 = bspl(linspace(a,b,grN+1),order,knots,0);
J0 = B0'*B0*(b-a)/grN;

B = cell(n,1);
for i = 1:n
    iab = id==indiv(i) & t>=a & t<=b;
    if any(iab)
        B{i} = bspl(t(iab),order,knots,0);
    end
end

y = zeros(q,n);
xhat = NaN(size(x));
w = zeros(n,1);
err = 1;
iter = 0;
while err>1e-3 && iter<maxiter
    [U,S] = svd(H'*J0*H);   %Don't delete the S!!!!!!!
    H = H*U;
    sigma20 = sigma2;
    iter = iter + 1;
    A1 = zeros(p,p);
    A2 = zeros(p,1);
    A3 = zeros(p*q,p*q);
    A4 = zeros(p*q,1);
    ss = 0;
    N = 0;
    for i = 1:n
        iab = id==indiv(i) & t>=a & t<=b;
        if any(iab)
            xx = x(iab);
            m = length(xx);
            V0_1 = inv(eye(q) + H'*(B{i}'*B{i})*H/sigma2);
            S0_1 = eye(m)/sigma2 - B{i}*H*V0_1*H'*B{i}'/sigma2^2;
            y(:,i) = H'*B{i}'*S0_1*(xx-B{i}*theta);
            d2 = (xx-B{i}*theta)'*S0_1*(xx-B{i}*theta);
            Eu = (nu+m)/(nu+d2);
            w(i) = Eu;
            Ezz = V0_1 + Eu*y(:,i)*y(:,i)';
            Euz = Eu*y(:,i);
            A1 = A1 + Eu*(B{i}'*B{i});
            A2 = A2 + Eu*B{i}'*(xx-B{i}*H*y(:,i));
            A3 = A3 + kron(Ezz,B{i}'*B{i});
            A4 = A4 + Eu*kron(y(:,i),B{i}')*(xx-B{i}*theta);
            ss = ss + Eu*norm(xx-B{i}*theta)^2 + trace(B{i}*H*Ezz*H'*B{i}') ...
                - 2*Euz'*H'*B{i}'*(xx-B{i}*theta);
            N = N + m;
            xhat(iab) = B{i}*theta + B{i}*H*y(:,i);
        end
    end
    theta = A1\A2;
    H = reshape(A3\A4,p,q);
    stdy = std(y,0,2);
    H = H*diag(stdy);
    y = diag(1./(stdy.*(stdy>0)+(stdy==0)))*y;
    sigma2 = ss/N;
    err = abs(sigma2/sigma20-1);
    if strcmp(dsp,'on')
        disp(['Iter: ' num2str(iter) ...
            ', S^2 = ' num2str(sigma2) ', Error = ' num2str(err)])
    end
end

[U,S] = svd(H'*J0*H);
lambdas = diag(S);
invsqlmb = zeros(size(lambdas));
invsqlmb(lambdas>eps) = 1./sqrt(lambdas(lambdas>eps));
etas = H*U*diag(invsqlmb);
y = U'*y;

loglik = 0;
for i = 1:n
    iab = id==indiv(i) & t>=a & t<=b;
    if any(iab)
        xx = x(iab);
        m = length(xx);
        S0 = sigma2*eye(m) + B{i}*etas*diag(lambdas)*etas'*B{i}';
        loglik = loglik - 0.5*m*log(nu*pi) + gammaln(0.5*(nu+m)) ...
            - gammaln(0.5*nu) - 0.5*log(det(S0)) ...
            - 0.5*(nu+m)*log(1+(xx-B{i}*theta)'*(S0\(xx-B{i}*theta))/nu);
    end
end
df = p + p*q + q + 1 - (q + q*(q-1)/2);
AIC = -loglik + df;
BIC = -loglik + (log(n)/2)*df;
