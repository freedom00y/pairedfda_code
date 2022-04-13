function [msqres,qval,m] = residuals(id,t,x,a,b,xhat,sigma2)
%function [msqres,qval,m] = residuals(id,t,x,a,b,xhat,sigma2)

indiv = unique(id);
n = length(indiv);
sqres = zeros(n,1);
msqres = zeros(n,1);
qval = zeros(n,1);
m = zeros(n,1);
for i = 1:n
    iab = id==indiv(i) & t>=a & t<=b;
    xx = x(iab);
    xxhat = xhat(iab);
    m(i) = length(xx);
    sqres(i) = sum((xx-xxhat).^2);
    msqres(i) = mean((xx-xxhat).^2);
    qval(i) = chi2cdf(sqres(i)/sigma2,m(i));
end