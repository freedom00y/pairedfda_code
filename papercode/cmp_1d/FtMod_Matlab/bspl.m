function y = bspl(x,k,t,r)
%function y = bspl(x,k,t,r)
%
%B-spline basis functions and their derivatives
%
%INPUT:
%   x   (m x 1 or 1 x m)    Input grid.
%   k   (scalar)            Spline order.
%   t   (n x 1 or 1 x n)    Knots, must be a strictly increasing sequence
%                             and must INCLUDE interval endpoints.
%   r	(scalar)            Order of derivative.
%
%OUTPUT:
%   y   (m x n+k-2)         Basis function (or derivative) values at X
%
% Version: May 2010

if nargin<4
    error('Not enough input arguments')
end
if size(t,1)>1
    t = t';
end

m = length(x);
n = length(t);
y = zeros(m,n+k-2);

if r==0
    
    tt = [t(1)*ones(1,k-1), t, t(n)*ones(1,k-1)];
    n = length(tt);
    b = zeros(1,k);
    dr = zeros(1,k-1);
    dl = zeros(1,k-1);
    for l = 1:m
        b(1) = 1;
        i = find(tt<=x(l),1,'last');
        if i==n, i = n-k; end
        for j = 1:k-1
            dr(j) = tt(i+j)-x(l);
            dl(j) = x(l)-tt(i+1-j);
            saved = 0;
            for r = 1:j
                term = b(r)/(dr(r)+dl(j+1-r));
                b(r) = saved + dr(r)*term;
                saved = dl(j+1-r)*term;
            end
            b(j+1) = saved;
        end
        y(l,i-k+1:i) = b;
    end
    
else
    
    tt = [repmat(t(1),1,k-2), t, repmat(t(n),1,k-2)];
    B = bspl(x,k-1,t,r-1);
    msp = ((k-1)./(ones(m,1)*(tt(k:n+2*(k-2))-tt(1:n+k-3)))).*B;
    y(:,1) =  - msp(:,1);
    y(:,2:n+k-3) = msp(:,1:n+k-4) - msp(:,2:n+k-3);
    y(:,n+k-2) =  msp(:,n+k-3);
    
end
