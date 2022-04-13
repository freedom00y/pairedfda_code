function Normalfit = fit_normal(id,t,x,a,b,nk,order,d,dsp,lsm)
%function Normalfit = fit_normal(id,t,x,a,b,nk,order,d,dsp,lsm)
% Fits reduced-rank Normal models sequentially, up to dimension D.
% Input arguments are as in EMnormal.m.
% Output is a struct with field names as in the output of EMnormal.m.
%   (NOTE: first element in each output field corresponds to d=0)
%
% External programs called: EMnormal0, EMnormal
%
% Version: April 2016

if nargin == 9
    lsm = 0;
end

theta = cell(1,d+1);
etas = cell(1,d+1);
lambdas = cell(1,d+1);
sigma2 = zeros(1,d+1);
y = cell(1,d+1);
xhat = cell(1,d+1);
AIC = zeros(1,d+1);
BIC = zeros(1,d+1);

p = nk+order;

if strcmp(dsp,'on')
    disp('Fitting initial mean-only model')
    disp(' ')
    pause(1)
end
[theta{1},sigma2(1),AIC(1),BIC(1)] = EMnormal0(id,t,x,a,b,nk,order,lsm);
if d>=1,
    if strcmp(dsp,'on')
        disp(' ')
        disp(['Fitting ' num2str(1) '-component model'])
        disp(' ')
        pause(1)
    end
    [theta{2},etas{2},lambdas{2},sigma2(2),y{2},xhat{2},AIC(2),BIC(2)] = ...
        EMnormal(id,t,x,a,b,theta{1},ones(p,1),sigma2(1),sigma2(1)/10,nk,...
        order,300,dsp,lsm);
    for i = 3:d+1
        if strcmp(dsp,'on')
            disp(' ')
            disp(['Fitting ' num2str(i-1) '-component model'])
            disp(' ')
            pause(1)
        end
        [theta{i},etas{i},lambdas{i},sigma2(i),y{i},xhat{i},AIC(i),BIC(i)] = ...
            EMnormal(id,t,x,a,b,theta{i-1},[etas{i-1},ones(p,1)],...
            [lambdas{i-1};lambdas{i-1}(end)/10],sigma2(i-1)/10,nk,order,...
            300,dsp,lsm);
    end
end

Normalfit.theta = theta;
Normalfit.etas = etas;
Normalfit.lambdas = lambdas;
Normalfit.sigma2 = sigma2;
Normalfit.y = y;
Normalfit.xhat = xhat;
Normalfit.AIC = AIC;
Normalfit.BIC = BIC;