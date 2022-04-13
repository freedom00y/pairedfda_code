function Tfit = fit_t(id,t,x,a,b,nk,order,nu,d,dsp)
%function Tfit = fit_t(id,t,x,a,b,nk,order,nu,d,dsp)
% Fits reduced-rank T models sequentially, up to dimension D.
% Input arguments are as in EMt.m.
% Output is a struct of dim D+1 with field names as in the output of EMt.m.
%   (NOTE: first element in each output field corresponds to d=0)
% Also a fine time grid TT and basis matrix BB are output (for plotting).
%
% External programs called: EMt0, EMt

theta = cell(1,d+1);
etas = cell(1,d+1);
lambdas = cell(1,d+1);
sigma2 = zeros(1,d+1);
y = cell(1,d+1);
xhat = cell(1,d+1);
w = cell(1,d+1);
AIC = zeros(1,d+1);
BIC = zeros(1,d+1);

p = nk+order;

disp('Fitting initial mean-only model')
disp(' ')
pause(1)
[theta{1},sigma2(1),w{1},AIC(1),BIC(1)] = EMt0(id,t,x,a,b,[],[],nk,order,nu,300,dsp);
if d>=1,
    disp(' ')
    disp(['Fitting ' num2str(1) '-component model'])
    disp(' ')
    pause(1)
    [theta{2},etas{2},lambdas{2},sigma2(2),y{2},xhat{2},w{2},AIC(2),BIC(2)] = ...
        EMt(id,t,x,a,b,theta{1},ones(p,1),sigma2(1),sigma2(1)/10,nk,order,nu,300,dsp);
    for i = 3:d+1
        disp(' ')
        disp(['Fitting ' num2str(i-1) '-component model'])
        disp(' ')
        pause(1)
        [theta{i},etas{i},lambdas{i},sigma2(i),y{i},xhat{i},w{i},AIC(i),BIC(i)] = ...
            EMt(id,t,x,a,b,theta{i-1},[etas{i-1},ones(p,1)],...
            [lambdas{i-1};lambdas{i-1}(end)/10],sigma2(i-1)/10,nk,order,nu,300,dsp);
    end
end

Tfit.theta = theta;
Tfit.etas = etas;
Tfit.lambdas = lambdas;
Tfit.sigma2 = sigma2;
Tfit.y = y;
Tfit.xhat = xhat;
Tfit.w = w;
Tfit.AIC = AIC;
Tfit.BIC = BIC;
Tfit.tt = linspace(a,b,200);
Tfit.BB = bspl(Tfit.tt,order,linspace(a,b,nk+2),0);