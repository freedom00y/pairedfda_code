clear
clc
addpath("FtMod_Matlab");
M=0;
for i=1:100
dirpath = strcat("data_s1/data",num2str(i),".mat");

load(dirpath)
id = data(:,1);
t  = data(:,2);
x  = data(:,3);
z  = data(:,4);
a  = 0;
b  = 1;
theta = [];
sigma2 = [];
nk = 10; % Need to check
order = 4;
nu = 2; % It didn't estimate the dof
maxit = 100;
dsp="off";

%% Analyze Y
% Mean curve
[thetaC0,sigma2C0,wC0] = EMt0(id,t,x,a,b,theta,sigma2,nk,order,nu,maxit,dsp);

tt = linspace(0,1,101);
knots1 = linspace(0,1,12);
B1 = bspl(tt,4,knots1,0);
%plot(tt,B1*theta)

% First PC
theta = thetaC0;
etas  = ones(nk+order,1);
lambdas = 2*sigma2C0;
sigma2 = sigma2C0;
[thetaC1,etasC1,lambdasC1,sigma2C1,yC1,xhatC1,wC1] = ...
EMt(id,t,x,a,b,theta,etas,lambdas,sigma2,nk,order,nu,maxit,dsp);
%subplot(121),plot(tt,B1*thetaC1)
%subplot(122),plot(tt,B1*etasC1)

% Second PC
[thetaC2,etasC2,lambdasC2,sigma2C2,yC2,xhatC2,wC2] = ...
EMt(id,t,x,a,b,thetaC1,[etasC1 ones(nk+order,1)], ...
[lambdasC1;lambdasC1(end)/2],sigma2C1,nk,order,nu,maxit,dsp);


mu = B1*thetaC2;
f  = B1*etasC2;
yhat = B1 * (thetaC2 + etasC2 * yC2);
%subplot(121),plot(tt,mu) %mean
%subplot(122),plot(tt,f) %PC

%% Analyze Z
% Mean curve
[thetaC0,sigma2C0,wC0] = EMt0(id,t,z,a,b,theta,sigma2,nk,order,nu,maxit,dsp);

% First PC
theta = thetaC0;
etas  = ones(nk+order,1);
lambdas = 2*sigma2C0;
sigma2 = sigma2C0;
[thetaC1,etasC1,lambdasC1,sigma2C1,yC1,xhatC1,wC1] = ...
EMt(id,t,z,a,b,theta,etas,lambdas,sigma2,nk,order,nu,maxit,dsp);
%subplot(121),plot(tt,B1*thetaC1)
%subplot(122),plot(tt,B1*etasC1)

% Second PC
[thetaC2,etasC2,lambdasC2,sigma2C2,yC2,xhatC2,wC2] = ...
EMt(id,t,z,a,b,thetaC1,[etasC1 ones(nk+order,1)], ...
[lambdasC1;lambdasC1(end)/2],sigma2C1,nk,order,nu,maxit,dsp);

nu = B1*thetaC2;
g  = B1*etasC2;
zhat = B1 * (thetaC2 + etasC2 * yC2);

%subplot(121),plot(tt,nu) %mean
%subplot(122),plot(tt,g) %PC


%% fitted curve
outpath = strcat("data_s1/res",num2str(i),".mat");
save(outpath,"mu","nu","f","g","yhat","zhat")
end











