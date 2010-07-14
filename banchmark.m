clear;clc;
% parSet(rho,lambdam,lambdah,r,w,muh,mum,theta,sigma,a)
% default parameter 
% par = parSet(0.1,0.2,0.2,0.05,0.15,0.001,0.02,0.1,1,(0:0.1:20)');
par = parSet(0.04,0.3,0.2,0.05,0.15,0.001,0.02,0,1,(0:0.005:5)');

maxInteration = 10000;
maxError = 10^-4;

numberOfIteration_mr = 0;

% Define the utility function
if par.sigma==1
    utility = @log;
else
    utility = @(x)(x.^(1-1./par.sigma)-1)./(1-1./par.sigma);
end

% initial value of value function.
v0=repmat(0,length(par.a),1);
h0=repmat(0,length(par.a),1);

aMatrix=repmat(par.a,1,length(par.a));

c=(1+par.r)./(1-par.mum).*aMatrix -aMatrix';
% c=max(c,0);
c(c<0)=nan;


for i=1:maxInteration
    [v1 h0]=max(utility(c)+par.beta.*(1-par.mum).*(repmat(v0',length(par.a),1)),[],2);
    error = max(abs(v0-v1));
    disp(['iteration '  int2str(i)  ': '  num2str(error)]);
    if error < maxError;
        numberOfIteration_mr = i;
        break;
    end
    v0=v1;
end

h0 = par.a(h0);
h1 = (1+par.r)./(1-par.mum).*par.a -h0;
figure;
% subplot(1,2,1)
% plot(par.a,h0)
% subplot(1,2,2)
% plot(par.a,h1)

ttt = par.a .* (1+par.r) ./ (1-par.mum) .* (1-par.beta.*(1-par.mum) );
plot(par.a,ttt,par.a,h1)

v_mr = v0;
clear v v1 c aMatrix i;
% plot(par.a,v_mr);

%  figure
%  compute the slope;
% vpi_mr = gradient(v_mr,par.a);
%  
% 
% % figure;
% subplot(2,1,1)
% plot(par.a,v_mr,par.a,dis_sol(par));
% e= v_mr-dis_sol(par);
% subplot(2,1,2)
% 
% adjError= e./dis_sol(par);
% plot(par.a,e);
% 
% mean(abs(e(2:end)))
% mean(abs(adjError(2:end)))
% max(abs(adjError(2:end)))
