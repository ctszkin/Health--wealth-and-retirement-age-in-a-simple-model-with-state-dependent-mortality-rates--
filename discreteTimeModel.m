clear;clc;
tic
% parSet(rho,lambdam,lambdah,r,w,muh,mum,theta,sigma,a)
% default parameter 
% par = parSet(0.1,0.2,0.2,0.05,0.15,0.001,0.02,0.1,1,(0:0.1:20)');
par = parSet(0.05,0.3,0.2,0.08,2.127,0.001,0.02,0.1,1,(0:0.5:100)');
educ=0:1;
maxInteration = 10000;
maxError = 10^-4;

numberOfIteration_mr = 0;
numberOfIteration_hr = 0;
numberOfIteration_mw = 0;
numberOfIteration_hw = 0;


% Define the utility function
if par.sigma==1
    utility = @log;
else
    utility = @(x)(x.^(1-1./par.sigma)-1)./(1-1./par.sigma);
end

%% morbid & retired stage

% initial value of value function.
v0=repmat(0,length(par.a),1);
h0 = repmat(0,length(par.a),1);

aMatrix=repmat(par.a,1,length(par.a));


c=(1+par.r)./(1-par.mum).*aMatrix -aMatrix';
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

v_mr(:,1) = v0;
h_mr(:,1) = h0;
clear v0 v1 c aMatrix i h0;

%% morbid and working stage

% method 1
% for edu=0:1
%     v0 = repmat(0,length(par.a),1);
%     h0 = repmat(0,length(par.a),1);
% 
%     aMatrix=repmat(par.a,1,length(par.a));
% 
% 
%     c=(1+par.r)./(1-par.mum).*aMatrix -aMatrix' + wage(edu);
%     c(c<0)=nan;
% 
%     for i=1:maxInteration
%         [v_t h0]=max(utility(c)-par.lambdam +par.beta.*(1-par.mum).*(repmat(v0',length(par.a),1)),[],2);
%         v1 = max(v_t,v_mr);
%         error = max(abs(v0-v1));
%         disp(['iteration '  int2str(i)  ': '  num2str(error)]);
%         if error < maxError;
%             numberOfIteration_mw = i;
%             break;
%         end;
%         v0=v1;
%     end;
% 
%     v_mw(:,edu+1) = v0;
%     h_mw(:,edu+1) = h0;
%     clear h0 v0 v1 c aMatrix i;
% end;

% method 2 ( easiser to extend )

v0 = repmat(0,[length(par.a) length(educ)]); 
h0 = repmat(0,[length(par.a) length(educ)]);

at_Matrix = permute(repmat(par.a,1,length(par.a)),[1 3 2] );
at_plus1_matrix = permute(repmat(par.a,1,length(par.a)),[2 3 1] );
edu = repmat(educ,[ length(par.a) 1 length(par.a) ] );


c= repmat( (1+par.r)./(1-par.mum).*at_Matrix -at_plus1_matrix, [1 2 1 ]  ) + wage(edu); % aMatrix is at+1
c(c<0)=nan;

for i=1:maxInteration
    [v_t h0]=max(utility(c)-par.lambdam +par.beta.*(1-par.mum).* repmat(permute(v0,[3 2 1]), [length(par.a) 1 1]) ,[],3  ); %v0' at+1
       
    v1 = max(v_t,repmat(v_mr,[1 2]));
    error = max(abs(v0-v1));
    disp(['iteration '  int2str(i)  ': '  num2str(error)]);
    if error < maxError;
        numberOfIteration_mw = i;
        break;
    end;
    v0=v1;
end;

v_mw = v0;
h_mw = h0;
clear h0 v0 v1 c aMatrix i;





%% morbid and schooling stage


v0 = repmat(0,length(par.a),2);
h0 = repmat(0,length(par.a),2);

aMatrix=repmat(par.a,1,length(par.a));

c=(1+par.r)./(1-par.mum).*aMatrix -aMatrix' + wage(edu);
c(c<0)=nan;

for i=1:maxInteration
    [v_t h0]=max(utility(c)-par.lambdam +par.beta.*(1-par.mum).*(repmat(v0',length(par.a),1)),[],2);
    v1 = max(v_t,v_mr);
    error = max(abs(v0-v1));
    disp(['iteration '  int2str(i)  ': '  num2str(error)]);
    if error < maxError;
        numberOfIteration_mw = i;
        break;
    end;
    v0=v1;
end;

v_mw(:,edu+1) = v0;
h_mw(:,edu+1) = h0;
clear h0 v0 v1 c aMatrix i;




%% healthy and retired

v0 = repmat(0,length(par.a),1);
h0 = repmat(0,length(par.a),1);


aMatrix=repmat(par.a,1,length(par.a));

c=(1+par.r)./(1-par.muh).*aMatrix -aMatrix';
c(c<0)=nan;

for i=1:maxInteration
    [v1 h0]=max(utility(c)+par.beta.*(1-par.muh-par.theta).*(repmat(v0',length(par.a),1))+par.beta.*par.theta.*(repmat(v_mr',length(par.a),1)) ,[],2);
    error = max(abs(v0-v1));
    disp(['iteration '  int2str(i)  ': '  num2str(error)]);
    if error < maxError;
        numberOfIteration_hr = i;
        break;
    end
    v0=v1;
end

v_hr = v0;
h_hr = h0;
clear h0 v0 v1 c aMatrix i;


%% healthy working
v0 = repmat(0,length(par.a),1);
h0 = repmat(0,length(par.a),1);


aMatrix=repmat(par.a,1,length(par.a));

c=(1+par.r)./(1-par.muh).*aMatrix -aMatrix' + par.w;
c(c<0)=nan;

for i=1:maxInteration
    [v_t h0]=max(utility(c)-par.lambdah +par.beta.*(1-par.muh-par.theta).*(repmat(v0',length(par.a),1))+ par.beta.*par.theta.*(repmat(v_mw',length(par.a),1)),[],2);
    v1 = max(v_t,v_hr);
    error = max(abs(v0-v1));
    disp(['iteration '  int2str(i)  ': '  num2str(error)]);
    if error < maxError;
        numberOfIteration_hw = i;
        break;
    end
    v0=v1;
end

v_hw = v0;
h_hw = h0;
clear h0 v0 v1 v_t c aMatrix i;


%% ms %%


v0 = repmat(0,length(par.a),1);
h0 = repmat(0,length(par.a),1);


aMatrix=repmat(par.a,1,length(par.a));

c=(1+par.r)./(1-par.mum).*aMatrix -aMatrix' + par.w;
c(c<0)=nan;

for i=1:maxInteration
    [v_t h0]=max(utility(c)-par.lambdam +par.beta.*(1-par.mum).*(repmat(v0',length(par.a),1)),[],2);
    v1 = max(v_t,v_mr);
    error = max(abs(v0-v1));
    disp(['iteration '  int2str(i)  ': '  num2str(error)]);
    if error < maxError;
        numberOfIteration_mw = i;
        break;
    end
    v0=v1;
end

v_mw = v0;
h_mw = h0;
clear h0 v0 v1 c aMatrix i;


%% end of ms %%




% plot(par.a,v_mw,par.a,v_mr,par.a,v_hr,par.a,v_hw);
legend('morbid working','morbid retired','healthy retired','healthy working','Location','southeast');

% compute the critial value;

 [temp critial_pt_h] = max(v_hr==v_hw);
 [temp critial_pt_m] = max(v_mr==v_mw);
 
%  compute the slope;
vpi_mr = gradient(v_mr,par.a);
vpi_mw = gradient(v_mw,par.a);
vpi_hr = gradient(v_hr,par.a);
vpi_hw = gradient(v_hw,par.a);
 
% difference in slope
d_vpi_m = vpi_mr-vpi_mw;
d_vpi_h = vpi_hr-vpi_hw;

 [temp critial_pt_m_slope] = max(vpi_mr==vpi_mw);
 [temp critial_pt_h_slope] = max(vpi_hr==vpi_hw);
 
 
%  Compute the policy function
%  require critial_pt_h and critial_pt_m exist(~=1)



api_mr = par.a(h_mr);
c_mr = (1+par.r)./(1-par.mum).*par.a -api_mr;

api_mw = par.a(h_mw);
c_mw = (1+par.r)./(1-par.mum).*par.a - api_mw +par.w;
c_mw(critial_pt_m:end) = c_mr(critial_pt_m:end);
api_mw(critial_pt_m:end) = api_mr(critial_pt_m:end);

api_hr = par.a(h_hr);
c_hr = (1+par.r)./(1-par.muh).*par.a -api_hr;

api_hw = par.a(h_hw);
c_hw = (1+par.r)./(1-par.muh).*par.a - api_hw +par.w;
c_hw(critial_pt_h:end) = c_hr(critial_pt_h:end);
api_hw(critial_pt_h:end) = api_hr(critial_pt_h:end);


%   summary
toc;clc;
disp('Summary');
disp([' ']);

disp(['Time used: ' num2str(toc)]);
disp(['Number of iteration']);
disp(['morbid retired: ' num2str(numberOfIteration_mr) ]);
disp(['morbid working: ' num2str(numberOfIteration_mw)]);
disp(['healthy retired: ' num2str(numberOfIteration_hr)]);
disp(['healthy working: ' num2str(numberOfIteration_hw)]);

disp([' ']);
disp(['Critial point for retirement'])
disp(['Healthy case: ' num2str(par.a(critial_pt_h)) ]);
disp(['Morbid case: ' num2str(par.a(critial_pt_m)) ]);
 
disp([' ']);
disp(['Critial point for retirement(slope)'])
disp(['Healthy case: ' num2str(par.a(critial_pt_h_slope)) ]);
disp(['Morbid case: ' num2str(par.a(critial_pt_m_slope)) ]);

save valueFn v_mr v_hr v_mw v_hw api_mr api_hr api_mw api_hw c_mr c_hr c_mw c_hw critial_pt_m critial_pt_h par; 


% subplot(1,2,1)
% plot(par.a,c_mr,par.a,c_mw)
% subplot(1,2,2)
% plot(par.a,c_hr,par.a,c_hw)
% 
% subplot(1,2,1)
% plot(par.a,api_mr,par.a,api_mw,par.a,par.a)
% subplot(1,2,2)
% plot(par.a,api_hr,par.a,api_hw,par.a,par.a)
% 
