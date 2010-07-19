clear;clc;
tic
% parSet(rho,lambdam,lambdah,r,wage,max_edu,muh,mum,theta,sigma,a)
% default parameter 
% par = parSet(0.1,0.2,0.2,0.05,0.15,0.001,0.02,0.1,1,(0:0.1:20)');
par = parSet(0.05,0.3,0.2,0.08,0.068,2,0.001,0.02,0.1,1,(0:0.5:100)',(0:3)');

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
v0 = repmat(0,length(par.a),1);
h0 = repmat(0,length(par.a),1);

aMatrix=repmat(par.a,1,length(par.a));


c=(1+par.r)./(1-par.mum).*aMatrix -aMatrix';  % aMatrix <-> at  ; aMatrix' <-> at+1 
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


%% morbid and working stage (method 1)

% method 1
% loop over different value of education. Easy to code but difficult to
% extend

for j=length(par.edu):-1:1
    v0 = repmat(0,length(par.a),1);
    h0 = repmat(0,length(par.a),1);

    aMatrix=repmat(par.a,1,length(par.a));

    c=(1+par.r)./(1-par.mum).*aMatrix -aMatrix' + wage(par.edu(j),par);
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

    v_mw1(:,j) = v0;
    h_mw1(:,j) = h0;
    clear h0 v0 v1 c aMatrix i;
end;


%% morbid and working stage (method 2)

% 1st dimension: asset (state variable)
% 2nd dimension: education (state variable)
% 3rd dimension: at+1 (control state variable)
% so the dimension would be [length(par.a) length(par.edu) length(par.a)]
v0 = repmat(0,[length(par.a) length(par.edu)]); 
h0 = repmat(0,[length(par.a) length(par.edu)]);

% create the 3rd dimension for control state varible. Since it is still
% par.a, so we use repmat(par.a,...) to createthe matrix
at_Matrix = permute(repmat(par.a,1,length(par.a)),[1 3 2] );
at_plus1_matrix = permute(repmat(par.a,1,length(par.a)),[2 3 1] );

% create an edu with correct dimensionality 
edu = repmat(par.edu',[ length(par.a) 1 length(par.a) ] );


c= repmat( (1+par.r)./(1-par.mum).*at_Matrix -at_plus1_matrix, [1 length(par.edu) 1 ]  ) + wage(edu,par); 
c(c<0)=nan;

for i=1:maxInteration
    %
    vt_plus1 = repmat(permute(v0,[3 2 1]), [length(par.a) 1 1]) ;
    [v_t h0]=max(utility(c)-par.lambdam +par.beta.*(1-par.mum).* vt_plus1 ,[],3  ); 
       
    v1 = max(v_t,repmat(v_mr,[1 length(par.edu)]));
    error = abs(v0-v1);
    error = max(error(:));
    disp(['iteration '  int2str(i)  ': '  num2str(error)]);
    if error < maxError;
        numberOfIteration_mw = i;
        break;
    end;
    v0=v1;
end;

v_mw2 = v0;
h_mw2 = h0;
clear h0 v0 v1 c aMatrix i;


%% morbid and schooling stage (method 1)

% method 1
% loop over different value of education. Easy to code but difficult to
% extend
tic;
for j=length(par.edu):-1:1
    v0 = repmat(0,length(par.a),1);
    h0 = repmat(0,length(par.a),1);

    aMatrix=repmat(par.a,1,length(par.a));

    c=(1+par.r)./(1-par.mum).*aMatrix -aMatrix';
    c(c<0)=nan;

    for i=1:maxInteration
        if par.edu(j)==par.edu(end)
            vt_plus1 = (repmat(v0',length(par.a),1));
        else 
            vt_plus1 = (repmat(v_ms1(:,j+1)',length(par.a),1));
        end

        [v_t h0]=max(utility(c)-par.lambdam +par.beta.*(1-par.mum).*vt_plus1 ,[],2);
        v1 = max(v_t,v_mw1(:,j) );
        error = max(abs(v0-v1));
        disp(['iteration '  int2str(i)  ': '  num2str(error)]);
        if error < maxError;
            numberOfIteration_mw = i;
            break;
        end;
        v0=v1;
        
    end;

    v_ms1(:,j) = v0;

    h_ms1(:,j) = h0;

    clear h0 v0 v1 c aMatrix i;
end;
toc;

%% morbid and schooling stage (method 2)

v0 = repmat(0,[length(par.a) length(par.edu)]); 
h0 = repmat(0,[length(par.a) length(par.edu)]);

at_Matrix = permute(repmat(par.a,1,length(par.a)),[1 3 2] );
at_plus1_matrix = permute(repmat(par.a,1,length(par.a)),[2 3 1] );
edu = repmat(par.edu,[ length(par.a) 1 length(par.a) ] );

c= repmat( (1+par.r)./(1-par.mum).*at_Matrix -at_plus1_matrix, [1 length(par.edu) 1 ]) ; 
c(c<0)=nan;


for i=1:maxInteration
    
    vt_plus1 = repmat(permute(v0(:,[2:end end]),[3 2 1]), [length(par.a) 1 1]) ;
    
    [v_t h0]=max(utility(c)-par.lambdam +par.beta.*(1-par.mum).*vt_plus1  ,[],3  );
       
    v1 = max(v_t,v_mw2);
    
    error = abs(v0-v1);
    error = max(error(:));

    disp(['iteration '  int2str(i)  ': '  num2str(error)]);
    if error < maxError;
        numberOfIteration_mw = i;
        break;
    end;
    v0=v1;
end;

v_ms2 = v0;
h_ms2 = h0;
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


%% healthy and working stage (method 1)

% method 1
% loop over different value of education. Easy to code but difficult to
% extend

for j=length(par.edu):-1:1
    v0 = repmat(0,length(par.a),1);
    h0 = repmat(0,length(par.a),1);

    aMatrix=repmat(par.a,1,length(par.a));

    c=(1+par.r)./(1-par.muh).*aMatrix -aMatrix' + wage(par.edu(j),par);
    c(c<0)=nan;

    for i=1:maxInteration
        [v_t h0]=max(utility(c)-par.lambdah +par.beta.*(1-par.muh-par.theta).*(repmat(v0',length(par.a) ,1)) ...
 +   par.beta * par.theta *   (repmat(v_mw1(:,j)',length(par.a) ,1))  ,[],2) ;
        v1 = max(v_t,v_hr);
        error = max(abs(v0-v1));
        disp(['iteration '  int2str(i)  ': '  num2str(error)]);
        if error < maxError;
            numberOfIteration_mw = i;
            break;
        end;
        v0=v1;
    end;

    v_ms1(:,j) = v0;
    h_ms1(:,j) = h0;
    clear h0 v0 v1 c aMatrix i;
end;


%% healthy and working stage (method 2)

% 1st dimension: asset (state variable)
% 2nd dimension: education (state variable)
% 3rd dimension: at+1 (control state variable)
% so the dimension would be [length(par.a) length(par.edu) length(par.a)]
v0 = repmat(0,[length(par.a) length(par.edu)]); 
h0 = repmat(0,[length(par.a) length(par.edu)]);

% create the 3rd dimension for control state varible. Since it is still
% par.a, so we use repmat(par.a,...) to createthe matrix
at_Matrix = permute(repmat(par.a,1,length(par.a)),[1 3 2] );
at_plus1_matrix = permute(repmat(par.a,1,length(par.a)),[2 3 1] );

% create an edu with correct dimensionality 
edu = repmat(par.edu',[ length(par.a) 1 length(par.a) ] );


c= repmat( (1+par.r)./(1-par.muh).*at_Matrix -at_plus1_matrix, [1 length(par.edu) 1 ]  ) + wage(edu,par); 
c(c<0)=nan;

for i=1:maxInteration
    %
    vt_plus1 = repmat(permute(v0,[3 2 1]), [length(par.a) 1 1]) ;
    vt_mw_plus1 = repmat(permute(v_mw2,[3 2 1]), [length(par.a) 1 1]) ;
    
    [v_t h0]=max(utility(c)-par.lambdah +par.beta.*(1-par.muh-par.theta).* vt_plus1 + par.beta*par.theta* vt_mw_plus1 ,[],3  ); 
       
    v1 = max(v_t,repmat(v_hr,[1 length(par.edu)]));
    error = abs(v0-v1);
    error = max(error(:));
    disp(['iteration '  int2str(i)  ': '  num2str(error)]);
    if error < maxError;
        numberOfIteration_mw = i;
        break;
    end;
    v0=v1;
end;

v_ms2 = v0;
h_ms2 = h0;
clear h0 v0 v1 c aMatrix i;


v_ms1-v_ms2


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
