clear;clc;
tic
% wage=1 
% return to educaiton 
% parSet(rho,lambdam,lambdah,r,wage,max_edu,muh,mum,theta,sigma,a,edu)
% default parameter 
% par = parSet(0.1,0.2,0.2,0.05,0.15,0.001,0.02,0.1,1,(0:0.1:20)');
par = parSet(0.05,0.3,0.2,0.08,0.2,2,0.001,0.02,0.1,1,(-20:0.5:200)',(0:3)');

maxInteration = 10000;
maxError = 10^-6;

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
%     disp(['iteration '  int2str(i)  ': '  num2str(error)]);
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
tic
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
%         disp(['iteration '  int2str(i)  ': '  num2str(error)]);
        if error < maxError;
            numberOfIteration_mw = i;
            break;
        end;
        v0=v1;
    end;

    v_mw(:,j) = v0;
    h_mw(:,j) = h0;
    clear h0 v0 v1 c aMatrix i;
end;
toc

% tic
% %% morbid and working stage (method 2)
% 
% % 1st dimension: asset (state variable)
% % 2nd dimension: education (state variable)
% % 3rd dimension: at+1 (control state variable)
% % so the dimension would be [length(par.a) length(par.edu) length(par.a)]
% v0 = repmat(0,[length(par.a) length(par.edu)]); 
% h0 = repmat(0,[length(par.a) length(par.edu)]);
% 
% % create the 3rd dimension for control state varible. Since it is still
% % par.a, so we use repmat(par.a,...) to createthe matrix
% at_Matrix = permute(repmat(par.a,1,length(par.a)),[1 3 2] );
% at_plus1_matrix = permute(repmat(par.a,1,length(par.a)),[2 3 1] );
% 
% % create an edu with correct dimensionality 
% edu = repmat(par.edu',[ length(par.a) 1 length(par.a) ] );
% 
% 
% c= repmat( (1+par.r)./(1-par.mum).*at_Matrix -at_plus1_matrix, [1 length(par.edu) 1 ]  ) + wage(edu,par); 
% c(c<0)=nan;
% 
% for i=1:maxInteration
%     %
%     vt_plus1 = repmat(permute(v0,[3 2 1]), [length(par.a) 1 1]) ;
%     [v_t h0]=max(utility(c)-par.lambdam +par.beta.*(1-par.mum).* vt_plus1 ,[],3  ); 
%        
%     v1 = max(v_t,repmat(v_mr,[1 length(par.edu)]));
%     error = abs(v0-v1);
%     error = max(error(:));
% %     disp(['iteration '  int2str(i)  ': '  num2str(error)]);
%     if error < maxError;
%         numberOfIteration_mw = i;
%         break;
%     end;
%     v0=v1;
% end;
% 
% v_mw2 = v0;
% h_mw2 = h0;
% clear h0 v0 v1 c aMatrix i;
% toc

%% morbid and schooling stage (method 1)

% method 1
% loop over different value of education. Easy to code but difficult to
% extend
tic
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
            vt_plus1 = (repmat(v_ms(:,j+1)',length(par.a),1));
        end

        [v_t h0]=max(utility(c)-par.lambdam +par.beta.*(1-par.mum).*vt_plus1 ,[],2);
        v1 = max(v_t,v_mw(:,j) );
        error = max(abs(v0-v1));
%         disp(['iteration '  int2str(i)  ': '  num2str(error)]);
        if error < maxError;
            numberOfIteration_mw = i;
            break;
        end;
        v0=v1;
        
    end;

    v_ms(:,j) = v0;

    h_ms(:,j) = h0;

    clear h0 v0 v1 c aMatrix i;
end;
toc

% %% morbid and schooling stage (method 2)
% 
% tic
% v0 = repmat(0,[length(par.a) length(par.edu)]); 
% h0 = repmat(0,[length(par.a) length(par.edu)]);
% 
% at_Matrix = permute(repmat(par.a,1,length(par.a)),[1 3 2] );
% at_plus1_matrix = permute(repmat(par.a,1,length(par.a)),[2 3 1] );
% edu = repmat(par.edu,[ length(par.a) 1 length(par.a) ] );
% 
% c= repmat( (1+par.r)./(1-par.mum).*at_Matrix -at_plus1_matrix, [1 length(par.edu) 1 ]) ; 
% c(c<0)=nan;
% 
% 
% for i=1:maxInteration
%     
%     vt_plus1 = repmat(permute(v0(:,[2:end end]),[3 2 1]), [length(par.a) 1 1]) ;
%     
%     [v_t h0]=max(utility(c)-par.lambdam +par.beta.*(1-par.mum).*vt_plus1  ,[],3  );
%        
%     v1 = max(v_t,v_mw2);
%     
%     error = abs(v0-v1);
%     error = max(error(:));
% 
% %     disp(['iteration '  int2str(i)  ': '  num2str(error)]);
%     if error < maxError;
%         numberOfIteration_mw = i;
%         break;
%     end;
%     v0=v1;
% end;
% 
% v_ms2 = v0;
% h_ms2 = h0;
% clear h0 v0 v1 c aMatrix i;
% 
% 
% toc;

%% healthy and retired

v0 = repmat(0,length(par.a),1);
h0 = repmat(0,length(par.a),1);

           
aMatrix=repmat(par.a,1,length(par.a));

c=(1+par.r)./(1-par.muh).*aMatrix -aMatrix';
c(c<0)=nan;

for i=1:maxInteration
    [v1 h0]=max(utility(c)+par.beta.*(1-par.muh-par.theta).*(repmat(v0',length(par.a),1))+par.beta.*par.theta.*(repmat(v_mr',length(par.a),1)) ,[],2);
    error = max(abs(v0-v1));
%     disp(['iteration '  int2str(i)  ': '  num2str(error)]);
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
tic
for j=length(par.edu):-1:1
    v0 = repmat(0,length(par.a),1);
    h0 = repmat(0,length(par.a),1);

    aMatrix=repmat(par.a,1,length(par.a));

    c=(1+par.r)./(1-par.muh).*aMatrix -aMatrix' + wage(par.edu(j),par);
    c(c<0)=nan;

    for i=1:maxInteration
        [v_t h0]=max(utility(c)-par.lambdah +par.beta.*(1-par.muh-par.theta).*(repmat(v0',length(par.a) ,1)) ...
 +   par.beta * par.theta *   (repmat(v_mw(:,j)',length(par.a) ,1))  ,[],2) ;
        v1 = max(v_t,v_hr);
        error = max(abs(v0-v1));
%         disp(['iteration '  int2str(i)  ': '  num2str(error)]);
        if error < maxError;
            numberOfIteration_mw = i;
            break;
        end;
        v0=v1;
    end;

    v_hw(:,j) = v0;
    h_hw(:,j) = h0;
    clear h0 v0 v1 c aMatrix i;
end;
toc

% %% healthy and working stage (method 2)
% tic
% % 1st dimension: asset (state variable)
% % 2nd dimension: education (state variable)
% % 3rd dimension: at+1 (control state variable)
% % so the dimension would be [length(par.a) length(par.edu) length(par.a)]
% v0 = repmat(0,[length(par.a) length(par.edu)]); 
% h0 = repmat(0,[length(par.a) length(par.edu)]);
% 
% % create the 3rd dimension for control state varible. Since it is still
% % par.a, so we use repmat(par.a,...) to createthe matrix
% at_Matrix = permute(repmat(par.a,1,length(par.a)),[1 3 2] );
% at_plus1_matrix = permute(repmat(par.a,1,length(par.a)),[2 3 1] );
% 
% % create an edu with correct dimensionality 
% edu = repmat(par.edu',[ length(par.a) 1 length(par.a) ] );
% 
% 
% c= repmat( (1+par.r)./(1-par.muh).*at_Matrix -at_plus1_matrix, [1 length(par.edu) 1 ]  ) + wage(edu,par); 
% c(c<0)=nan;
% 
% for i=1:maxInteration
%     %
%     vt_plus1 = repmat(permute(v0,[3 2 1]), [length(par.a) 1 1]) ;
%     vt_mw_plus1 = repmat(permute(v_mw2,[3 2 1]), [length(par.a) 1 1]) ;
%     
%     [v_t h0]=max(utility(c)-par.lambdah +par.beta.*(1-par.muh-par.theta).* vt_plus1 + par.beta*par.theta* vt_mw_plus1 ,[],3  ); 
%        
%     v1 = max(v_t,repmat(v_hr,[1 length(par.edu)]));
%     error = abs(v0-v1);
%     error = max(error(:));
% %     disp(['iteration '  int2str(i)  ': '  num2str(error)]);
%     if error < maxError;
%         numberOfIteration_mw = i;
%         break;
%     end;
%     v0=v1;
% end;
% 
% v_hw2 = v0;
% h_hw2 = h0;
% clear h0 v0 v1 c aMatrix i;
% toc







%% healthy and schooling stage (method 1)

% method 1
% loop over different value of education. Easy to code but difficult to
% extend
tic;
for j=length(par.edu):-1:1
    v0 = repmat(0,length(par.a),1);
    h0 = repmat(0,length(par.a),1);

    aMatrix=repmat(par.a,1,length(par.a));

    c=(1+par.r)./(1-par.muh).*aMatrix -aMatrix';
    c(c<0)=nan;

    for i=1:maxInteration
        if par.edu(j)==par.edu(end)
            vt_plus1 = (repmat(v0',length(par.a),1));
            vt_ms_plus1 = (repmat(v_ms(:,j)',length(par.a),1));
        else 
            vt_plus1 = (repmat(v_hs(:,j+1)',length(par.a),1));
            vt_ms_plus1 = (repmat(v_ms(:,j+1)',length(par.a),1));
        end

        [v_t h0]=max(utility(c)-par.lambdah +par.beta.*(1-par.muh-par.theta ).*vt_plus1 + par.beta*par.theta*vt_ms_plus1 ,[],2);
        v1 = max(v_t,v_hw(:,j) );
        error = max(abs(v0-v1));
%         disp(['iteration '  int2str(i)  ': '  num2str(error)]);
        if error < maxError;
            numberOfIteration_mw = i;
            break;
        end;
        v0=v1;
        
    end;

    v_hs(:,j) = v0;

    h_hs(:,j) = h0;

    clear h0 v0 v1 c aMatrix i;
end;
toc

% %% morbid and schooling stage (method 2)
% tic
% v0 = repmat(0,[length(par.a) length(par.edu)]); 
% h0 = repmat(0,[length(par.a) length(par.edu)]);
% 
% at_Matrix = permute(repmat(par.a,1,length(par.a)),[1 3 2] );
% at_plus1_matrix = permute(repmat(par.a,1,length(par.a)),[2 3 1] );
% edu = repmat(par.edu,[ length(par.a) 1 length(par.a) ] );
% 
% c= repmat( (1+par.r)./(1-par.muh).*at_Matrix -at_plus1_matrix, [1 length(par.edu) 1 ]) ; 
% c(c<0)=nan;
% 
% 
% for i=1:maxInteration
%     
%     vt_plus1 = repmat(permute(v0(:,[2:end end]),[3 2 1]), [length(par.a) 1 1]) ;
%     vt_ms_plus1 = repmat(permute(v_ms2(:,[2:end end]),[3 2 1]), [length(par.a) 1 1]) ;
%     
%     [v_t h0]=max(utility(c)-par.lambdah +par.beta.*(1-par.muh).*vt_plus1 + par.beta*par.theta*vt_ms_plus1  ,[],3  );
%        
%     v1 = max(v_t,v_hw2);
%     
%     error = abs(v0-v1);
%     error = max(error(:));
% 
% %     disp(['iteration '  int2str(i)  ': '  num2str(error)]);
%     if error < maxError;
%         numberOfIteration_mw = i;
%         break;
%     end;
%     v0=v1;
% end;
% 
% v_hs2 = v0;
% h_hs2 = h0;
% clear h0 v0 v1 c aMatrix i;
% 
% 
% toc;


% compute the critial value;
[temp critical_point_morbid_schooling]=max(v_ms==v_mw)
[temp critical_point_morbid_retirement] = max(v_mw==repmat(v_mr,[1 4]))
  
[temp critical_point_healthy_schooling]=max(v_hs==v_hw)
[temp critical_point_healthy_retirement] = max(v_hw==repmat(v_hr,[1 4]))
  
critical_point= [critical_point_morbid_schooling;...
                    critical_point_morbid_retirement;...
                    critical_point_healthy_schooling;...
                    critical_point_healthy_retirement]


save valueFn v_hr v_hw v_hs v_mr v_mw v_ms par;



