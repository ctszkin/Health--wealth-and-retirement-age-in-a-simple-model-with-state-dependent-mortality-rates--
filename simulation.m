clc;clear;
load valueFn

numberOfAgent = 1000;


health_status = stateGen(par,numberOfAgent);
lifeSpan = sum(health_status,2);

a=nan(max(lifeSpan),numberOfAgent);
c=nan(max(lifeSpan),numberOfAgent);

retirementAge = nan(1,numberOfAgent);
a(1,:)=0;
% c(1,:)=0;

retirement_counter = zeros(numberOfAgent,1);
% 0 if the agent retire before death, otherwise 1

for i=1:numberOfAgent
    % healthy period
    for j=1:(health_status(i,1))
        a(j+1,i) = api_hw(par.a==a(j,i));
        c(j,i) = c_hw(par.a==a(j,i));
        
        % In order to trace the record explicitly, we can use the cmd 'disp' 
        % to ask the matlab to print out some of the information. 
        % Here's an example: 
        % comment them if you dont need it
        % copy the code below to morbid part if necessary

%         disp(['Agent: ' num2str(i) '  Period: ' num2str(j)]);
%         disp(['state variable: a(' num2str(j) ')=' num2str(a(j,i)) ]);
%         disp(['control variable: a(' num2str(j+1) ')=' num2str(a(j+1,i)) '  c(' num2str(j,i) ')=' num2str(c(j,i))   ]);
%         disp(' ');
        
        % record the time for retirement if a >= critial pt of a 
        % once retired, set retirement_counter(i)=1
        if ( a(j,i)>=par.a(critial_pt_h) && retirement_counter(i)==0 )
            retirementAge(i)=j;
            retirement_counter(i)=1;
        end
    end

    % morbid period
    for j=(health_status(i,1)+1):(lifeSpan(i))
        a(j+1,i) = api_mw(par.a==a(j,i));
        c(j,i) = c_mw(par.a==a(j,i));
        if ( a(j,i)>=par.a(critial_pt_m) && retirement_counter(i)==0 )
            retirementAge(i)=j;
            retirement_counter(i)=1;
        end
    end
    c(lifeSpan(i),i)=c_mw(par.a==a(lifeSpan(i),i)) ;
end


% simulationresult
% Another additional elements inlcude, lifeSpan

% Basic elements: c a retirementAge healthyTime morbidTime
% c and a is a period x agent matrix
% retirementAge healthyTime morbidTime is a vector with length equal to the
% number of agent.
simResult.par=par;
simResult.numberOfAgent=numberOfAgent;

simResult.healthyTime = health_status(:,1);
simResult.morbidTime = health_status(:,2);
simResult.c = c;
simResult.a = a;
simResult.retirementAge=retirementAge;
simResult.trimmedRetirementAge=simResult.retirementAge(~isnan(simResult.retirementAge));
simResult.numberOfAgentDieB4Retirement = simResult.numberOfAgent - length(simResult.trimmedRetirementAge);
simResult.percentageOfAgentDieB4Retirement = simResult.numberOfAgentDieB4Retirement/simResult.numberOfAgent;

% Additional elements:
simResult.lifeSpan = sum(health_status,2);

simResult.avg_c=nan(max(lifeSpan),1);
simResult.avg_a=nan(max(lifeSpan),1);
for i=1:max(lifeSpan)
    simResult.avg_c(i) = mean(c(i,~isnan(c(i,:))));
    simResult.avg_a(i) = mean(a(i,~isnan(a(i,:))));
end

% in terms of percentage
simResult.pop=nan(max(lifeSpan),1);
for i=1:max(lifeSpan)
    simResult.pop(i) = sum(lifeSpan>=i)./numberOfAgent;
end

% Summary

% average healthy period and average morbid period
disp(['Average healthy period : ' num2str(mean(simResult.healthyTime)) ])
disp(['Average morbid period : ' num2str(mean(simResult.morbidTime)) ])

% average life
disp(['Average ife-span : ' num2str(mean(simResult.lifeSpan)) ])

% average retirement age
disp(['Average retirement age (Drop out those who die b4 retirement): ' num2str(mean(simResult.trimmedRetirementAge ))]);
disp(['Percentage of agent who die b4 retirement: ' num2str(simResult.percentageOfAgentDieB4Retirement)]);



% avg consumption & avg wealth & population  at period t;
t = 50;
disp(['Average consumption at time ' int2str(t) ' : ' num2str(simResult.avg_c(t))]);
disp(['Average wealth at time ' int2str(t) ' : ' num2str(simResult.avg_a(t))]);
disp(['Population at time ' int2str(t) ' : ' num2str(simResult.pop(t))]);





