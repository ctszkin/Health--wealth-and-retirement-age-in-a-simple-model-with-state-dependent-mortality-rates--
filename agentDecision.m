  function result = agentDecision(state,model)
        
        HEALTHY=0;
        MORBID=1;
        DEAD=2;
        STUDYING=0;
        WORKING=1;
        RETIRED=2;
        


%         HEALTHY='Healthy';
%         MORBID='Morbid';
%         DEAD='Dead';
%         STUDYING='Studying';
%         WORKING='Working';
%         RETIRED='Retired';
%         



        lifeSpan=sum(state);

        result.age = (1:(lifeSpan+1))';
        
        result.health = nan(lifeSpan+1,1);
        result.health(1: state(1))=HEALTHY;
        result.health(state(1)+1: state(2)+state(1))=MORBID;
        result.health(lifeSpan+1)=DEAD;
        
       result.employmentStatus= nan(lifeSpan+1,1);
       result.yearOfEducation= nan(lifeSpan+1,1);
       result.asset = nan(lifeSpan+1,1);
       result.status= nan(lifeSpan+1,1);
       
        
        [temp postitionOfZeroInAsset]=max(model.par.a==0);

        
        % period 1
        result.asset(1)=postitionOfZeroInAsset;
        result.yearOfEducation(1)=0;

        if result.health(1)==HEALTHY  
            if result.asset(1) <  model.criticalPoint.healthy_schooling2working(1)
                result.status(1) = STUDYING;
                result.asset(2) = model.policyFn.healthy_schooling(result.asset(1),1);
                result.yearOfEducation(2) = 1;
            elseif result.asset(1) >= model.criticalPoint.healthy_working2retired(1)
                result.status(1) = RETIRED;
                result.asset(2) =  model.policyFn.healthy_retired(result.asset(1));
                result.yearOfEducation(2) = 0;
            else
                result.status(1) = WORKING;
                result.asset(2) = model.policyFn.healthy_working(result.asset(1),1);
                result.yearOfEducation(2) = 0;
            end 
        end
        if result.health(1)==MORBID 
            if result.asset(1) < model.criticalPoint.morbid_schooling2working(1)
                result.status(1) = STUDYING;
                result.asset(2) = model.policyFn.morbid_schooling(result.asset(1),1);
                result.yearOfEducation(2) = 1;
             elseif result.asset(1) >= model.criticalPoint.morbid_working2retired(1)
                result.status(1) = RETIRED;
                 result.asset(2) =  model.policyFn.morbid_retired(result.asset(1));
                result.yearOfEducation(2) = 0;
            else
                result.status(1) = WORKING;
                result.asset(2) = model.policyFn.morbid_working(result.asset(1),1);
                result.yearOfEducation(2) = 0;
            end 
        end
        
        
        % period 1+
        for j=2:lifeSpan            
            if result.health(j)==HEALTHY  
                if result.asset(j) < model.criticalPoint.healthy_schooling2working(result.yearOfEducation(j)+1)
                    result.status(j) = STUDYING;
                    result.asset(j+1) = model.policyFn.healthy_schooling(result.asset(j),result.yearOfEducation(j)+1);
                    result.yearOfEducation(j+1) = result.yearOfEducation(j)+1;
                elseif result.asset(j) >= model.criticalPoint.healthy_working2retired(result.yearOfEducation(j)+1) 
                    result.status(j) = RETIRED;
                    result.asset(j+1) = model.policyFn.healthy_retired(result.asset(j));
                    result.yearOfEducation(j+1) = result.yearOfEducation(j);
                else
                    result.status(j) = WORKING;
                    result.asset(j+1) = model.policyFn.healthy_working(result.asset(j),result.yearOfEducation(j)+1);
                    result.yearOfEducation(j+1) = result.yearOfEducation(j);
                end 
            end
            if result.health(j)==MORBID 
                if result.asset(j) < model.criticalPoint.morbid_schooling2working(result.yearOfEducation(j)+1)
                    result.status(j) = STUDYING;
                    result.asset(j+1) = model.policyFn.morbid_schooling(result.asset(j),result.yearOfEducation(j)+1);
                    result.yearOfEducation(j+1) = result.yearOfEducation(j)+1;
                elseif result.asset(j) >= model.criticalPoint.morbid_working2retired(result.yearOfEducation(j)+1)
                    result.status(j) = RETIRED;
                    result.asset(j+1) = model.policyFn.morbid_retired(result.asset(j));
                    result.yearOfEducation(j+1) = result.yearOfEducation(j);
                else
                    result.status(j) = WORKING;
                    result.asset(j+1) = model.policyFn.morbid_working(result.asset(j));
                    result.yearOfEducation(j+1) = result.yearOfEducation(j);
                end 
            end
        end
        result.asset =model.par.a(result.asset)
 end
    


