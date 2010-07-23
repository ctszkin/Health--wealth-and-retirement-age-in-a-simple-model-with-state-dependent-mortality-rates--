function result = agentDecision(state)      
        lifeSpan=sum(state);
%       result: age health status e a c
        result=nan(lifeSpan,6)
        result(:,1)=1:lifeSpan
        result(1: state(1),2)=0
        result(state(1)+1: state(2)+state(1),2)=1


        
        [temp pos_of_a_zero]=max(par.a==0)
        
        result(1,5)=pos_of_a_zero;
        result(1,4)=0;
        
        % period 1
        if result(1,2)==0  % if healthy
            if result(1,5) > critical_point_healthy_schooling(1)
                result(2,5) = h_hs(result(1,5),1);
                result(1,3) = 0;
                result(2,4) = 1;
            elseif result(1,5) >= critical_point_healthy_retirement(1)
                result(2,5) = h_hr(result(1,5));
                result(1,3) = 2;
                result(2,4) = 0;
            else
                result(2,5) = h_hw(result(1,5),1);
                result(1,3) = 1;
                result(2,4) = 0;
            end 
        end
        if result(1,2)==1  % if morbid
            if result(1,5) > critical_point_morbid_schooling(1)
                result(2,5) = h_ms(result(1,5),1);
                result(1,3) = 0;
                result(2,4) = 1;
            elseif result(1,5) >= critical_point_morbid_retirement(1)
                result(2,5) = h_mr(result(1,5));
                result(1,3) = 2;
                result(2,4) = 0;
            else
                result(2,5) = h_mw(result(1,5),1);
                result(1,3) = 1;
                result(2,4) = 0;
            end 
        end
        
        
        % period 1+
        for j=2:size(result,1)            
            if result(j,2)==0  % if healthy
                if result(j,5) > critical_point_healthy_schooling(result(j,4)+1)
                    result(j+1,5) = h_hs(result(j,5),result(j,4)+1);
                    result(j,3) = 0;
                    if result(j,4)+1 < par.max_edu
                        result(j+1,4) = result(j,4)+1;
                    else
                         result(j+1,4) = result(j,4);
                    end
                elseif result(j,5) >= critical_point_healthy_retirement(result(j,4)+1)
                    result(j+1,5) = h_hr(result(j,5));
                    result(j,3) = 2;
                    result(j+1,4) = result(j,4);
                else
                    result(j+1,5) = h_hw(result(j,5),result(j,4)+1);
                    result(j,3) = 1;
                    result(j+1,4) = result(j,4);
                end 
            end
            if result(j,2)==1  % if morbid
                if result(j,5) > critical_point_morbid_schooling(result(j,4)+1)
                    result(j+1,5) = h_ms(result(j,5),result(j,4)+1);
                    result(j,3) = 0;
                    if result(j,4)+1 < par.max_edu
                        result(j+1,4) = result(j,4)+1;
                    else
                         result(j+1,4) = result(j,4);
                    end
                elseif result(1,5) >= critical_point_morbid_retirement(result(j,4)+1)
                    result(j+1,5) = h_mr(result(j,5));
                    result(j,3) = 2;
                    result(j+1,4) = result(j,4);
                else
                    result(j+1,5) = h_mw(result(j,5),result(j,4)+1);
                    result(j,3) = 1;
                    result(j+1,4) = result(j,4);
                end 
            end
        end
        
        
    end