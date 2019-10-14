function psup = psup(random_belief,current_action,Obs_z,Obs_theta,Power_obs_variance,Angle_obs_variance,rho,action_at_time_minus_1 ) %distribution of next data point given I know the present belief and present action
    
  %SOSOSOS  ==NaN does not work!!!! Need to use isnan() command

 %find the previous action
  if (random_belief(1) == 1/14) & (random_belief(2) == 1/14) &(random_belief(3) == 1/14) &(random_belief(4) == 1/14) &(random_belief(5) == 1/14) &(random_belief(6) == 1/14) &(random_belief(8) == 1/14) &(random_belief(8) == 1/14) &(random_belief(9) == 1/14) &(random_belief(10) == 1/14) &(random_belief(11) == 1/14) &(random_belief(12) == 1/14) &(random_belief(13) == 1/14) &(random_belief(14) == 1/14) 
        u_k_minus_1 =action_at_time_minus_1 ; % for the first repetition
  else
        if (random_belief(1) ~= 0) | (random_belief(2) ~= 0) |(random_belief(3) ~= 0) |(random_belief(4) ~= 0) |(random_belief(5) ~= 0) |(random_belief(6) ~= 0) |(random_belief(7) ~= 0) |(random_belief(8) ~= 0) 
            u_k_minus_1 = 0; % 0 for unshed 1 for non shed
        elseif (random_belief(9) ~= 0) | (random_belief(10) ~= 0) |(random_belief(11) ~= 0) |(random_belief(12) ~= 0) |(random_belief(13) ~= 0) |(random_belief(14) ~= 0) 
            u_k_minus_1 = 1;
        else
            u_k_minus_1 = 2;
        end
  end

%for debugging
%     action=u;
%     Obs_z = v;
%     Obs_theta = q;
    
    
            future_belief(1:14)=0;
            for i = 1:1:14
                [Power,Gamma,action] = identify(i);
                Theta = asin(Power/Gamma);
                for k = 1:1:14
                    [Previous_Power,Previous_Gamma,Previous_Action] = identify(k);
                    Previous_Theta = asin(Previous_Power/Previous_Gamma);
                    future_belief(i) = future_belief(i) + (random_belief(k))*(transition(Previous_Power,Previous_Gamma,Previous_Action,Power,Gamma,action,current_action-1,u_k_minus_1,rho));
                    %(transition(Previous_Power,Previous_Gamma,Previous_Action,Power,Gamma,action,current_action-1,u_k_minus_1,rho))
                end
            end
            
            breakdown_prob = 1-sum(future_belief); % the remaining probability is the probability of breaking down
            
            if ((isnan(Obs_z)) & (isnan(Obs_theta) )) %the probability of getting NaN on next is the probability of breaking down on next step
                psup = breakdown_prob;
            else %the probability of receiving regular observations (Non -NaN)
            
                z_prob = 0;
                theta_prob =0;
                %scale according to probabilities of the pairs
                    for i =1:1:14
                        [Power,Gamma] = identify(i); %identify the power and gamma for this i
                        Angle= asin(Power/Gamma);
                        z_prob = z_prob + normpdf(Obs_z,Power,Power_obs_variance)*future_belief(i);
                        %for this theta make sure the first norm in the
                        %theta_prob does not produce a complex number
                        theta_prob = theta_prob +  normpdf(Obs_theta,Angle,Angle_obs_variance)*future_belief(i) ; %average over power and gamma choices
                    end
                psup=theta_prob*z_prob;
            end

end