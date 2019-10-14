function phi = phi(previous_belief,current_action,Power_Obs,angle_obs,rho,power_obs_variance,Angle_obs_variance,action_at_time_minus_1 ) %this is the recursion function

    phi(1:14) = 0;
    
    if (previous_belief(6) == 1)
        u_k_minus_1 = action_at_time_minus_1; % for the first repetition
    else
        if (previous_belief(1) ~= 0) | (previous_belief(2) ~= 0) |(previous_belief(3) ~= 0) |(previous_belief(4) ~= 0) |(previous_belief(5) ~= 0) |(previous_belief(6) ~= 0) |(previous_belief(7) ~= 0) |(previous_belief(8) ~= 0) 
            u_k_minus_1 = 0; % 0 for unshed 1 for non shed
        elseif (previous_belief(9) ~= 0) | (previous_belief(10) ~= 0) |(previous_belief(11) ~= 0) |(previous_belief(12) ~= 0) |(previous_belief(13) ~= 0) |(previous_belief(14) ~= 0) 
            u_k_minus_1 = 1;
        else
            u_k_minus_1 = 2;
        end
    end
    
                for j=1:1:14
                      [Power,Gamma,action] = identify(j);
                      Theta = asin(Power/Gamma);
                      for k=1:1:14 % to calculate the sum
                         [Previous_Power,Previous_Gamma,Previous_Action] = identify(k);
                         Previous_Theta = asin(Previous_Power/Previous_Gamma);
                         phi(j) = phi(j) +previous_belief(k)*transition(Previous_Power,Previous_Gamma,Previous_Action,Power,Gamma,action,current_action-1,u_k_minus_1,rho)*normpdf(Power_Obs,Power,power_obs_variance)*normpdf(angle_obs,Theta,Angle_obs_variance);
                        % transition(Previous_Power,Previous_Gamma,Previous_Action,Power,Gamma,action,current_action-1,u_k_minus_1,rho)
                      end   % -1 is needed because current action is in [1,3] but we want to be at [0,2]
                end
                phi = phi/sum(phi); %normalize
end


