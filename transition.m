function transition = transition(Previous_Power,Previous_Gamma,Previous_Action,Power,Gamma,Action,current_action,u_k_minus_1,rho)
%transition probability from Previous_Power,Previous_Gamma,Previous_Action to Power,Gamma,Action

%calculateting prob(gamma_k+1 | gamma_k)
if ((Previous_Gamma==0.5 )&(Gamma==0.5)) 
    gamma_trans = 1;
elseif ((Previous_Gamma==0.5)&(Gamma==1))
    gamma_trans = 0;
elseif ((Previous_Gamma==1)&(Gamma==1))
    gamma_trans = 1 - rho;
else
    gamma_trans = rho ;
end

%calculating prob(P_{k+1)|P_k,u_k) \\\\ n,k = 1 or 2 or 3 depending on the 
%action u_k = 
%b,n,m ---> j,k,l
if abs(Previous_Power-Power)>0.4
    power_trans = 0;
else
    if Action == 0 % if u_k = non-shed
        
        if Previous_Power == Power %stay same power
            if ((Previous_Power==-0.8)|(Previous_Power==0.8))
                power_trans = 1/2;
            else
                power_trans = 1/3;
            end
        else % change power
%             power_trans = 1/3;%wrong
            if abs(Previous_Power)==0.8
                power_trans=1/2;
            else
                power_trans=1/3;
            end
        end
        
    elseif Action==1 % if u_k = shed
        if ((Power==-0.8)|(Power==0.8))
            power_trans = 0;
        else
            if Previous_Power==Power
                if Previous_Power==0
                    power_trans = 1/3;
                else
                    power_trans =1/2;
                end
            else
                if ((Previous_Power==-0.8)|(Previous_Power==0.8))
                    power_trans = 1;
                elseif ((Previous_Power==-0.4)|(Previous_Power==0.4))
                    power_trans = 1/2;
                else
                    power_trans = 1/3;
                end
            
            end
        end
    else
        1
        %if i shutdown keep only power = 0 ?? Probably not needed
        if Power==0
            power_trans = 1;
        else
            power_trans = 0;
        end
    end
end


%calculate Prob(uk|uk-1) // compare the actions with certain random actions
if ((Previous_Action == u_k_minus_1) &(Action==current_action)) %since in random actions 0 = unshed 1 =shed
    action_trans =1;
else 
    action_trans=0;
end

transition = action_trans*power_trans*gamma_trans;

    
    

    

end