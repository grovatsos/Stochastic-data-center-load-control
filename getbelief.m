function [belief,action_at_time_zero]=getbelief(rho)

breakdown=0;
Power_obs_variance = 0.5; %variance of the power observations, this may change
Angle_obs_variance = 0.5; %variance of the angle observations, this may change

shutdown_rho=0.005; %not 0.5 since shutdown is rare

 %parameter for gamma going from 1 to 0.5



belief(1:14) = 0;
belief(6)=1;
max_sample_number = geornd(shutdown_rho)+3; %number of belief samples before shutdown // ADD 2 BECAUSE I PRACTICALY DONT SHUTDOWN AT TIME INSTANT ONE
%may be smaller due to breakdowns


%
%========================generate the belief samples  / sample 1
%is already initialized================
%above=======================

%generate the gammas 1---> 0.5 can happen. Otherway around can not



%changepoint =floor(max_sample_number/3) ; 
changepoint =geornd(rho)+3;
%get the initial power and initial gamma from the initial belief%%%

Power(1:1:max_sample_number)=0;

r= randi([1 14],1,1);
r=9;
if r ==1 
    Power(1) = -0.4;
    action_at_time_zero = 0;
    initial_gamma = 0.5;    
elseif r==2
    Power(1) = 0;
    action_at_time_zero = 0;
    initial_gamma = 0.5; 
elseif r ==3
    Power(1) = 0.4;
    action_at_time_zero = 0;
    initial_gamma = 0.5; 
elseif r==4
    Power(1) = -0.4;
    action_at_time_zero = 1;
    initial_gamma = 0.5; 
elseif r ==5 
    Power(1) =  0;
    action_at_time_zero = 1;
    initial_gamma = 0.5; 
elseif r==6
    Power(1) = 0.4;
    action_at_time_zero = 1;
    initial_gamma = 0.5; 
elseif r ==7
    Power(1) = -0.8;
    action_at_time_zero = 0;
    initial_gamma = 1; 
elseif r==8
    Power(1) =  -0.4;
    action_at_time_zero = 0;
    initial_gamma = 1; 
elseif r ==9 
    Power(1) =  0;
    action_at_time_zero = 0;
    initial_gamma = 1; 
elseif r==10
    Power(1) = 0.4;
    action_at_time_zero = 0;
    initial_gamma = 1; 
elseif r ==11 
    Power(1) = 0.8;
    action_at_time_zero = 0;
    initial_gamma = 1; 
elseif r==12
    Power(1) = -0.4;
    action_at_time_zero = 1;
    initial_gamma = 1; 
elseif r ==13 
    Power(1) = 0;
    action_at_time_zero = 1;
    initial_gamma = 1; 
else 
    Power(1) = 0.4;
    action_at_time_zero = 1;
    initial_gamma = 1; 
end



%generate random_actions/the last one will be the shutdown
%0 = non_shed///1=shed//2 shutdown
for i=1:1:max_sample_number
    random_actions(i)=binornd(1,0.5);
end
% random_actions(max_sample_number-499:max_sample_number-1)=1; 

random_actions(max_sample_number)=2;
% random_actions=ones(1,max_sample_number);
%specify the gamma values at each step
if initial_gamma == 0.5
    gamma(1:max_sample_number)=0.5;
else
    if changepoint > max_sample_number
        gamma(1:max_sample_number) = 1;
    else
        gamma(1:changepoint - 1) = 1 ;
        gamma(changepoint:max_sample_number) = 0.5 ;
    end
end



%We assume only continue action / not shutdown / also we exclue shutdowns
%Generate the power for the other points (outside of initial)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%need to be careful with shed load action
% if P=0.8, and shed load, it is reduced to 0.4, and then changes with the
% given statistics


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=2:1:max_sample_number
    if random_actions(i-1) == 0 % if unshed at previous time instant //
            if ( Power(i-1) ==  -0.8 )
                 if (binornd(1,0.5) == 0) %power stays the same
                     if gamma(i)==1
                        Power(i) = Power(i-1);
                     else
                        Power(i)=-0.4;
                     end
                 else
                    Power(i) = -0.4;
                 end   
            elseif ( Power(i-1) ==  0.8 )
                 if (binornd(1,0.5) == 0) %power stays the same
                    if gamma(i)==1
                      Power(i) = Power(i-1);
                    else
                        Power(i)=0.4;
                    end
                 else
                    Power(i) = 0.4;
                 end
            elseif gamma(i)==1
                x = mnrnd(1,[1/3 1/3 1/3]);
                if x(1) == 1
                    Power(i) = Power(i-1) + 0.4;
                end
                if x(2) == 1
                    Power(i) = Power(i-1);
                end
                if x(3) == 1
                    Power(i) = Power(i-1) -0.4;
                end
            else
                if binornd(1,0.5)==1
                    Power(i)=Power(i-1);
                elseif Power(i-1)~=0
                    Power(i)=0;
                else
                    Power(i)=0.4*(binornd(1,0.5)*2-1);
                end
                    
            end
    elseif random_actions(i-1) == 1 %shed
        if ( Power(i-1) ==  -0.8 )
            Power(i) = -0.4;
        elseif ( Power(i-1) ==  0.8 )
            Power(i) = 0.4;
        elseif ( Power(i-1) ==  -0.4 )
             if (binornd(1,0.5) == 0) %power stays the same
                Power(i) = Power(i-1);
             else
                Power(i) = 0;
             end 
         elseif ( Power(i-1) ==  0.4 )
             
             if (binornd(1,0.5) == 0) %power stays the same
                Power(i) = Power(i-1);
             else
                Power(i) = 0;
             end 
        else %Power(i-1) =0
            
            x = mnrnd(1,[1/3 1/3 1/3]);
            if x(1) == 1
                Power(i) = Power(i-1)+ 0.4;
            end
            if x(2) == 1
                Power(i) = Power(i-1);
            end
            if x(3) == 1
                Power(i) = Power(i-1) -0.4;
            end
            
        end
        
        
    else%shut down // probably not needed
        Power(i)=0;
    end
    if ((gamma(i) == 0.5) &( abs(Power(i)) == 0.8  )) %find when the shutdown happens    
        breakdown = 1
        break
    end
end

%truncate the Power, gamma and actions if breakdown occurs as to simulate a
%real system
sample_number = i ;%final sample number
if breakdown == 1 %if that holds breakdown occured at time instant i
    gamma=gamma(1:i);
    random_actions=random_actions(1:i);
    Power=Power(1:i);
    
end

%generate the observations Z_k which are noisy version of theta_k
%and P_k
Obs_z(1:sample_number)=0;
Obs_theta(1:sample_number)=0;

for i=1:1:sample_number
    Obs_z(i) = normrnd(Power(i),Power_obs_variance);
    Theta(i) = asin(Power(i) / gamma(i));
    Obs_theta(i) = normrnd( Theta(i) , Angle_obs_variance); % do we need to add i.i.d. gaussian noise in each component
end

if breakdown ==1 %If i finished by breakdown the last observations are not acquired
   Theta(length(Theta))=[];
    Power(length(Power))=[];
    Obs_z(length(Obs_z))=[];
    Obs_theta(length(Obs_theta))=[];
    random_actions(length(random_actions))=[];
    sample_number = sample_number - 1 ;
    gamma(length(gamma)) = [];
end
    
    

%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%

% use the recursion to generate the belief points
%expand the belief vector over time
belief_initial=belief;
clear belief %ALWAYS CLEAR VARIABLES LIKE THIS WHEN YOU GIVE THEM VALUE AGAIN
belief(1:14,1:sample_number)=0;
belief(:,1)=belief_initial ;



for i=2:1:sample_number
    for j=1:1:14
          for k=1:1:14 % to calculate the sum
              if i  == 2 
                 [Power_now,Gamma_now,Action_now] = identify(j); %transition from k to j
                 Theta_now = asin(Power_now/Gamma_now);
                 [Previous_Power,Previous_Gamma,Previous_Action] = identify(k);
                 belief(j,i) = belief(j,i) +belief(k,i-1)*transition(Previous_Power,Previous_Gamma,Previous_Action,Power_now,Gamma_now,Action_now,random_actions(i-1),action_at_time_zero,rho)*normpdf(Obs_z(i),Power_now,Power_obs_variance)*normpdf(Obs_theta(i),Theta_now,Angle_obs_variance);
              else
                 [Power_now,Gamma_now,Action_now] = identify(j); %transition from k to j
                 Theta_now = asin(Power_now/Gamma_now);
                 [Previous_Power,Previous_Gamma,Previous_Action] = identify(k);
                 belief(j,i) = belief(j,i) +belief(k,i-1)*transition(Previous_Power,Previous_Gamma,Previous_Action,Power_now,Gamma_now,Action_now,random_actions(i-1),random_actions(i-2),rho)*normpdf(Obs_z(i),Power_now,Power_obs_variance)*normpdf(Obs_theta(i),Theta_now,Angle_obs_variance);
              end
          end   
    end
    belief(:,i)=belief(:,i)/(sum(belief(:,i)));      
end
belief;
belief(:,1)=[];
end
