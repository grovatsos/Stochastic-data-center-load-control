function [ hyper_f_shed,hyper_f_unshed ] = pesues( hyper_shed,hyper_unshed,belief_P )


 
% hyper_planes_non_shed = rand(3,14)*10;
 hyper_planes_shed =hyper_shed;
 
 
 %rand(3,14)*10;
 hyper_planes_non_shed =hyper_unshed;

 breakdown=0;
Power_obs_variance = 1; %variance of the power observations, this may change
Angle_obs_variance = 1; %variance of the angle observations, this may change

shutdown_rho=0.001; %not 0.5 since shutdown is rare
rho = 0.01;  %parameter for gamma going from 1 to 0.5

Cb=100; %breakdown cost
Cs_0=3; %shut down cost
Cs_1=20; %shut down cost
Cd_0=0.2;
Cd_1=0.2;

belief(1:14) = 0;
belief(6)=1;



%%%%%%%%%Check the zeros in the belief vectors

%=============== UP TO HERE I HAVE ALL THE BELIEF SAMPLES================


%initialize the g vectors // g_0 for nonshed , g_1 for shed , g_2 for
%shutdown
g_0(1:14)=0;
g_1(1:14)=0;
g_2(1:14)=0;

%0 is unshed
g_0(4) = Cb*(rho/2);
g_0(8) = Cb*(rho/2);
g_0(5) = Cb*(rho/3);
g_0(7) = Cb*(rho/3);
g_0(1) = Cb*(1/3);
g_0(3) = Cb*(1/3);
g_0(12) = Cb*(rho/3);
g_0(14) = Cb*(rho/3);
g_0(9) = Cb*(1/3);
g_0(11) = Cb*(1/3);
% 1 is shed
g_1(1)=Cd_1;
g_1(2)=Cd_1;
g_1(3)=Cd_1;
g_1(9)=Cd_1;
g_1(10)=Cd_1;
g_1(11)=Cd_1;
g_1(4)=(Cd_0*(1-rho))  + ( Cd_1*rho);
g_1(5)=(Cd_0*(1-rho))  + ( Cd_1*rho);
g_1(6)=(Cd_0*(1-rho))  + ( Cd_1*rho);
g_1(7)=(Cd_0*(1-rho))  + ( Cd_1*rho);
g_1(8)=(Cd_0*(1-rho))  + ( Cd_1*rho);
g_1(12)=(Cd_0*(1-rho))  + ( Cd_1*rho);
g_1(13)=(Cd_0*(1-rho))  + ( Cd_1*rho);
g_1(14)=(Cd_0*(1-rho))  + ( Cd_1*rho);
%2 is shut down
g_2(9)=Cs_0;
g_2(10)=Cs_0;
g_2(11)=Cs_0;
g_2(12)=Cs_1;
g_2(13)=Cs_1;
g_2(14)=Cs_1;
g_2(1)=(Cd_0+Cd_1+Cb+Cs_0+Cs_1);
g_2(2)=(Cd_0+Cd_1+Cb+Cs_0+Cs_1);
g_2(3)=(Cd_0+Cd_1+Cb+Cs_0+Cs_1);
g_2(4)=(Cd_0+Cd_1+Cb+Cs_0+Cs_1);
g_2(5)=(Cd_0+Cd_1+Cb+Cs_0+Cs_1);
g_2(6)=(Cd_0+Cd_1+Cb+Cs_0+Cs_1);
g_2(7)=(Cd_0+Cd_1+Cb+Cs_0+Cs_1);
g_2(8)=(Cd_0+Cd_1+Cb+Cs_0+Cs_1);

%%%% MAKE SURE THE CODE FROM HERE AND BELOW IS CORRECT %%%%%%%%%%%%%%%%%

 % one random hyperplane for EACH action for a total of two random
 % hyperplanes initially
 % this will probably need to be changed if I see this thing as a function
 
 hyper_planes_non_shed_future = 0;
 hyper_planes_shed_future = 0;
 hyper_planes_future = 0;
 

 hyper_planes_shutdown = [Cd_0+Cd_1+Cb+Cs_0+Cs_1 Cd_0+Cd_1+Cb+Cs_0+Cs_1 Cd_0+Cd_1+Cb+Cs_0+Cs_1 Cd_0+Cd_1+Cb+Cs_0+Cs_1 Cd_0+Cd_1+Cb+Cs_0+Cs_1 Cd_0+Cd_1+Cb+Cs_0+Cs_1 Cd_0+Cd_1+Cb+Cs_0+Cs_1 Cd_0+Cd_1+Cb+Cs_0+Cs_1 Cs_0 Cs_0 Cs_0 Cs_1 Cs_1 Cs_1]; % each row is a hyperplane. Start with a single random hyperplane for SHUTDOWN action

 
 if hyper_planes_non_shed ==0
    hyper_planes = vertcat(hyper_planes_shed,hyper_planes_shutdown);%all hyperplanes independent of action
 elseif hyper_planes_shed == 0
     hyper_planes = vertcat(hyper_planes_non_shed,hyper_planes_shutdown);
 else
     hyper_planes = vertcat(hyper_planes_non_shed,hyper_planes_shed,hyper_planes_shutdown);
 end
 
 random_matrix = [2:1:sample_number];
 h=0;
while 1 %number of repetitions (1000 just a large number // not necessarily that many repetitions are achieved)
     %this loop runs until a stopping conditions is satisfied

    if length(random_matrix) == 0
        break
    end
    
    random_selection = datasample(random_matrix,1);
    random_matrix = random_matrix(find(random_matrix ~= random_selection));
    
    random_belief = belief_samples(:,random_selection); % sample a random belief out of all sampled beliefs
    
    size_hyper_planes_non_shed = size(hyper_planes_non_shed) ;
    size_hyper_planes_shed = size(hyper_planes_shed) ;
    size_hyper_planes_shutdown = size(hyper_planes_shutdown) ;
    
    %CHANGE THEM LATER
    
    if (length(hyper_planes_non_shed) ~= 1 ) & (length(hyper_planes_shed) ~= 1)
         no_planes= size_hyper_planes_non_shed(1) + size_hyper_planes_shed(1) + size_hyper_planes_shutdown(1) ; % determine the total number of current hyper_planes(all actions summed)
    elseif (length(hyper_planes_non_shed) == 1 ) & (length(hyper_planes_shed) ~= 1)
         no_planes= size_hyper_planes_shed(1) + size_hyper_planes_shutdown(1) ; % determine the total number of current hyper_planes(all actions summed)
    elseif (length(hyper_planes_non_shed) ~= 1 ) & (length(hyper_planes_shed) == 1)
        no_planes= size_hyper_planes_non_shed(1) + size_hyper_planes_shutdown(1) ; % determine the total number of current hyper_planes(all actions summed)
    else
        no_planes= size_hyper_planes_shutdown(1) ; % determine the total number of current hyper_planes(all actions summed)
    end

    %====================================================================
    right_term(1:2,1:14)=0; %start calculating a_u^p where the first coordinate is for the action
    for u = 1:1:2 % go through all actions // meaning calculate a_u^p for every u // here we have on action continue
      % u=1 is for nonshed /// u=2 for shed // shutdown hyperplane if fixed 
        %calculate the integral by treating it like a sum
        step_size=0.2;
        Power_obs_range = [-1.6:step_size:1.6];%(-1,1) since my power is mostly around these values
        Angle_obs_range = [-1.6:step_size:1.6];%q = -1:0.001:1 %discritize the angle part of the integral
        for v= Power_obs_range
            for q= Angle_obs_range 

            minimizants(1:no_planes) = 0;
            for i = 1:1:no_planes % calculate the minimum part

                minimizants(i) = dot(phi(random_belief,u,v,q,rho,Power_obs_variance,Angle_obs_variance,action_at_time_zero) , hyper_planes(i,:)); 

            end
            
            [t,argmin] = min(minimizants); %find the arg min       %%=========SOS CHECK PSUP=====
            right_term(u,:) = right_term(u,:) + (step_size^2)*psup(random_belief,u,v,q,Power_obs_variance,Angle_obs_variance,rho,action_at_time_zero )*hyper_planes(argmin,:);
            %right_term(u,:) = right_term(u,:) + minimizants(argmin); 

            end
        end
        toc
        v=NaN;
        q=NaN;
        right_term(u,:)  = right_term(u,:)  +  (step_size^2)*psup(random_belief,u,v,q,Power_obs_variance,Angle_obs_variance,rho,action_at_time_zero )*Cb ; % SInce the min is equal to Cb in this case
        
        if u ==1 
            right_term(u,:) = right_term(u,:) +g_0;
        else
            right_term(u,:) = right_term(u,:) +g_1;
        end
    end
    %==========================================================================================
    
    %find the new candidate hyperplane to add
    %the candidates are between the non-shed and the shed since the
    %shutdown hyperplane is fixed

    hyper_1 = dot(random_belief,right_term(1,:)); 
    hyper_2 = dot(random_belief,right_term(2,:));


     if  (hyper_1 <= hyper_2)   %we have a candidate hyper for non-shed action
         canditate_hyper = right_term(1,:);
         add_to = 1; %meaning it should be added to the hyperplanes attached to non shed
     else
         canditate_hyper = right_term(2,:);
         add_to=2; %meaning it should be added to the hyperplanes attached to shed
     end

    %decide whether to add the new hyperplane

    %calculate J_k(p) which is caclulated as the min of inner product of p
    %to every possible accumulated hyperplane except from the last one

    J_k= dot(random_belief,hyper_planes(1,:));
    if no_planes >1
        for i=2:1:no_planes
            if dot(random_belief ,hyper_planes(i,:)) <= J_k
                J_k = dot(random_belief,hyper_planes(i,:));
            end

        end
    end


    %check the add condition
    
    size_hyper_planes_non_shed_future = size( hyper_planes_non_shed_future) ;
    size_hyper_planes_shed_future = size(hyper_planes_shed_future) ;
    size_hyper_planes_shutdown = size(hyper_planes_shutdown) ;
    size_hyper_planes_future = size(hyper_planes_future);
     
    if dot(random_belief,canditate_hyper) <= J_k
        if add_to == 1 %Add new hyperplane to the non-shed hyperplanes
            if hyper_planes_non_shed_future ==0 %check if this hyperplane class is empty
               hyper_planes_non_shed_future = canditate_hyper;
            else
               if BelongsTo(hyper_planes_non_shed_future,canditate_hyper) == 0 %candidate hyper does not exist in future hyperplanes
                  hyper_planes_non_shed_future(size_hyper_planes_non_shed_future(1) +1 , :) = canditate_hyper;
               end
            end
        else %add to shed hyper_planes
            if hyper_planes_shed_future ==0 %check if this hyperplane class is empty
               hyper_planes_shed_future = canditate_hyper;
            else
               if BelongsTo(hyper_planes_shed_future,canditate_hyper) == 0 %candidate hyper does not exist in future hyperplanes
                  hyper_planes_shed_future(size_hyper_planes_shed_future(1) +1 , :) = canditate_hyper;
               end
            end   
        end
    else 
        %else add the hyperplane that maximizes dot(random_belief,hyperplane)
        for i =1:1:no_planes
            maxz(i) = dot(random_belief,hyper_planes(i,:));
            
        end    
        [w,argm] = min(maxz); 
        %find to which action this hyperplane corresponds
        if argm~= no_planes  %i dont want to readd the shutdown hyperplane  
            if hyper_planes_non_shed == 0 %disregard non-shed
                old_add =2;
                
            elseif hyper_planes_shed == 0 
                
                old_add =1;
            elseif (length(hyper_planes_non_shed) ~= 1) & (length(hyper_planes_shed) ~= 1)
                size_hyper_planes_shed = size(hyper_planes_shed) ;
                size_hyper_planes_non_shed = size(hyper_planes_non_shed) ;
                size_hyper_planes_shutdown = size(hyper_planes_shutdown) ;
                value_shed(1:size_hyper_planes_shed(1)) = 0;
                value_non_shed(1:size_hyper_planes_non_shed(1)) = 0;
                value_shutdown(1:size_hyper_planes_shutdown(1)) = 0;
                for g=1:1:size_hyper_planes_shed(1)
                    value_shed(g) = dot(random_belief,hyper_planes_shed(g,:));
                end
                for g=1:1:size_hyper_planes_non_shed(1)
                    value_non_shed(g) = dot(random_belief,hyper_planes_non_shed(g,:));
                end

                if (min(value_non_shed) <  min(value_shed))
                    old_add = 1 ; 
                else
                    old_add=2;
                end
            end

            
                if old_add ==1
                    if hyper_planes_non_shed_future ==0 %check if this hyperplane class is empty
                        hyper_planes_non_shed_future = hyper_planes(argm,:);
                    else
                            if BelongsTo(hyper_planes_non_shed_future,hyper_planes(argm,:)) == 0 %candidate hyper does not exist in future hyperplanes
                                   hyper_planes_non_shed_future(size_hyper_planes_non_shed_future(1) +1 , :) =hyper_planes(argm,:);
                            end
                    end
                else
                    if hyper_planes_shed_future ==0 %check if this hyperplane class is empty
                        hyper_planes_shed_future = hyper_planes(argm,:);
                    else
                        
                        if BelongsTo(hyper_planes_shed_future,hyper_planes(argm,:)) == 0 %candidate hyper does not exist in future hyperplanes
                                   hyper_planes_hed_future(size_hyper_planes_shed_future(1) +1 , :) =hyper_planes(argm,:);
                        end
                    end  
                end
        end
            
    end
    
    %concatinate the future hyperplanes
    
    if (length(hyper_planes_non_shed_future) ==1)  %if one of the future hyperplane classes is zero (meaning empty)
        hyper_planes_future =vertcat(hyper_planes_shed_future,hyper_planes_shutdown);
    elseif ((length(hyper_planes_shed_future) ==1))
        hyper_planes_future =vertcat(hyper_planes_non_shed_future,hyper_planes_shutdown);
    else  
        hyper_planes_future =vertcat(hyper_planes_non_shed_future,hyper_planes_shed_future,hyper_planes_shutdown);
    end
    
    size_hyper_planes_non_shed_future = size(hyper_planes_non_shed_future) ;
    size_hyper_planes_shed_future = size(hyper_planes_shed_future) ;
    size_hyper_shutdown = size(hyper_planes_shutdown) ;
    size_hyper_planes_future = size(hyper_planes_future);
    
    no_planes_future = size_hyper_planes_future(1);
    %stopping condition
    %stop when no belief can be improved
    counter = 0 ;
    
    
    for m =1:1:sample_number
        
        %calculate J_k for all samples
        J_k= dot(belief_samples(:,m),hyper_planes(1,:));
        if no_planes >1
            for k=2:1:no_planes
                if dot(belief_samples(:,m),hyper_planes(k,:)) <= J_k
                    J_k = dot(belief_samples(:,m),hyper_planes(k,:));
                end

            end
        end
        
        %calculate J_k+1 for all samples
        
        J_k1= dot(belief_samples(:,m),hyper_planes_future(1,:));
        if no_planes_future>1
            for k=2:1:no_planes_future
                if dot(belief_samples(:,m),hyper_planes_future(k,:)) <= J_k1
                    J_k1 = dot(belief_samples(:,m),hyper_planes_future(k,:));
                end

            end
        end
        
        if J_k1>J_k
            counter = counter+1;
        end
        %??????????????????????????????????????????????????????????????????????????????????????????????????????
    end
    
    if counter ==0
            break
    end
    %
 end
hyper_f_unshed=hyper_planes_non_shed_future
hyper_f_shed=hyper_planes_shed_future









end

