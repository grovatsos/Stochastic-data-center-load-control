
%[ (3 states) uk_1 = NS , gamma_k =0.5 | (5 states) uk_1 = NS , gamma_k
%=1|(3 states) uk_1 = SH , gamma_k =0.5|(3 states) uk_1 = SH , gamma_k
%=1]

function [hyper_planes_future] = FullShedPerseus(hyper_planes,belief,action_at_time_zero)

    belief_samples = belief; %belief_samples = B', belief =B
    %hyper_planes =[123.4000  123.4000  123.4000  123.4000  123.4000  123.4000  123.4000  123.4000    3.0000    3.0000    3.0000   20.0000   20.0000   20.0000
    %1.1827    2.1654    3.1481    4.1308    5.1135    6.0962    7.0789    8.0616    9.0443   10.0270   11.0097   11.9924   12.9751   13.9578
    %1.0000    2.0000    3.0000    4.0000    5.0000    6.0000    7.0000    8.0000    9.0000   10.0000   11.0000   12.0000   13.0000   14.0000];

    Power_obs_variance = 0.5; %variance of the power observations, this may change
    Angle_obs_variance = 0.5; %variance of the angle observations, this may change
    rho = 0.01;  %parameter for gamma going from 1 to 0.5
    Cb=100; %breakdown cost
    Cs_0=3; %shut down cost
    Cs_1=20; %shut down cost
    Cd_0=0.2;
    Cd_1=0.2;

    g_0(1:14)=0;
    g_1(1:14)=0;
    g_2(1:14)=0;


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

     hyper_planes_future=0;
     hyper_planes_shutdown = [Cd_0+Cd_1+Cb+Cs_0+Cs_1 Cd_0+Cd_1+Cb+Cs_0+Cs_1 Cd_0+Cd_1+Cb+Cs_0+Cs_1 Cd_0+Cd_1+Cb+Cs_0+Cs_1 Cd_0+Cd_1+Cb+Cs_0+Cs_1 Cd_0+Cd_1+Cb+Cs_0+Cs_1 Cd_0+Cd_1+Cb+Cs_0+Cs_1 Cd_0+Cd_1+Cb+Cs_0+Cs_1 Cs_0 Cs_0 Cs_0 Cs_1 Cs_1 Cs_1]; % each row is a hyperplane. Start with a single random hyperplane for SHUTDOWN action

     while length(belief_samples) ~= 0    
        sample_number =size(belief_samples,2) ; 
        random_selection = randi([1 sample_number],1,1);
        random_belief = belief_samples(:,random_selection); % sample a random belief out of all sampled beliefs
        size_hyper_planes = size(hyper_planes) ;
        no_planes= size_hyper_planes(1) ; % determine the total number of current hyper_planes(all actions summed)
        %====================================================================
        right_term(1:3,1:14)=0; %start calculating a_u^p where the first coordinate is for the action
        for u = 1:1:2 % go through all actions // meaning calculate a_u^p for every u // here we have on action continue
          % u=1 is for nonshed /// u=2 for shed // shutdown hyperplane if fixed 
            %calculate the integral by treating it like a sum
            step_size=0.1;
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
%             v=NaN;
%             q=NaN;%???????????????????????????????????????????????????????????????????????????????????????????????
%             right_term(u,:)  = right_term(u,:)  +  (step_size^2)*psup(random_belief,u,v,q,Power_obs_variance,Angle_obs_variance,rho,action_at_time_zero )*0 ; % SInce the min is equal to Cb in this case
            if u ==1 
                right_term(u,:) = right_term(u,:) +g_0;
            else
                right_term(u,:) = right_term(u,:) +g_1;
            end
        end
        right_term(3,:) = hyper_planes_shutdown;
        %==========================================================================================
        %find the new candidate hyperplane to add
        %the candidates are between the non-shed and the shed since the
        %shutdown hyperplane is fixed
         hyper_1 = dot(random_belief,right_term(1,:)); 
         hyper_2 = dot(random_belief,right_term(2,:));
         hyper_3 = dot(random_belief,hyper_planes_shutdown);
         if  ((hyper_1 <= hyper_2)&(hyper_1 <= hyper_3))   %we have a candidate hyper for non-shed action
             canditate_hyper = right_term(1,:);
         elseif ((hyper_2 <= hyper_1)&(hyper_2 <= hyper_3))
             canditate_hyper = right_term(2,:);
         else
             canditate_hyper = right_term(3,:);
         end
         canditate_hyper
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
        size_hyper_planes_future = size(hyper_planes_future);
        if dot(random_belief,canditate_hyper) <= J_k
            if hyper_planes_future ==0 %check if this hyperplane class is empty
                   hyper_planes_future = canditate_hyper;
            else
                   if BelongsTo(hyper_planes_future,canditate_hyper) == 0 %candidate hyper does not exist in future hyperplanes
                       hyper_planes_future(size_hyper_planes_future(1) +1 , :) = canditate_hyper;
                   end
            end
            %remove all belief points whole value was already improved by
            %the new hyperplane
            f=1;
            while f<= sample_number
                J_k= dot(belief_samples(:,f),hyper_planes(1,:));
                if no_planes >1
                    for k=2:1:no_planes
                        if dot(belief_samples(:,f),hyper_planes(k,:)) <= J_k
                            J_k = dot(belief_samples(:,f),hyper_planes(k,:));
                        end%

                    end
                end
                if dot(belief_samples(:,f),canditate_hyper) <= J_k
                    belief_samples(:,f)=[];
                    sample_number =size(belief_samples,2) ;  
                else
                    f=f+1;
                end
            end 
        else 
            size_hyper_planes_future = size(hyper_planes_future);
            no_planes_future = size_hyper_planes_future(1);
            belief_samples(:,random_selection) = [];
            sample_number =size(belief_samples,2) ; 
            %else add the hyperplane that maximizes dot(random_belief,hyperplane)
            for i =1:1:no_planes
                maxz(i) = dot(random_belief,hyper_planes(i,:));
                [w,argm] = min(maxz); 
            end    
                %find to which action this hyperplane corresponds
            if hyper_planes_future ==0 %check if this hyperplane class is empty
                   hyper_planes_future = hyper_planes(argm,:);
            else
                if BelongsTo(hyper_planes_future,hyper_planes(argm,:)) == 0 %candidate hyper does not exist in future hyperplanes
                    hyper_planes_future(size_hyper_planes_future(1) +1 , :) =hyper_planes(argm,:);
                end
            end
        end
        size_hyper_planes_future = size(hyper_planes_future);
        no_planes_future = size_hyper_planes_future(1);
        sample_number
     end  
end

