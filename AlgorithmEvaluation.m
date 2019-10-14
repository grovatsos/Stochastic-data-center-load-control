clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SIMULATION RESULTS

% THE CODE RIGHT NOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Cb=100; %breakdown cost
Cs_0=3; %shut down cost
Cs_1=20; %shut down cost
Cd_0=0.2;
Cd_1=0.2;
rho=0.01;
Power_obs_variance = 0.5;
Angle_obs_variance = 0.5;

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
    


%hyperplanes we are going to use for our decision rule

hyper_planes =[
  123.4000  123.4000  123.4000  123.4000  123.4000  123.4000  123.4000  123.4000    3.0000    3.0000    3.0000   20.0000   20.0000   20.0000
    9.3900    4.3570    7.1122    5.7853    8.8517    9.0372    5.7672    0.7145    8.2531    4.3165    4.2284    1.7112    4.2825    3.4388
    2.5666   29.1593    4.3333    4.5749   46.6880   15.7906    6.7152    3.8673   49.1418   16.7656   47.8468    8.7505   26.2523   28.8697];

 

   
  
  repetitions =1000
  cost(1:repetitions)=0;
 for  a=1:1:repetitions
   a
   
   
    size_hyper_planes = size(hyper_planes);
    no_planes = size_hyper_planes(1);
    %initialize the power gamma and actions at time 1 and zero
    for h = 1:1:100000

        if h==1  % VERIFY THE END FOR THIS IF IS CORRECT
            belief(1:14) = 0;
            belief(6)=1;
            Power = 0;
            gamma = 1;
            initial_action = 0; % action we take at time 0 (is unshed)
            %0 is for non_shed, 1 is for shed, 2 for shutdown
            Obs_z = normrnd(Power,Power_obs_variance);
            Theta = asin(Power / gamma);
            Obs_theta = normrnd( Theta , Angle_obs_variance); % do we need to add i.i.d. gaussian noise in each component
            action=[];
        else
            %calculate gamma at time h
            if gamma(h-1) == 0.5
                gamma =[gamma 0.5];
            else
                if binornd(1,rho)==1
                    gamma =[gamma 0.5];
                else
                    gamma =[gamma 1];
                end
            end

            %calculate the Power at time k

            if action(h-1) == 0 % if unshed at previous time instant //
                if ( Power(h-1) ==  -0.8 )
                     if (binornd(1,0.5) == 0) %power stays the same
                        Power = [Power Power(h-1)];
                     else
                        Power = [Power -0.4];
                     end   
                elseif ( Power(h-1) ==  0.8 )
                     if (binornd(1,0.5) == 0) %power stays the same
                        Power = [Power Power(h-1)];
                     else
                        Power = [Power 0.4];
                     end
                else
                    x = mnrnd(1,[1/3 1/3 1/3]);
                    if x(1) == 1
                        Power = [Power Power(h-1)+0.4];%% Power(h) = Power(h-1) + 0.4;
                    end
                    if x(2) == 1
                        Power = [Power Power(h-1)];
                    end
                    if x(3) == 1
                        Power = [Power Power(h-1)-0.4];
                    end
                end
            elseif action(h-1) == 1 %shed
                if ( Power(h-1) ==  -0.8 )
                    Power = [Power -0.4];
                elseif ( Power(h-1) ==  0.8 )
                    Power = [Power 0.4];
                elseif ( Power(h-1) ==  -0.4 )
                     if (binornd(1,0.5) == 0) %power stays the same
                        Power = [Power Power(h-1)];
                     else
                        Power = [Power 0];
                     end 
                 elseif ( Power(h-1) ==  0.4 )

                     if (binornd(1,0.5) == 0) %power stays the same
                        Power = [Power Power(h-1)];
                     else
                        Power = [Power 0];
                     end 
                else %Power(i-1) =0

                    x = mnrnd(1,[1/3 1/3 1/3]);
                    if x(1) == 1
                        Power = [Power Power(h-1)+0.4];
                    end
                    if x(2) == 1
                        Power = [Power Power(h-1)];
                    end
                    if x(3) == 1
                        Power = [Power Power(h-1)-0.4];
                    end

                end
            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %check for breakdown

            if ((gamma(h) == 0.5) &( abs(Power(h)) == 0.8  )) %find when the shutdown happens    
                breakdown = 1;
                cost(a) = cost(a) + Cb ; 
                break
            end

            %calculate angle

            Theta = [Theta asin(Power(h)/gamma(h))];

            %calculate the observations

            Obs_z = [Obs_z normrnd(Power(h),Power_obs_variance)];
            Obs_theta = [ Obs_theta  normrnd( Theta(h) , Angle_obs_variance)];

            %calculate the belief

            belief = [belief ; phi(belief(h-1,:),action(h-1)+1,Obs_z(h),Obs_theta(h),rho,Power_obs_variance,Angle_obs_variance,initial_action) ];
            %+1 has to do because phi wants the action to be in[1,3]

        end

        %decide which action to take
        minimizers(1:3)=0;
        for u=1:1:2
            step_size=0.05;
                Power_obs_range = [-1.6:step_size:1.6];%(-1,1) since my power is mostly around these values
                Angle_obs_range = [-1.6:step_size:1.6];%q = -1:0.001:1 %discritize the angle part of the integral
                for v= Power_obs_range
                    for q= Angle_obs_range 
                        %calculate the V(bao)

                        bao =  phi(belief(h,:),u,v,q,rho,Power_obs_variance,Angle_obs_variance,initial_action) ; 

                        J_k= dot(bao,hyper_planes(1,:));
                            if no_planes >1
                                for k=2:1:no_planes
                                    if dot(bao,hyper_planes(k,:)) <= J_k
                                        J_k = dot(bao,hyper_planes(k,:));
                                    end

                                end
                            end


                        minimizers(u) = minimizers(u) + (step_size^2)*psup(belief(h,:),u,v,q,Power_obs_variance,Angle_obs_variance,rho,initial_action )*J_k;

                    end
                end
                if u ==1 
                    minimizers(u) = minimizers(u) +dot(g_0,belief(h,:));
                elseif u==2
                    minimizers(u) = minimizers(u) +dot(g_1,belief(h,:));
                end



        end
        minimizers(3) = minimizers(3) +dot(g_2,belief(h,:)); %is that correct ?
        minimizers
        [w,action_to_take] = min(minimizers); 

        %differentiate between whether I am at unshed or shed because I can
        %take different actions at each case.

        if h==1
            previous_action = initial_action; %find the previous action
        else
            previous_action =action(h-1);
        end

        if previous_action == 0% If I am at unshed
            if (minimizers(2)<= minimizers(1))
                action= [action 1];
            else
                action= [action 0];
            end
        else % if I am at shed
            action= [action action_to_take-1];
                 %If i shutdown I need to add a  different cost depending on whether gamma = 1 or gamma =0.5
                 if action(h)==2
                        if gamma(h) == 1
                            cost(a) = cost(a) + Cs_1;
                        else
                            cost(a) = cost(a) + Cs_0;
                        end

                    break
                 end

        end

        %If i do not shutdown It means I may want to add Cd in the costt

        if (( gamma(h) == 1 ) & ( previous_action == 1 ))
            cost(a) = cost(a) + Cd_0;
        end

        if (( gamma(h) == 0.5 ) & ( previous_action == 1 ))
            cost(a) = cost(a) + Cd_1;
        end
     belief(h,:)   ;
     h       
     action(h)

    end
    
 performance(a)=sum(cost)/a
 end
final_performance=performance(repetitions)