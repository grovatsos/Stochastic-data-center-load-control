%this functionr returns Pk,Gamma_k and u_k-1 for the box j
%0 = unshed ,1 =shed
function [Power,Gamma,Action] = identify(j)
    if j == 1
        Power = -0.4;
        Gamma = 0.5;
        Action = 0 ;
    elseif j== 2
        Power = 0;
        Gamma = 0.5;
        Action = 0 ;
    elseif j== 3
        Power = 0.4;
        Gamma = 0.5;
        Action = 0 ;
    elseif j== 4
        Power = -0.8;
        Gamma = 1;
        Action = 0 ;
    elseif j== 5
        Power = -0.4;
        Gamma = 1;
        Action = 0 ;
    elseif j== 6
        Power = 0;
        Gamma = 1;
        Action = 0 ;
    elseif j== 7
        Power = 0.4;
        Gamma = 1;
        Action = 0 ;
    elseif j== 8
        Power = 0.8;
        Gamma = 1;
        Action = 0 ;
    elseif j== 9
        Power = -0.4;
        Gamma = 0.5;
        Action = 1;
    elseif j== 10
        Power = 0;
        Gamma = 0.5;
        Action = 1;
    elseif j== 11
        Power = 0.4;
        Gamma = 0.5;
        Action = 1;
    elseif j== 12
        Power = -0.4;
        Gamma = 1;
        Action = 1;
    elseif j== 13
        Power = 0;
        Gamma = 1;
        Action = 1;
    else
        Power = 0.4;
        Gamma = 1;
        Action = 1;
    end
end