%% Multiple effect evaporator algorithm for 'n' effects

% This is just a basic and amateur algorithm used for an Unit Operations discipline project in a Chemical Enginnering course
%
% I know it's very limited and basic with a lot of programming issues, but it get the job done (sorta).
% Feel free to make improvements/additions.
%
% This algorithm is based on the one present at the chapter 8 of the book:
% Transport Processes and Separation Process Principles, GENKOPLIS, Forth Edition, 2003, page 527

%% Known Bugs/Limitations:
% 1 - Low concentration of solid fraction at the feed (<0.1%) generates absurd values
% 2 - Breaks at 1 effect (n = 1). Which is equivalent to one effect.
% 3 - Only applies to coffee extracts (which was my project, by the way) 
% 4 - Ideal thermodinamic models only 
% 5 - SI units only 

format long
clear all
clc

global L V n

%% Project Variables
% Matrix allocation
n   = input('How many effects to use: ') ;   %# of effects
V   = zeros(n+1,5);  % Vapor flow
L   = zeros(n+1,5);  % Liquid flow
DT  = zeros(1,n);    % Temperature difference along the effects
U   = zeros(1,n);    % Global Heat Transfer Coefficients
Q = zeros(1,n);      % Heat transfer at the effects

%Antoine Parameters
At = 11.6834;
Bt = 3816.44;
Ct = -46.13;

%Steam heat capacity coefficients 
Av = 3.47;
Bv = 1.45E-3;
Dv = 0.121E5;

L(1,1)    = input('Feed current (kg/h): ');              % Feed mass flow at first effect
L(1,2)    = input('Solid fraction at feed (m/m): ');     % Solid fraction at feed
L(n+1,2)  = input('Solid fraction at output(m/m): ');    % Solid fraction at output (last effect)
L(1,3)    = input('Feed temperature (ºC): ');            % Feed flow temperature
V(1,2)    = input('Heating vapor temperature(ºC): ');    % Heating vapor temperature
V(n+1,2)  = input('Condensed vapor temperature (ºC): '); % Condensed vapor at the last effect  
Tst = 100; 

flag1 = input('Do you like to inform the Global Heat Coeficients manually for each effect? (s/n): ','s'); # Manual global heat transfer coefficients insertion
flag1 = lower(flag1);

%% Global Heat Coefficients

if (flag1 ~= 's')
    stp = 2500 - 1000;
    stp = stp/n;
    for (k = 1:n)
        U(k) = 2500 - (k-1)*stp;
    end
else 
    for (k = 1:n)
        fprintf("Coefficient for the %i th effect",k);
        U(k) = input(': ');
    end
end

iU = 0;
for (i = 1:n);
    iU = iU + (1/U(i));
end

%% Effect temperature first approximation  
DeltaT = V(1,2) - V(n+1,2);
for (i=1:n)
    DT(i) = (DeltaT*(1/U(i)))/iU;
end

%% Initial Mass Balance
L(n+1,1) = L(1,1)*L(1,2)/(L(n+1,2));
Vt = (L(1,1) - L(n+1,1))/n;

for (i = 2:n)
    L(i,1) = L(i-1,1) - Vt;
    L(i,2) = (L(i-1,1)*L(i-1,2))/L(i,1);
end

V(:,1) = Vt;

flag2 = 0;

%% Main Loop
% Since it's an iterative algorithm, it needs to run multiple times to make the effect heat transfer area converge. 
% I know, I could just make a while loop and set a threshold for the difference of the areas among the iterations.
% It was my first by program, so take easy :|

for i = 1:5
    %% Vapor current temperature (except 1 and 'n'th)

    for (i = 1:n)
        V(i+1,2) = V(i,2) - DT(i);
    end
   
    if (V(2,2) > Tst || flag2 == 1)
       
        flag2 = 1;
        V(2,2) = Tst;
        DT(1) = V(1,2) - V(2,2);
        DeltaT = V(2,2) - V(n+1,2);
        
        iU = 0;

        for (i = 2:n)
            iU = iU + (1/U(i));
        end
        
        for (i=2:n)
            DT(i) = (DeltaT*(1/U(i)))/iU; 
        end    
    end
        
        
    %% Saturation Pressure
    for (i = 2:n+1)
        V(i,5) = 1E2 * exp(At - (Bt/((V(i,2))+273+Ct)));
        L(i,5) = 0.8474E-2 * (L(i,2)*100)^0.9895 * exp(2.570E-2 * (L(i,2)*100))*V(i,5)^0.1163;
        L(i,3) = V(i,2) + L(i,5);
    end
 
    %% Latent heat and entalphy
    for i = 1:n+1
        L(i,4) = (L(i,3)*(1549*L(i,2)+4176*(1-L(i,2)))+((L(i,3)^2)/2)*(1.96*L(i,2)-0.0909*(1-L(i,2)))...
            +((L(i,3)^3)/3)*(-0.00549*L(i,2)+0.00547*(1-L(i,2))))/1000; % kJ/kg
    end
        
    %% Steam Heat Capacity 
    for i = 2:n+1
        V(i,3) = 2501.3 + 0.4618 * (Av*(V(i,2)) + (Bv/2)*((V(i,2)+273)^2-273^2) - Dv*((V(i,2)+273)^-1 - 273^-1 )); %kJ/kg
    end
    
    %% Latent Heat for Heating Steam 
    for i = 1:n
       V(i,4) =  2257 * ((1 - (V(i,2)+273)/(647.1))/(1-0.577))^0.38;
    end
    
    %% First Values Approximation for fsolve function
    x0 = zeros(2*n,1);
    
    for i = 1:n+1
        x0(i) = V(i,1);
    end
    
    for i = n+2:2*n
        x0(i) = L(i-n,1);
    end
    
    %options = optimoptions('fsolve','MaxIterations',2500)
    %options = optimoptions('fsolve','FunctionTolerance',1e-12)
    
    [x, fval, exitflag] = fsolve(@evaporator, x0, [] , L, V, n);
      
    %% The result from the fsolve function (the vapor and liquid mass flow) is atributed to the matrixes    
    for i = 1:n+1
        V(i,1) = x(i);
    end
    
    for i = n+2:2*n
        L(i-n,1) = x(i);
    end
    
    %% New mass balance based on the new data
    for (i = 2:n+1)
        L(i,1) = L(i-1,1) - V(i,1);
        L(i,2) = (L(i-1,1)*L(i-1,2))/L(i,1);
    end
    
    %% Heat and effect heat transfer area calculation
    for (i = 1:n)
        Q(i) = V(i,4)*V(i,1)*(1/3600);
        A(i) = 1000*Q(i)/((DT(i) - L(i+1,5))*U(i));
    end
    
    %% DeltaT recalculation
    if (flag2 == 0)
        
        Am = sum(A)/n;
        
        for (i = 1:n)
             DT(i) = DT(i)*A(i)/Am;
        end
    else
        
        Am = sum(A(2:n))/(n-1);
        
        for (i = 2:n)        
            DT(i) = DT(i)*A(i)/Am;
            
        end
    end
        
end

fprintf("Generating log file...\n\n")
logname = strcat('evap-log-',datestr(datetime('now','InputFormat','ddMMyyyy')),'.txt');

fileID = fopen(logname,'w');

fprintf(fileID,"------ LOG ------ \n\n");

fprintf(fileID,"--PROJECT PARAMETERS--\n");
fprintf(fileID,"Feed - %2.2f kg/h \n", L(1,1));
fprintf(fileID,"Solid fraction at feed - %2.2f \n", L(1,2));
fprintf(fileID,"Solids at the output - %2.2f \n", L(n+1,2));
fprintf(fileID,"Steam temperatura - %2.2f ºC \n", V(1,2));
fprintf(fileID,"Condensed steam at output - %2.2f ºC \n", V(n+1,2));
fprintf(fileID,"Heating steam temperature - %2.2f kg/ħ \n", V(1,1));
fprintf(fileID,"Economy - %2.4f \n\n",(sum (V(2:n+1,1)))/V(1,1));


for (i = 2:n+1)
   
    fprintf(fileID,"-- %iº EFFECT -- \n", i-1);
    fprintf(fileID,"Global Heat Transfer Coefficient - %2.2f W/m².ºC \n", U(i-1));
    fprintf(fileID,"Liquid flow - %2.2f kg/h \n", L(i,1));
    fprintf(fileID,"Steam flow - %2.2f kg/h \n", V(i,1));
    fprintf(fileID,"Liquid flow exit temperature - %2.2f ºC \n", L(i,3));
    fprintf(fileID,"Steam flow exit temperature - %2.2f ºC \n", V(i,2));
    fprintf(fileID,"Solid fraction - %2.5f (m/m) \n", L(i,2));
    fprintf(fileID,"Effect heat transfer area - %2.2f m² \n\n", A(i-1));
    
end

fclose(fileID);
fprintf("Voilá!");

function G = evaporator(x, L, V, n) 
     
    G = zeros(2*n,1);

    for i = 1:n
        if(i == 1)
            G(1) = x(n+2)*L(2,4) + x(2)*V(2,3) - L(1,1)*L(1,4) - x(1)*V(1,4);
            G(2) = x(2) + x(n+2) - L(1,1);
        elseif (i == n)
            G(2*n-1) = L(n+1,1)*L(n+1,4) + x(n+1)*V(n+1,3) - x(2*n)*L(n,4) - x(n)*V(n,4);
            G(2*n) = L(n+1,1) + x(n+1) - x(2*n);
        else
            G(2*i-1) = x(i+1+n)*L(i+1,4) + x(i+1)*V(i+1,3) - x(n+i)*L(i,4) - x(i)*V(i,4);
            G(2*i) = x(i+1) + x(i+1+n) - x(i+n);
        end
    end

end

% Doing some random modifications to test git resources
% Added another comment line for learning 

