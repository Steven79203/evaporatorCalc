%% Protótipo - Algorítmo pra cálculo de evaporadores de 'n' efeitos

%% Observações
% 1 - Frações de sólidos muito baixas na alimentação geram erros absurdos
% (x < 0.1)
% 2 - Adicionar apenas um efeito (n = 1), o que é equivalente a um evaporador simples,
% termina o programa em erro. 
% 3 - Aplica-se apenas a soluções de extrato de café
% 4 - 


format long
clear all
clc

global L V n

%% Variáveis iniciais de projeto

n   = input('Informe o número de efeitos: ') ;         %Nº de efeitos
V   = zeros(n+1,5);
L   = zeros(n+1,5);
DT  = zeros(1,n);
U   = zeros(1,n);
Q = zeros(1,n);

% L(1,1)    = 50000;     %Vazão mássica de alimentação
% L(1,2)    = 0.01;  	   %fração de sólidos na corrente F
% L(n+1,2)  = 0.6;       %Fração de sólidos na saída
% L(1,3)    = 50;        %Temperatura de F
% V(1,2)    = 100;       %Temperatura do vapor vivo
% V(n+1,2)  = 50;        %Temperatura do condensado no 'n-ésimo' efeito

L(1,1)    = input('Corrente de entrada kg/h: ');     
L(1,2)    = input('Fração de sólidos na entrada: '); 
L(n+1,2)  = input('Fração de sólidos na saída: ');    
L(1,3)    = input('Temperatura de entrada ºC: ');       
V(1,2)    = input('Temperatura do vapor vivo ºC: ');    
V(n+1,2)  = input('Temperatura do condensado ºC: ');    
Tst = 100; %% Temperatura do primeiro efeito

flag1 = input('Informar manualmente os coeficientes globais? (s/n): ','s');
flag1 = lower(flag1);

%% Matriz (n+1)x(m)
% 
%           1   2   3    4     5
% L  = [] % L, XM, TL,   HL,  EPE (5xn)
% V  = [] % V, TV, HV, lambda, P  (5xn)

%% Coeficientes globais 

if (flag1 ~= 's')
    
    stp = 2500 - 1000;
    stp = stp/n;
    
    for (k = 1:n)
        U(k) = 2500 - (k-1)*stp;
    end

else 
    for (k = 1:n)
        fprintf("Coeficiente global do %iº efeito",k);
        U(k) = input(': ');
    end
end
    
%% Somatório do inverso dos coeficiente U.

iU = 0;

for (i = 1:n);
    iU = iU + (1/U(i));
end


%% Distribuição da temperaturas entre os efeitos

DeltaT = V(1,2) - V(n+1,2);

for (i=1:n)
    DT(i) = (DeltaT*(1/U(i)))/iU;
end
 

%% Balanço de massa inicial

L(n+1,1) = L(1,1)*L(1,2)/(L(n+1,2));
Vt = (L(1,1) - L(n+1,1))/n;

for (i = 2:n)
    L(i,1) = L(i-1,1) - Vt;
    L(i,2) = (L(i-1,1)*L(i-1,2))/L(i,1);
end

V(:,1) = Vt;

flag2 = 0;

%% INICIO DO LOOP( FOR OR WHILE )
for i = 1:5
    %% Temperatura das correntes de vapor (exceto 1 e n-ésima)

    for (i = 1:n)
        V(i+1,2) = V(i,2) - DT(i);
    end
    
    %% Caso TV2 > 100, TV2 é fixo em 100ºC e as temperaturas dos efeitos
    %% são recalculadas levando em conta TV2 e TVn+1 
    
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
        
        
    %% Pressões de saturação, EPE e nova temp. 
    
    %Parâmetros de Antoine
    At = 11.6834;
    Bt = 3816.44;
    Ct = -46.13;
    
    for (i = 2:n+1)
        V(i,5) = 1E2 * exp(At - (Bt/((V(i,2))+273+Ct)));
        L(i,5) = 0.8474E-2 * (L(i,2)*100)^0.9895 * exp(2.570E-2 * (L(i,2)*100))*V(i,5)^0.1163;
        L(i,3) = V(i,2) + L(i,5);
    end
 

   %% Entalpias e calores latentes das correntes
   
    for i = 1:n+1
        
        L(i,4) = (L(i,3)*(1549*L(i,2)+4176*(1-L(i,2)))+((L(i,3)^2)/2)*(1.96*L(i,2)-0.0909*(1-L(i,2)))...
            +((L(i,3)^3)/3)*(-0.00549*L(i,2)+0.00547*(1-L(i,2))))/1000; % kJ/kg
    end
    
    % Parâmetros de capacidade calorífica do vapor d'água
    Av = 3.47;
    Bv = 1.45E-3;
    Dv = 0.121E5;
    
    for i = 2:n+1
        
        V(i,3) = 2501.3 + 0.4618 * (Av*(V(i,2)) + (Bv/2)*((V(i,2)+273)^2-273^2) - ...
                Dv*((V(i,2)+273)^-1 - 273^-1 )); %kJ/kg
    end
    
    %% Calor latente (Correntes 1:n)
    
    for i = 1:n
       V(i,4) =  2257 * ((1 - (V(i,2)+273)/(647.1))/(1-0.577))^0.38;
    end
    
    %% Vetor "chute inicial" p/ fsolve
    
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
      
    %% Realocação das variáveis
    
    for i = 1:n+1
        V(i,1) = x(i);
    end
    
    for i = n+2:2*n
        L(i-n,1) = x(i);
    end
    
    %%  Redistribuição de massa
    
    for (i = 2:n+1)
        L(i,1) = L(i-1,1) - V(i,1);
        L(i,2) = (L(i-1,1)*L(i-1,2))/L(i,1);
    end
    
    %% Cálculo das energias e das áreas
    
    for (i = 1:n)
        Q(i) = V(i,4)*V(i,1)*(1/3600);
        A(i) = 1000*Q(i)/((DT(i) - L(i+1,5))*U(i));
    end
    
    
    %% Redistribuição das DTs
    
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

fprintf("Cálculos finalizados. Gerando relatório...\n\n")
logname = strcat('evap-log-',datestr(datetime('now','InputFormat','ddMMyyyy')),'.txt');

fileID = fopen(logname,'w');

fprintf(fileID,"------ RELATORIO ------ \n\n");

fprintf(fileID,"--PARÂMETROS DE PROJETO--\n");
fprintf(fileID,"Alimentação - %2.2f kg/h \n", L(1,1));
fprintf(fileID,"Sólidos na entrada - %2.2f \n", L(1,2));
fprintf(fileID,"Sólidos na sáida - %2.2f \n", L(n+1,2));
fprintf(fileID,"Temperatura do vapor vivo - %2.2f ºC \n", V(1,2));
fprintf(fileID,"Temperatura do condensado - %2.2f ºC \n", V(n+1,2));
fprintf(fileID,"Vazão de vapor vivo - %2.2f kg/ħ \n", V(1,1));
fprintf(fileID,"Economia - %2.4f \n\n",(sum (V(2:n+1,1)))/V(1,1));


for (i = 2:n+1)
   
    fprintf(fileID,"-- %iº EFEITO -- \n", i-1);
    fprintf(fileID,"Coeficiente global - %2.2f W/m².ºC \n", U(i-1));
    fprintf(fileID,"Vazão da corrente líquida - %2.2f kg/h \n", L(i,1));
    fprintf(fileID,"Vazão da corrente de vapor - %2.2f kg/h \n", V(i,1));
    fprintf(fileID,"Temperatura do líquido na saída - %2.2f ºC \n", L(i,3));
    fprintf(fileID,"Temperatura do vapor na saída - %2.2f ºC \n", V(i,2));
    fprintf(fileID,"Fração de sólidos na sáida - %2.5f (m/m) \n", L(i,2));
    fprintf(fileID,"Área de troca térmica - %2.2f m² \n\n", A(i-1));
    
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
