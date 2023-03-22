#!/usr/bin/env python 

import numpy as np
from scipy.optimize import fsolve

def evaporator(x,L,V,n):
    G = np.zeros(2*n,1)

    for i in range(n):
        if i == 0:
            G[0] = x[n+2] * L[1,3] * x[1]*V[1,2] - L[0,0]*L[0,3] - x[0]*V[0,3]
            G[1] = x[1] + x[n+2] - L[0,0]
        elif i == n:
            G[2*n-1] = L[n+1,0]*L[n+1,3] + x[n+1]*V[n+1,2] - x[2*n]*L[n,3] - x[n]*V[n,3]
            G[2*n] = L[n+1,0] + x[n+1] - x[2*n]
        else
            G[2*i-1] = x[i+1+n]*L[i+1,4] + x[i+1]*V[i+1,2] - x[n+1]*L[i,3] - x[i]*V[i,3]
            G[2*n] = x[i+1] + x[i+1+n] - x[i+n]

    return G


n = int(input("Nº de efeitos a ser usado"))

V = np.zeros(n+1, 5)
L = np.zeros(n+1, 5)
DT = np.zeros(1,n)
U = np.zeros(1,n)
Q = np.zeros(1,n)


# Parâmetros da Antoine da água
At = 11.6834
Bt = 3816.44
Ct = -46.13

# Coeficientes para cálculo da capacidade calorífica
Av = 3.47
Bv = 1.45E-3
Dv = 0.121E5

# Vazão de alimentação (kg/h)
L(0,0) = float(input("Alimentação em kg/g"))
L(0,1) = float(input("Fração de sólidos na corrente de entrada (m/m)"))
L(n,1) = float(input("Fração de sólidas deseja na saída (m/m)"))
L(0,2) = float(input("Temperatura da corrente de entrada (ºC)"))

V(0,1) = float(input("Temperatura do vapor de aqueciment (ºC)"))
V(n+1,1) = float(input("Temperatura do condensado final (ºC)"))
Tst = 100

# Definir coeficientes globais de transferência 
set_hc = lower(input("Informar os coeficientes globais de transferência de calor? [s/n]"))

if (flag == 's'):
    step = (2500 - 1000)/n
    for k in range(n):
        U[k] = 2500 - (k-1)*step
else:
    for k in range(n):
        print("Coeficiente no estágio",k)
        U[k] = float(input(":"))


iU = 0
for i in range(n):
    iU = iU + (1/U[i])

# Primeira aproximação das temperaturas dos efeitos
DeltaT = V(0,1) - V(n+1,1)
for i in range(n):
    DT[i] = (DeltaT*(1/U[i]))/iU


# Balanço de massa inicial
L[n+1,0] = L[0,0]*L[1,1]/L[n+1,1]
Vt = L[0,0] - L[n+1,0]/n

for i in range(2,n):
    L[i,0] = L[i-1,0] - Vt
    L[i,1] = (L[i-1,0]*L[i-1,1])/L[i,0]

V[:,1] = Vt

flag2 = 0


for i in range(5):
    for i in range(n):
        V[i+1,1] = V[i,1] - DT[i]

    if (V(1,1) > Tst) or (flag2 == 1):
        flag2 = 1
        V[1,1] = Tst
        DT[0] = V[0,1] - V[1,1]
        DeltaT = V[1,1] - V[n,2]

        iU = 0

        for i in range(2,n):
            iU = iU + (1/U[i])

        for i in range(2,n):
            DT[i] = (DeltaT*(1/U[i]))/iU

    # Pressão de saturação
    for i in range(2:n):
        V[i,4] = 1E2 * np.exp(At - (Bt/((V[i,1]) + 273 + Ct)))
        L[i,4] = 0.847E-2 * (L[i,2]*100)**0.9895 * exp(2.570E-2 * (L[i,1]*100)) * V[i,5]**0.1163
        L[i,2] = V[i,1] + L[i,4]
    
    # Entalpia e Calor latente
    for i in range(n):
        L[i,3] = (L[i,3]*(1549*L[i,1]+4176*(1-L[i,1]))+((L[i,2]**2)/2)*(1.96*L[i,1]-0.0909*(1-L[i,2])) + ((L[i,2]**3/3)*(-0.00549*L[i,2]+0.00547*(1-L[i,1]))))/1000


    # Capacidade calorífica do vapor
    for i in range(n):
        v[i,3] = 2257 * ((1-(V[i,1]+273)/647.1))/(1-0.577))**0.38


    # Valores de aproximação para o fsolve
    x0 = np.zeros(2*n,1)

    for i in range(n+1):
        x0[i] = V[i,0]

    for i in range(n+2,2*n):
        x0[i] = L[i-n],1]



    sol = fsolve(evaporator, x0, args=(L,V,n))

    for i in range(n):
        V[i,0] = sol.x[i]

    for i in range(n+2:2*n):
        L[i-n,0] = x[i]

    # Novo balanço de massa 
    for i in range(2:n):
        L[i,0] = L[i-1,0] - V[i,0]
        L[i,1] = (L[i-1,1]*L[i-1,1))/L[i,0]

    # Calor e área de transfência efetiva de calor
    for i in range(n):
        Q[i] = V[i,3] * V[i,0] * (1/3600)
        A[i] = 1000*Q[i]/((DT[i]-L[i+1,4])*U[i])


    #DeltaT recalculation
    if flag2 == 0:
        Am = np.sum(A)/n
        
        for i in range(n):
            DT[i] = DT[i] * A[i]/Am

    else:
        Am = np.sum(A[2:n])/(n-1)

        for i in range(2,n):
            DT[i] = DT[i] * A[i]/Am

print("Gerando arquivo de Log")

log_name = "evaporator_log"

log_file = open(logname,'w')

log_file.write("--------LOG----------\n\n")

params = f"""
--PARÂMETROS DE PROJETO--
-> Alimentação: {L[0,0]}
-> Fração de sólidos na alimentação: {L[0,1]}
-> Sólidos na saída: {L[n,1]}
-> Temperatura do vapor: {V[0,1]}
-> Vapor condensado na śaida: {V[n,1]}
-> Temperatura do vapor de aquecimento {V[0,0]}
"""

log_file.write(params)


for i in range(2,n):
    effect = f"""
    -- EFEITO Nº {i} --
    Coeficiente global        - {U[i-1]}
    Vazão de líquido          - {L[i,0]}
    Vazão de vapor            - {V[i,0]}
    Temperatura de L na saída - {L[i,2]}
    Temperatura de V na saída - {V[i,1]}
    Fração de sólidos         - {L[i,1]}
    Área de transferência     - {A[i-1]}
    \n
    """

    log_file.write(effect)

log_file.close()


