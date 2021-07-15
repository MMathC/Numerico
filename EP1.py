# -*- coding: utf-8 -*-
"""
Created on Sat Jun 12 19:36:58 2021

@author: matheus
"""
import matplotlib.pyplot as plt
import numpy as np
import math

def autovalor(matrizT): # retorna os autovalores de uma matriz triangular
    autovalores = []
    sizeLine = np.shape(matrizT)[0]
    sizeColumn = np.shape(matrizT)[1]
    while sizeLine >= 1:
        while sizeColumn >= 1:
            if matrizT[-sizeLine][-sizeColumn] != 0:
                autovalores.append(matrizT[-sizeLine][-sizeColumn])
                sizeLine-=1
                sizeColumn-=1
    return autovalores

def det(autovalores):
    det = 1
    for i in range(len(autovalores)):
        det *= autovalores[i]
    return det


def identidade(n): # retorna a matriz identidade com ordem n x n 
    V = []
    v = []
    line = 0
    column = 0
    for i in range(n):
        for j in range(n):
            if line == i and column == j and line == column:
                v.append(float(1))
                column+=1
            else:
                v.append(0)
        V.append(v)
        v = []
        line+= 1
    return np.array(V)


def tGivens (A, i, j, k):
    alfak = A[i][k]
    betak = A[j][k]
    ck = alfak/math.sqrt((alfak**2)+(betak**2))
    sk = -betak/math.sqrt((alfak**2)+(betak**2))
    return ck , sk

def ajuste(A, i, j, k, ck, sk):
    lp = 1
    while lp >= 0:
        if i == k:
            A[i][k] = ck
            i+=1
        if j!=k:
            if k>j:
                A[j][k] = -sk
            else:
                A[j][k] = sk
                k += 1
                j-=1
        lp -= 1
    return A, i-1, k

def sgn(d): # Ajuste do sinal
    if d >= 0:
        return 1
    else:
        return -1
    
def heuristicaWilkinson(alfan_1, alfan, betan_1, desloc): # heurística de Wilkinson
    if desloc == 's':
        dk = (alfan_1 - alfan)/2
        mik = alfan + dk - sgn(dk) * math.sqrt(dk**2 + betan_1**2)
        return mik
    return 0

def deslocamento(mik,n): # matriz deslocamento
    return mik * identidade(n)
    
def rotGivens(A, Vi, desloc):
    n = len(A)
    Qr = identidade(n)
    mik = heuristicaWilkinson(A[n-2][n-2],A[n-1][n-1],A[n-1][n-2],desloc)
    A = A - deslocamento(mik * identidade(n),n)
    ck, sk = tGivens(A,0,1,0)
    auxi = 0
    auxj = 0
    Q1 = identidade(n)
    Q1, auxi, auxj = ajuste(Q1, auxi, auxi+1, auxj, ck, sk)
    R = Q1@A
    Q_ts = Q1.T
    for r in range(n-2):
        ck, sk = tGivens(R, auxi, auxi+1, auxj)
        Qr, auxi, auxj = ajuste(Qr, auxi, auxi+1, auxj, ck, sk)
        Q_ts = Q_ts @ Qr.T
        R = Qr@R
        Qr = identidade(n)
    autovetores = Vi @ Q_ts
    print("autovetores: \n",autovetores,"\n")
    Ak = (R @ Q_ts) + deslocamento(mik * identidade(n),n)
    autovalores = autovalor(Ak)
    return Ak, autovalores, autovetores

def diminui(A): # diminui a ordem da matriz em 1. ex: A -> n x n ; Anova -> n-1 x n-1
    Anova = identidade(len(A)-1)
    for i in range(len(A)-1):
        for j in range (len(A)-1):
            Anova[i][j] = A[i][j]
    return Anova

def diminuivet(A): # diminui o tamanho do vetor em 1. ex A -> n ; Anova -> n-1
    Anova = np.array([])
    for i in range(len(A)-1):
        Anova = np.append(Anova,A[i])
    return Anova

def tiraColuna(A):
    Anova = []
    for i in range(len(A)):
        Anova.append([])
        for j in range(len(A[0])-1):
            Anova[i].append(A[i][j])
    return np.array(Anova)

def verificaBeta(A): # verifica o beta menos 1 da matriz parametro
    if abs(A[len(A)-1][len(A)-2]) < 1e-6:
        #print("beta_: ",A[len(A)-1][len(A)-2])
        Aux = diminui(A)
        menor = True
        return menor, Aux, A
    else:
        menor = False
        return menor, A, A

def matrizPrincipal(Anova, Amain): # Atualiza a matriz principal com as matrizes auxiliares diminuidas
    if len(Anova)==len(Anova[len(Anova)-1]):
        for i in range(len(Anova)):
            for j in range(len(Anova[0])):
                Amain[i][j] = Anova[i][j]
    else:
        for i in range(len(Anova)):
            for j in range(len(Anova[len(Anova)-1])):
                Amain[i][j] = Anova[i][j]
    return Amain

def vetorPrincipal(Anova, Amain): # Atualiza o vetor principal com os vetores auxiliares diminuidos
    for i in range(len(Anova)):
        Amain[i] = Anova[i]
    return Amain

def calculosAnaliticos(n): # Realiza o calculo analítico dos autovalores e autovetores
    lambdavec = []
    vj = identidade(n)
    for j in range(-n,0):
        lambdaj = 2*(1-math.cos((-j)*math.pi/(n+1))) # Auto-valores analiticos
        lambdavec.append(lambdaj)
    for j in range (n):
        for i in range (n):
            vj[i][j] = math.sin(abs(i+1)*abs(j+1)*math.pi/(n+1)) # Auto-vetores analiticos
    return lambdavec, normalizacao(vj)

def ordemAutovetores(autovetores): # Altera a ordem dos autovetores para estarem de acordo com a ordem dos autovalores
    n = len(autovetores)
    matrizAux = identidade(n)
    for i in range(n):
        for j in range(-n,0):
            matrizAux[i][j] = autovetores[i][((-j)-1)]
    return matrizAux
    
def QR(A, desloc):
    n = len(A)
    Vi = identidade(n)
    A, autovalores, autovetores = rotGivens(A,Vi,'n')
    n_ = n
    matrix, autovalores, autovetores = rotGivens(A,autovetores,desloc)
    iteracoes = 2
    menor, Aux, A = verificaBeta(matrix)

    while menor == False:
        Aux, autovalores, autovetores = rotGivens(Aux,autovetores,desloc)
        iteracoes+=1
        menor, Aux, A = verificaBeta(Aux)
        if menor == True:
            n_-= 1
            autovetoresAux = tiraColuna(autovetores)
            autovetores = matrizPrincipal(autovetoresAux,autovetores)
            print("autovetores inside 1: \n",autovetores,"\n")
            autovaloresAux = diminuivet(autovalores)
            autovalores = vetorPrincipal(autovaloresAux, autovalores) 
            
            
    Aux, autovaloresAux, autovetoresAux = rotGivens(Aux,autovetoresAux,desloc)
    iteracoes+=1
    menor, Aux, Ai = verificaBeta(Aux)

    while n_ >= 2:
        if menor == False:
            Aux, autovaloresAux, autovetoresAux = rotGivens(Aux,autovetoresAux,desloc)
            iteracoes+=1
            menor, Aux, Ai = verificaBeta(Aux)
            print("autovetoresAux inside: \n",autovetoresAux,"\n")

        if menor == True:
            if n_> 2:
                autovetores = matrizPrincipal(autovetoresAux,autovetores)
                autovalores = vetorPrincipal(autovaloresAux, autovalores)
                #print("autovalores: \n",autovalores,"\n")
                autovetoresAux = tiraColuna(autovetoresAux)
                #print("autovaloresAux: \n",autovaloresAux,"\n")
                #print("autovetores: \n",autovetores,"\n")
                autovaloresAux = diminuivet(autovaloresAux)
                #print("autovetoresAux:\n ",autovetoresAux,"\n")
                autovetores = matrizPrincipal(autovetoresAux,autovetores)
                print("autovetores inside 2: \n",autovetores,"\n")
                n_ -= 1
                #autovetores = matrizPrincipal(autovetoresAux,autovetores)
                #autovalores = vetorPrincipal(autovaloresAux, autovalores)
                A = matrizPrincipal(Ai,A)
                Aux, autovaloresAux, autovetoresAux = rotGivens(Aux,autovetoresAux,desloc)
                iteracoes+=1
                menor, Aux, Ai = verificaBeta(Aux)
            else:
                A = matrizPrincipal(Ai,A)
                autovetores = matrizPrincipal(autovetoresAux,autovetores)
                print("autovetores inside 3: \n",autovetores,"\n")
                autovalores = vetorPrincipal(autovaloresAux, autovalores)
                n_-=1
    print("autovetores: \n",autovetores,"\n")
    return iteracoes, autovalores, autovetores, A

def normalizacao(autovetores): # Normaliza a matriz parametro
    for i in range(len(autovetores)):
        soma = 0
        for j in range(len(autovetores[0])):
            soma += autovetores[i][j]**2
        autovetores[i]=autovetores[i]/(soma)**(1/2)

    return autovetores

# ------------------------------------------------------------------------- #
#                                   Tarefa 
# ------------------------------------------------------------------------- #

# ------------------------------------------------------------------------- #
#                                   item a
# ------------------------------------------------------------------------- #

def tar1(n,desloc): 
    alfak = 2
    betak = -1
    A = identidade(n)
    for i in range(n):
        for j in range(n):
            if j == i:
                A[i][j] = alfak
            if j == i+1 or j == i-1:
                A[i][j] = betak
    iteracoes, autovalores, autovetores, A = QR(A,desloc)
    lambdavec, vj = calculosAnaliticos(n)
    autovetores = ordemAutovetores(autovetores)
    

    print("Número de iterações: ",iteracoes,"\n")
    print("Autovalores encontrados: \n",autovalores,"\n")
    print("Autovalores analíticos: \n",lambdavec,"\n")
    print("Autovetores encontrados: \n",autovetores,"\n")
    print("Autovetores analíticos: \n",vj,"\n")

# ------------------------------------------------------------------------- #
#                                   item b
# ------------------------------------------------------------------------- #
def ki(i): # equação de ki do item b
    ki = 40 + 2*i
    return ki

def kic(i): # equação de ki do item c
    ki = 40 + 2*(-1)**i
    return ki

def matrizki(escolha): # Cria a matriz dos k's do item b
    if escolha == 2:
        m = 5
    elif escolha == 3:
        m = 10
    ks = identidade(m)
    for i in range(len(ks)):
        for j in range(len(ks[0])):
            if j == i:
                ks[i][j] = ki(i+1) + ki(i+2)
            if j == i+1:
                ks[i][j] = -ki(i+1)
            elif j == i-1:
                ks[i][j] = -ki(j+1)
    return ks

def matrizkic(escolha): # Cria a matriz dos k's do item b
    if escolha == 2:
        m = 5
    elif escolha == 3:
        m = 10
    ks = identidade(m)
    for i in range(len(ks)):
        for j in range(len(ks[0])):
            if j == i:
                ks[i][j] = kic(i+1) + kic(i+2)
            if j == i+1:
                ks[i][j] = -kic(i+1)
            elif j == i-1:
                ks[i][j] = -kic(j+1)
    return ks

def matrizAutovalores(autovalores):
    matrizAutovalores = identidade(len(autovalores))
    for i in range(len(autovalores)):
        for j in range(len(autovalores)):
            if i == j:
                matrizAutovalores[i][j] = autovalores[i]
    return matrizAutovalores

def tempo():
    time = np.array([])
    t = 0
    while t <= 10:
        time = np.append(time,t)
        t += 0.025
    return time

def plotar_massas(Xt): # Plota os Valores das massas
    n = len(Xt)

    fig, axs = plt.subplots(n)
    fig.suptitle('Posição das Massas [m] X Tempo [s]')
    t = np.linspace(0,10,int(10/0.025))
    for i in range(n):
        axs[i].plot(t,Xt[i])
        axs[i].legend(str(i + 1),loc='upper right')
    plt.show()

def yt(y0,W,t):
    Yt = np.array([])
    for i in range(len(y0)):
        Yt = np.append(Yt,y0[i]*math.cos(W[i]*t))
    return Yt

def tar2(escolha,n2,dt,tempo):
    num_inter = int(round(tempo/dt))
    A = (1/2) * matrizki(escolha)
    
    iteracoes, autovalores , autovetores,Q = QR(A,'s')
    
    autovalores = matrizAutovalores(autovalores)

    W = np.array([])
    for i in range(len(autovalores)):
        for j in range(len(autovalores[0])):
            if i == j:
                w = math.sqrt(autovalores[i][j])
                W = np.append(W,w)
    if n2 == 1:
        X0 = [-2,-3,-1,-3,-1]
    elif n2 == 2:
        X0 = [1,10,-4,3,-2]
    elif n2 ==3:
        maior = 0
        for i in range(len(autovalores)):
            if W[i]>maior:
                maior = W[i]
                ind = i
        X0 = autovetores[ind].copy()
    Q = autovetores.copy()
    Y0 = Q.T @ X0
    X_posicao = []
    T = []
    for i in range(num_inter):
        t = i*dt
        T.append(t)
        Yt = yt(Y0,W,t)
        Xt = Q @ Yt
        X_posicao.append(Xt)
    posi_massa = np.array(X_posicao).T
    fig ,sub = plt.subplots(nrows=len(posi_massa), figsize=(7,7))
    for i in range(len(posi_massa)):
        sub[i].plot(T,posi_massa[i],label="M{0}".format(i+1))
        sub[i].grid()
        sub[i].legend()
    plt.show()



# ------------------------------------------------------------------------- #
#                                   item c
# ------------------------------------------------------------------------- #
def tar3(escolha,n3,dt,tempo):

    num_inter = int(round(tempo / dt))
    A = (1 / 2) * matrizkic(escolha)

    iteracoes, autovalores, autovetores, Q = QR(A, 's')
    
    autovalores = matrizAutovalores(autovalores)
    W = np.array([])
    for i in range(len(autovalores)):
        for j in range(len(autovalores)):
            if i == j:
                w = math.sqrt(autovalores[i][j])
                W = np.append(W, w)

    
    if n3 == 1:
        X0 = [-2,-3,-1,-3,-1,-2,-3,-1,-3,-1]
    elif n3 == 2:
        X0 = [1,10,-4,3,-2,1,10,-4,3,-2]
    elif n3 ==3:
        maior = 0
        for i in range(len(autovalores)):
            if W[i]>maior:
                maior = W[i]
                ind = i

        X0 = autovetores[ind].copy()
    Q = autovetores.copy()
    Y0 = Q.T @ X0
    X_posicao = []
    T = []
    for i in range(num_inter):
        t = i * dt
        T.append(t)
        Yt = yt(Y0, W, t)
        Xt = Q @ Yt
        X_posicao.append(Xt)
    posi_massa = np.array(X_posicao).T
    fig, sub = plt.subplots(nrows=len(posi_massa), figsize=(7, 7))
    for i in range(len(posi_massa)):
        sub[i].plot(T, posi_massa[i], label="M{0}".format(i + 1))
        sub[i].grid()
        sub[i].legend()
    plt.show()




# ------------------------------------------------------------------------- #
#                                   interface
# ------------------------------------------------------------------------- #
def main():
    print("1) item a")
    print("2) item b")
    print("3) item c")
    escolha = int(input("Escolha qual item do exercício programa: "))
    if escolha == 1: # item a
        desloc = str(input("Com deslocamento (s/n): "))
        n = int(input("Escolha o tamanho n da matriz simétrica desejada: "))
        print("\n")
        print("Os resultados obtidos foram: \n")
        tar1(n,desloc)
    if escolha == 2: # item b
        print("1) X(0) = −2, −3, −1, −3, −1")
        print("2) X(0) = 1, 10, −4, 3, −2")
        print("3) X(0) correspondente ao modo de maior frequência")
        n2 = int(input("Escolha a condição inicial: "))
        tmax = float(input("Escolha o tempo maximo:(sugerido 10s)\n"))
        dt = float(input("Escolha o passo no tempo:(sugerido 0.0025)\n"))
        print("Os resultados obtidos foram: \n")
        tar2(escolha,n2,dt,tmax)
    if escolha == 3:
        print("1) X(0) = −2, −3, −1, −3, −1, −2, −3, −1, −3, −1")
        print("2) X(0) = 1, 10, −4, 3, −2, 1, 10, −4, 3, −2")
        print("3) X(0) correspondente ao modo de maior frequência")
        n3 = int(input("Escolha a condição inicial: "))
        tmax = float(input("Escolha o tempo maximo:(sugerido 10s)\n"))
        dt = float(input("Escolha o passo no tempo:(sugerido 0.0025)\n"))
        print("\n")
        print("Os resultados obtidos foram: \n")
        tar3(escolha,n3,dt,tmax)
main()
