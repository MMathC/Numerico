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

def sgn(d):
    if d >= 0:
        return 1
    else:
        return -1
    
def heuristicaWilkinson(alfan_1, alfan, betan_1, desloc):
    if desloc == 's':
        dk = (alfan_1 - alfan)/2
        mik = alfan + dk - sgn(dk) * math.sqrt(dk**2 + betan_1**2)
        return mik
    return 0

def deslocamento(mik,n):
    return mik * identidade(n)
    
def QR(A, Vi, desloc):
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
    Ak = (R @ Q_ts) + deslocamento(mik * identidade(n),n)
    autovalores = autovalor(Ak)
    
    return Ak, autovalores, autovetores

def diminui(A):
    Anova = identidade(len(A)-1)
    for i in range(len(A)-1):
        for j in range (len(A)-1):
            Anova[i][j] = A[i][j]
    return Anova

def diminuivet(A):
    Anova = np.array([])
    for i in range(len(A)-1):
        Anova = np.append(Anova,A[i])
    return Anova

def verificaBeta(A):
    if abs(A[len(A)-1][len(A)-2]) < 1e-6:
        Aux = diminui(A)
        menor = True
        return menor, Aux, A
    else:
        menor = False
        return menor, A, A

def matrizPrincipal(Anova, Amain):
    for i in range(len(Anova)):
        for j in range(len(Anova[0])):
            Amain[i][j] = Anova[i][j]
    return Amain

def vetorPrincipal(Anova, Amain):
    for i in range(len(Anova)):
        Amain[i] = Anova[i]
    return Amain

def calculosAnaliticos(n):
    lambdavec = []
    vj = identidade(n)
    for j in range(-n,0):
        lambdaj = 2*(1-math.cos((-j)*math.pi/(n+1))) # Auto-valores analiticos
        lambdavec.append(lambdaj)
    for j in range (n):
        for i in range (n):
            vj[i][j] = math.sin(abs(i+1)*abs(j+1)*math.pi/(n+1)) # Auto-vetores analiticos
    return lambdavec, vj
    
    
def erros(A, n, desloc):
    Vi = identidade(n)
    A, autovalores, autovetores = QR(A,Vi,'n')
    n_ = n
    print("0 autovetores: \n",autovetores,"\n")
    matrix, autovalores, autovetores = QR(A,autovetores,desloc)
    iteracoes = 2
    menor, Aux, A = verificaBeta(matrix)


    while menor == False:
        print("iteração: ",iteracoes)
        print("1 autovetores: \n",autovetores,"\n")
        Aux, autovalores, autovetores = QR(Aux,autovetores,desloc)
        print("1 autovetores: \n",autovetores,"\n")
        iteracoes+=1
        menor, Aux, A = verificaBeta(Aux)
        if menor == True:
            n_-= 1
            print("iteração: ",iteracoes)
            autovetoresAux = diminui(autovetores)
            print("2 autovetoresAux: \n",autovetoresAux,"\n")
            autovetores = matrizPrincipal(autovetoresAux,autovetores)
            print("2 autovetores: \n",autovetores,"\n")
            autovaloresAux = diminuivet(autovalores)
            autovalores = vetorPrincipal(autovaloresAux, autovalores) 
            
            
    Aux, autovaloresAux, autovetoresAux = QR(Aux,autovetoresAux,desloc)
    iteracoes+=1
    menor, Aux, Ai = verificaBeta(Aux)

    while n_ >= 2:
        if menor == False:
            print("iteração: ",iteracoes)
            print("3 autovetoresAux: \n",autovetoresAux,"\n")
            Aux, autovaloresAux, autovetoresAux = QR(Aux,autovetoresAux,desloc)
            iteracoes+=1
            menor, Aux, Ai = verificaBeta(Aux)
            print("3 autovetoresAux: \n",autovetoresAux,"\n")

        if menor == True:
            if n_> 2:
                print("iteração: ",iteracoes)
                autovetoresAux = diminui(autovetoresAux)
                print("4 autovetoresAux: \n",autovetoresAux,"\n")
                autovaloresAux = diminuivet(autovaloresAux)
                n_ -= 1
                autovetores = matrizPrincipal(autovetoresAux,autovetores)
                print("4 autovetores: \n",autovetores,"\n")
                autovalores = vetorPrincipal(autovaloresAux, autovalores)
                A = matrizPrincipal(Ai,A)
                Aux, autovaloresAux, autovetoresAux = QR(Aux,autovetoresAux,desloc)
                iteracoes+=1
                print("4 autovetoresAux: \n",autovetoresAux,"\n")
                menor, Aux, Ai = verificaBeta(Aux)
                print("menor: ",menor,"\n")
            else:
                print("iteração: ",iteracoes)
                autovetores = matrizPrincipal(autovetoresAux,autovetores)
                print("5 autovetores: \n",autovetores,"\n")
                autovalores = vetorPrincipal(autovaloresAux, autovalores)
                n_-=1
                
    return iteracoes, autovalores, autovetores

def normalizacao(autovetores):
    normalizado = 0
    normas = np.array([])
    for j in range(len(autovetores)):
        for i in range(len(autovetores[0])):
            normalizado += (autovetores[i][j])**2
        norma = math.sqrt(normalizado)
        normas = np.append(normas, norma)
    for j in range(len(autovetores)):
        for i in range(len(autovetores[0])):
            autovetores[j][i] = autovetores[j][i]/normas[j]
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
    iteracoes, autovalores, autovetores = erros(A,n,desloc)
    lambdavec, vj = calculosAnaliticos(n)
    matrizAux = identidade(n)
    for i in range(n):
        for j in range(-n,0):
            matrizAux[i][j] = autovetores[i][((-j)-1)]
    
    

    print("Número de iterações: ",iteracoes,"\n")
    print("Autovalores encontrados: \n",autovalores,"\n")
    print("Autovalores analíticos: \n",lambdavec,"\n")
    print("Autovetores encontrados: \n",matrizAux,"\n")
    print("Autovetores analíticos: \n",vj,"\n")

# ------------------------------------------------------------------------- #
#                                   item b
# ------------------------------------------------------------------------- #
def ki(i):
    ki = 40 + 2*i
    return ki

def matrizki():
    ks = identidade(6)
    for i in range(len(ks)):
        for j in range(len(ks[0])):
            if j == i:
                ks[i][j] = ki(i+1) + ki(i+2)
            if j == i+1:
                ks[i][j] = -ki(i+1)
            elif j == i-1:
                ks[i][j] = -ki(j+1)
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

#Funcao que Plota os Valores das massas
def plotar_massas(Xt):
    n = len(Xt)

    fig, axs = plt.subplots(n)
    fig.suptitle('Posição das Massas [m] X Tempo [s]')
    t = np.linspace(0,10,int(10/0.025))
    for i in range(n):
        axs[i].plot(t,Xt[i])
        axs[i].legend(str(i + 1),loc='upper right')
    plt.show()

def tar2(m):
    A = (1/m) * matrizki()
    iteracoes, autovalores, X = erros(A,len(A),'s')
    #iteracoes, autovalores, X = erros(X,-3,'s')
    #iteracoes, autovalores, X = erros(X,-1,'s')
    #iteracoes, autovalores, X = erros(X,-3,'s')
    #iteracoes, autovalores, X = erros(X,-1,'s')
    
    print(autovalores)
    print(X)
    a = matrizAutovalores(autovalores)
    print(a)
    n = len(autovalores)
    print(tempo())




# ------------------------------------------------------------------------- #
#                                   item c
# ------------------------------------------------------------------------- #






# ------------------------------------------------------------------------- #
#                                   interface
# ------------------------------------------------------------------------- #
def main():
    print("1 - item a")
    print("2 - item b")
    print("3 - item c")
    escolha = int(input("Escolha qual item do exercício programa: "))
    if escolha == 1:
        desloc = str(input("Com deslocamento (s/n): "))
        n = int(input("Escolha o tamanho n da matriz simétrica desejada: "))
        print("\n")
        print("Os resultados obtidos foram: \n")
        tar1(n,desloc)
    if escolha == 2:
        m = int(input("Massa: "))
        
        print("Os resultados obtidos foram: \n")
        tar2(m)
        
main()
