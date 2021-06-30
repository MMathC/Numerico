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
    if abs(alfak)>abs(betak):
        tau = - betak/alfak
        ck = 1/(math.sqrt(1+tau**2))
        sk = ck * tau
    else:
        tau = - alfak/betak
        sk = 1/(math.sqrt(1+tau**2))
        ck = sk * tau
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
    
def QR(A,n, Vi , desloc):
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

def diminui(A, n):
    Anova = identidade(n-1)
    for i in range(n-1):
        for j in range (n-1):
            Anova[i][j] = A[i][j]
    return Anova

def diminuivet(A,n):
    Anova = []
    for i in range(n-1):
        Anova.append(A[i])
    return Anova

def verificaBeta(A,n):
    if abs(A[n-1][n-2]) < 1e-6:
        Aux = diminui(A,n)
        menor = True
        return menor, Aux, A
    else:
        menor = False
        return menor, A, A

def matrizPrincipal(Anova, Amain, n_, n):
    for i in range(len(Anova)):
        for j in range(len(Anova[0])):
            Amain[i][j] = Anova[i][j]
    return Amain

def vetorPrincipal(Anova, Amain, n_):
    for i in range(n_):
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
            vj[j][i] = math.sin(abs(i+1)*abs(j+1)*math.pi/(n+1)) # Auto-vetores analiticos
    return lambdavec, vj
    
    
def erros(A, n, desloc):
    Vi = identidade(n)
    A, autovalores, autovetores = QR(A,n,Vi, 'n')
    n_ = n

    matrix, autovalores, autovetores = QR(A,n_,autovetores,desloc)
    iteracoes = 2
    menor, Aux, A = verificaBeta(matrix,n)


    while menor == False:
        Aux, autovalores, autovetores = QR(Aux,n_,autovetores,desloc)
        iteracoes+=1
        menor, Aux, A = verificaBeta(Aux,n_)
        if menor == True:
            n_-= 1
            autovetoresAux = diminui(autovetores,n)
            autovetores = matrizPrincipal(autovetoresAux,autovetores,n_,n)
            autovaloresAux = diminuivet(autovalores,n)
            autovalores = vetorPrincipal(autovaloresAux, autovalores,n_) 
            
            
    Aux, autovaloresAux, autovetoresAux = QR(Aux,n_,autovetoresAux,desloc)
    iteracoes+=1
    menor, Aux, Ai = verificaBeta(Aux,n_)

    while n_ >= 2:
        if menor == False:
            Aux, autovaloresAux, autovetoresAux = QR(Aux,n_,autovetoresAux,desloc)
            iteracoes+=1
            menor, Aux, Ai = verificaBeta(Aux,n_)

        if menor == True:
            if n_> 2:
                autovetoresAux = diminui(autovetoresAux,n_)
                autovaloresAux = diminuivet(autovaloresAux,n_)
                n_ -= 1
                autovetores = matrizPrincipal(autovetoresAux,autovetores,n_,n)
                autovalores = vetorPrincipal(autovaloresAux, autovalores,n_)
                A = matrizPrincipal(Ai,A,n_,n)
                Aux, autovaloresAux, autovetoresAux = QR(Aux,n_,autovetoresAux,desloc)
                iteracoes+=1
                menor, Aux, Ai = verificaBeta(Aux,n_)
            else:
                autovetores = matrizPrincipal(autovetoresAux,autovetores,n_,n)
                autovalores = vetorPrincipal(autovaloresAux, autovalores,n_)
                n_-=1
                

    return iteracoes, autovalores, autovetores

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
            matrizAux[i][j] = autovetores[i][((j-2*j)-1)]
    print("Número de iterações: ",iteracoes,"\n")
    print("Autovalores encontrados: \n",autovalores,"\n")
    print("Autovalores analíticos: \n",lambdavec,"\n")
    print("Autovetores encontrados: \n",matrizAux,"\n")
    print("Autovetores analíticos: \n",vj,"\n")

# ------------------------------------------------------------------------- #
#                                   item b
# ------------------------------------------------------------------------- #






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
    escolha = int(input("Escolha qual item do exercício programa rodar: "))
    if escolha == 1:
        desloc = str(input("Com deslocamento (s/n): "))
        n = int(input("Escolha o tamanho n da matriz simétrica desejada: "))
        tar1(n,desloc)
main()
