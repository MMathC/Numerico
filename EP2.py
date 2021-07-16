# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 16:19:35 2021

@author: macas
"""

import numpy as np
import math

# ------------------------------------------------------------------------- #
#                                   EP1 
# ------------------------------------------------------------------------- #
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
                autovalores = vetorPrincipal(autovaloresAux, autovalores)
                n_-=1
                
    autovetores = ordemAutovetores(autovetores)
    
    print("Número de iterações: ",iteracoes,"\n")
    print("Autovalores encontrados: \n",autovalores,"\n")
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
#                                   EP2 
# ------------------------------------------------------------------------- #

def normaAoQuadrado(a):
    c = 0
    for i in range(len(a)):
        c += a[i]**2 
    return c

def norma(a):
    return math.sqrt(normaAoQuadrado(a))

def e(n): # retorna vetor e
    e = np.array([])
    e = np.append(e,1)
    for i in range(n-1):
        e = np.append(e,0)
    return e

def wBarra(an):
    delta = an[0]/abs(an[0])
    wbarraT = an + delta*norma(an)*e(len(an))
    wbarra = []
    for i in range(len(an)):
        wbarra.append([wbarraT[i]])
    return np.array(wbarra)

def produtoVetorial(a,b):
    produto = a@b
    return produto

def Hw(wbarra,n):
    I = identidade(len(wbarra))
    Hwi = I - 2*produtoVetorial(wbarra,wbarra.T)/normaAoQuadrado(wbarra)
    Hw = identidade(n)
    diferenca = n - len(wbarra)
    for i in range(n):
        for j in range(n):
            if i == j and i < diferenca:
                Hw[i][j] = 1
            elif i < diferenca or j < diferenca and i !=j:
                Hw[i][j] = 0
            else:
                Hw[i][j] = Hwi[i-diferenca][j-diferenca]
    return Hw

def a(A,k):
    a = np.array([])
    for i in range(len(A)):
        for j in range(len(A[0])):
            if j == k and i > k:
                a = np.append(a,A[i][j])
    return a

def arrumaZeros(A):
    for i in range(len(A)):
        for j in range(len(A[0])):
            if abs(A[i][j]) < 1e-10:
                A[i][j] = 0
    return A

def TransformacaoHouseholder(A):
    n = len(A)
    Hw1 = Hw(wBarra(a(A,0)),n)
    print("Hw1: \n",Hw1,"\n")
    HT = Hw1
    i = n-1
    while i > 2:
        H = Hw1@A@Hw1
        print("H: \n",arrumaZeros(H),"\n")
        print("a: \n",a(H,1),"\n")
        Hwn = Hw(wBarra(a(H,1)),n)
        print(Hwn,"\n")
        H = Hwn@H@Hwn
        HT = HT@Hwn
        i-=1
    H = arrumaZeros(H)
    print("H: \n",H,"\n")
    print("HT: \n",HT,"\n")
    return H, HT
    
# ------------------------------------------------------------------------- #
#                                   Tarefa 
# ------------------------------------------------------------------------- #

def tar1(): 
    A = np.array([[2,4,1,1],[4,2,1,1],[1,1,1,2],[1,1,2,1]])
    H, HT = TransformacaoHouseholder(A)
    iteracoes, autovalores, autovetores, A = QR(A,'s')


