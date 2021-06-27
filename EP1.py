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

#def escalonador(A, b):
    
'''
def givens(i, j, col, A, k, b):
    # i: linha superior;
    # j: linha inferior, que contém o elemento que deseja zerar;
    # A: matriz tridiagonal simétrica
    # k: k-ésima coluna da matriz A, onde será aplicada a rotação de Givens;
    if abs(A[i][k]) > abs(A[j][k]):
        tau = -A[j][k] / A[i][k]
        c = 1 / ((1 + tau**2)**0.5)
        s = c * tau
    else:
        tau = -A[i][k] / A[j][k]
        s = 1 / ((1 + tau**2)**0.5)
        c = s * tau
        
    for r in range(0, col, 1):
        aux1 = c * A[i][r] + c * A[j][r]
        A[i][r] = aux1
    aux2 = b[i]
    b[i] = c * b[i] - s * b[i]
    b[j] = s * aux2 - c * b[j]
    return A, b'''
    
A = np.array([[4,3,0],[3,4,3],[0,3,4]])
Rn = np.array([[5, 24/5, 9/5],[0, 7/5, 12/5],[0,3,4]])
Q1 = Rn @ np.linalg.inv(A)
R = [[5, 24/5, 9/5],[0, 54.8/math.sqrt(247),76.8/math.sqrt(247)],[0,0,-8/math.sqrt(247)]]
Q2 = R @ np.linalg.inv(Q1 @ A)

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
    
def heuristicaWilkinson(alfakanterior, alfak, betakanterior):
    dk = (alfakanterior - alfak)/2
    
    mik = alfak - sgn(dk) * math.sqrt(dk**2 - betakanterior**2)
    
    return mik
    
def QR(A,n, Vi):
    Qr = identidade(n)
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
    Ak = R @ Q_ts
    autovalores = autovalor(Ak)
    return Ak, autovalores, autovetores

def erros(A, n):
    Vi = identidade(n)
    a, autovalores, autovetores = QR(A,n,Vi)
    print("autovetores: \n",autovetores,"\n")
    matrix, autovalores, autovetores = QR(a,n,autovetores)
    print("autovetores: \n",autovetores,"\n")
    lambdavec = []
    vj = identidade(n)
    for j in range(-n,0):
        lambdaj = 2*(1-math.cos((-j)*math.pi/(n+1))) # Auto-valores
        lambdavec.append(lambdaj)
        
    for j in range (-n,0):
        for i in range (-n,0):
            vj[j][i] = math.sin(abs(i)*abs(j)*math.pi/(n+1)) # Auto-vetores
    print("vj: \n",vj,"\n")
    
    erros = []
    iteracoes = 2
    i = 0
    j = 0
    while i < n:
        logic = abs(lambdavec[i]-autovalores[i])
        if logic > 1e-6:
            i = 0
            matrix, autovalores, autovetores = QR(matrix,n, autovetores)
            iteracoes+=1
        else:
            erros.append(logic)
            i+=1
            
    '''
    j = 0
    for i in range(n):
        while j < n:
            logic = abs(vj[i][j]-autovetores[i][j])
            if logic > 1e-6:
                j = 0
                matrix, autovalores, autovetores = QR(matrix,n)
                iteracoes+=1
            else:
                j+=1'''
            
            
    return iteracoes, autovalores, autovetores, lambdavec


    
    
    
# ------------------------------------------------------------------------- #
#                                   Tarefa 
# ------------------------------------------------------------------------- #

# ------------------------------------------------------------------------- #
#                                   item a
# ------------------------------------------------------------------------- #

def tar1(n):
    alfak = 2
    betak = -1
    A = identidade(n)
    for i in range(n):
        for j in range(n):
            if j == i:
                A[i][j] = alfak
            if j == i+1 or j == i-1:
                A[i][j] = betak
    iteracoes, autovalores, autovetores, lambdavec = erros(A,n)
    autovetores[:] = autovetores[3::-1]
    print("Número de iterações: ",iteracoes,"\n")
    print("Autovalores encontrados: \n",autovalores,"\n")
    print("Autovalores analíticos: \n",lambdavec,"\n")
    print("Autovetores encontrados: \n",autovetores,"\n")


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
        k = int(input("Escolha o valor de k: "))
        n = int(input("Escolha o tamanho n da matriz simétrica desejada: "))
        tar1(n)
main()
