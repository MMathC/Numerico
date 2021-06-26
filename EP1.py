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
        #print(v)
        V.append(v)
        v = []
        line+= 1
    return np.array(V)

#def escalonador(A, b):
    

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
    return A, b

'''
A = np.array([[2,1,0,0],[1,2,1,0],[0,1,2,1],[0,0,1,2]])
R = np.array([[5/math.sqrt(5), 4/math.sqrt(5), 1/math.sqrt(5), 0],[0, 3/math.sqrt(5), 2/math.sqrt(5), 0],[0,1,2,1],[0,0,1,2]])
Q1 = R @ np.linalg.inv(A)
print("Matriz Q: \n",np.round(Q1,5))'''

        
'''
A = np.array([[2,1,0,0],[1,2,1,0],[0,1,2,1],[0,0,1,2]])
print(A)
Q1, R1 = np.linalg.qr(A)
B = Q1.T @ A @ Q1
print(np.around(Q1@A,5),"\n")
print(np.around(R1,5),"\n")
print(np.around(B,5),"\n")
for i in range(20):
    Q1, R1 = np.linalg.qr(B)
    B = Q1.T @ B @ Q1
    
print(np.around(Q1@B,3),"\n")'''
#print(Q1@A)
    

A = np.array([[4,3,0],[3,4,3],[0,3,4]])
Rn = np.array([[5, 24/5, 9/5],[0, 7/5, 12/5],[0,3,4]])
Q1 = Rn @ np.linalg.inv(A)
R = [[5, 24/5, 9/5],[0, 54.8/math.sqrt(247),76.8/math.sqrt(247)],[0,0,-8/math.sqrt(247)]]
Q2 = R @ np.linalg.inv(Q1 @ A)

print("Matriz Q1: \n",np.round(Q1,5),"\n")
print("Matriz Rn: \n",np.round(Rn,5),"\n")
print("Matriz Q2: \n",np.round(Q2,5),"\n")
print("Matriz R: \n",np.round(Q2@Q1@A,5),"\n") #errado

print("Matriz A1: \n",np.round(Q2@Q1@A@Q1.T@Q2.T,6),"\n")
def tGivens (A, i, j, k):
    alfak = A[i][k]
    betak = A[j][k]
    print("alfak: ",alfak,"\n")
    print("betak: ",betak,"\n")
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
    print("i: ",i)
    print("j: ",j)
    print("k: ",k,"\n")
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
    print("i: ",i-1)
    print("j: ",j)
    print("k: ",k,"\n")
    return A, i-1, k
    
def QR(A,n):
    Qr = identidade(n)
    ck, sk = tGivens(A,0,1,0)
    auxi = 0
    auxj = 0
    Q1 = identidade(n)
    Q1, auxi, auxj = ajuste(Q1, auxi, auxi+1, auxj, ck, sk)
    R = Q1@A
    print("Q1: \n", np.around(Q1,6), "\n")
    print("R: \n", np.around(R,6), "\n")
    Q_ts = Q1.T
    print("A: \n", np.around(A,6), "\n")
    for r in range(n-2):
        ck, sk = tGivens(R, auxi, auxi+1, auxj)
        Qr, auxi, auxj = ajuste(Qr, auxi, auxi+1, auxj, ck, sk)
        print("Q2: \n",np.around(Q2,6),"\n")
        Q_ts = Q_ts @ Qr.T
        print("Qs: \n",np.around(Q_ts,6),"\n")
        R = Qr@R
        Qr = identidade(n)
        print("R: \n",np.around(R,6),"\n")
    Ak = R @ Q_ts
    print("Ak: \n",np.around(Ak,6))
    print(autovalor(Ak))
    autovalores = autovalor(Ak)
    return Ak, autovalores
# ------------------------------------------------------------------------- #
#                                   Tarefa 
# ------------------------------------------------------------------------- #

# ------------------------------------------------------------------------- #
#                                   item a
# ------------------------------------------------------------------------- #

def tar1(n):
    #vj = (math.sin(j*math.pi/(n+1))) # Auto-vetores
    alfak = 2
    betak = -1
    A = identidade(n)
    for i in range(n):
        for j in range(n):
            if j == i:
                A[i][j] = alfak
            if j == i+1 or j == i-1:
                A[i][j] = betak
    return A
n = 4
a, autovalores = QR(tar1(n),n)
matrix, autovalores = QR(a,n)
lambdavec = []
for j in range(-n,0):
    lambdaj = 2*(1-math.cos((j)*math.pi/(n+1))) # Auto-valores
    lambdavec.append(lambdaj)
print("auto-valores: \n",lambdavec)


erros = []
iteracoes = 2

i = 0
while i < n:
    logic = abs(lambdavec[i]-autovalores[i])
    print("logic: ",logic)
    if logic > 1e-6:
        i = 0
        matrix, autovalores = QR(matrix,n)
        iteracoes+=1
    else:
        erros.append(logic)
        i+=1
    
print("erros: \n",erros)
print("iterações: ",iteracoes)
    
'''
    for j in range(n):
        lambdaj = 2*(1-math.cos(j*math.pi/(n+1))) # Auto-valores
    print(A)
A = np.array([[2,-1,0,0],[-1,2,-1,0],[0,-1,2,-1],[0,0,-1,2]]) # Matriz tridiagonal simetrica
'''


'''   
print("autovalor: ",autovalor([[1,10,3],[0,99,22],[0,0,8]]))
print("det: ",det(autovalor([[1,10,3],[0,99,22],[0,0,8]])))

print("autovalor: ",autovalor([[1,10,3,4],[0,99,22,32],[0,0,8,29],[0,0,0,73]]))
print("det: ",det(autovalor([[1,10,3,4],[0,99,22,32],[0,0,8,29],[0,0,0,73]])))


A = np.array([[4,3,0],[3,4,3],[0,3,4]])
Q1 , R1 = np.linalg.qr(A)
print(np.around(Q1))
print(np.around(R1))
B = Q1 @ R1
print(np.around(B, 5))
print(np.around(np.linalg.inv(Q1)))'''
