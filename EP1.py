# -*- coding: utf-8 -*-
"""
Created on Sat Jun 12 19:36:58 2021

@author: matheus
"""
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
                v.append(1)
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
        
    for rG in range(0, col, 1):
        aux1 = c * A[i][rG] + c * A[j][rG]
        A[i][rG] = aux1
    aux2 = b[i]
    b[i] = c * b[i] - s * b[i]
    b[j] = s * aux2 - c * b[j]
    return A, b

A = np.array([[2,1,0,0],[1,2,1,0],[0,1,2,1],[0,0,1,2]])
R = np.array([[5/math.sqrt(5), 4/math.sqrt(5), 1/math.sqrt(5), 0],[0, 3/math.sqrt(5), 2/math.sqrt(5), 0],[0,1,2,1],[0,0,1,2]])
Q1 = R @ np.linalg.inv(A)
print("Matriz Q: \n",np.round(Q1,5))


c1 = 2/ math.sqrt(5)
s1 = -1/ math.sqrt(5)
Q = [[c1,-s1,0,0],[s1,c1,0,0],[0,0,1,0],[0,0,0,1]]

R1 = Q@A
print("Matriz R1: \n", np.around(R1,5))
'''for i in range(np.shape(A[0])):
    for j in range(np.shape(A[1])):
        b[i] = c1 * b[i] - s1 * b[i]
        b[j] = s1 * aux2 - c1 * b[j]'''
        


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
    
print(np.around(Q1@B,3),"\n")
#print(Q1@A)
    

A = np.array([[4,3,0],[3,4,3],[0,3,4]])
R = np.array([[5, 24/5, 9/5],[0, 7/5, 12/5],[0,3,4]])
Q1 = R @ np.linalg.inv(A)
print("Matriz Q: \n",np.round(Q1,5))

# ------------------------------------------------------------------------- #
#                                   Tarefa 
# ------------------------------------------------------------------------- #

# ------------------------------------------------------------------------- #
#                                   item a
# ------------------------------------------------------------------------- #

def tar1(n):
    #lambdaj = 2*(1-math.cos(j*math.pi/(n+1))) # Auto-valores
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

    print(A)
A = np.array([[2,-1,0,0],[-1,2,-1,0],[0,-1,2,-1],[0,0,-1,2]]) # Matriz tridiagonal simetrica



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
