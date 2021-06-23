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



def autovetores(n): # retorna a matriz identidade com ordem n x n 
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