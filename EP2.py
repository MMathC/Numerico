# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 16:19:35 2021

@author: macas
"""

import numpy as np
import math

A = np.array([-4.3166, 1, 3])
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
    print(Hw1,"\n")
    H = Hw1@A@Hw1
    print("H: \n",H,"\n")
    print("a: \n",a(H,1),"\n")
    Hwn = Hw(wBarra(a(H,1)),n)
    print(Hwn,"\n")
    H = Hwn@H@Hwn
    print(arrumaZeros(H),"\n")
    print(Hw1@Hwn,"\n")

    
def tiraColuna(A):
    Anova = []
    for i in range(len(A)):
        Anova.append([])
        for j in range(len(A[0])-1):
            Anova[i].append(A[i][j])
    return np.array(Anova)
