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
        #print("n: ",n)
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
    
def QR(A, HT,desloc):
    n = len(A)
    Vi = HT
    A, autovalores, autovetores = rotGivens(A,Vi,'n')
    n_ = n
    menor = False
    iteracoes = 1
    while menor == False:
        Aux, autovalores, autovetores = rotGivens(A,autovetores,desloc)
        iteracoes += 1
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

def Hw(wbarra,n):
    I = identidade(len(wbarra))
    Hwi = I - 2*(wbarra@wbarra.T)/normaAoQuadrado(wbarra)
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
            if abs(A[i][j]) < 1e-8:
                A[i][j] = 0
    return A

def TransformacaoHouseholder(A):
    n = len(A)
    iteracao = 0
    Hw1 = Hw(wBarra(a(A,iteracao)),n)
    iteracao+=1
    HT = Hw1
    i = n-1
    H = Hw1@A@Hw1
    while i > 2:
        Hwn = Hw(wBarra(a(H,iteracao)),n)
        H = Hwn@H@Hwn
        HT = HT@Hwn
        i-=1
        iteracao+=1
    H = arrumaZeros(H)
    print("H: \n",np.round(H,1),"\n")
    print("HT: \n",np.round(HT,7),"\n")
    return H, HT

def matrizAutovalores(autovalores):
    ident = identidade(len(autovalores))
    for i in range(len(ident)):
        for j in range(len(ident[0])):
            if ident[i][j] == 1:
                ident[i][j] = autovalores[i]
    return ident

def LerArquivo(nome):
    with open(nome,'r') as arq:
        content = arq.readlines()
        tamanho = int(content[0])
        newList = []
        line = 0
        for i in range(1,tamanho+1):
            newList.append([])
            linhas = content[i].split()
            for j in range(tamanho):
                newList[line].append(int(linhas[j]))
            line+=1
        A = np.array(newList)
    return A

def LerArquivoC(nome):
    with open(nome,'r') as arq:
        content = arq.readlines()
        primeiraLinha = content[0].split()
        segundaLinha = content[1].split()
        nosTotal = int(primeiraLinha[0])
        print("Numeros de nós total: ",nosTotal)
        nosNaoFixos = int(primeiraLinha[1])
        print("Número de nós não fixos: ",nosNaoFixos)
        numBarras = int(primeiraLinha[2])
        print("Número de barras: ",numBarras)
        ro = float(segundaLinha[0])
        print("Ro",ro)
        area = float(segundaLinha[1])
        print("A: ",area)
        E = float(segundaLinha[2])
        print("E: ",E)
        line = 0
        K = MatrizZerada(24)
        for linha in range(2,29):
            linhas = content[linha].split()
            for coluna in range(4):
                if coluna == 0:
                    i = int(linhas[coluna])
                    print("i1: ",i)
                elif coluna == 1:
                    j = int(linhas[coluna])
                    print("j1: ",j)
                elif coluna == 2:
                    angle = float(linhas[coluna])
                    print("angulo: ",deg2rad(angle))
                else:
                    L = float(linhas[coluna])
                    print("L1: ",L)
            k = MatrizDeRigidezIntermediaria(area,E,L,deg2rad(angle))
            print(k)
            K = MatrizDeRigidez(i,j,k,K)
            line+=1
        EscreverArquivo(np.round(K,3))
        return K
            
def deg2rad(angle):
    return angle*math.pi/180

def EscreverArquivo(A):
    nome = str(input("nome do arquivo: "))
    with open(nome,'w') as arq:
        for i in range(len(A)):
            for j in range(len(A[i])):
                if j == len(A[i])-1:
                    arq.write(str(A[i][j]))
                else:
                    arq.write(str(A[i][j])+" ")
            arq.write("\n")
    arq.close()
    
def AutovaloresAnaliticos():
    lambdai = []
    n = 20
    for i in range(1,n+1):
        autoval = 1/(2*(1-(math.cos(math.pi*(2*i - 1)/(2*n + 1)))))
        lambdai.append(autoval)
    return lambdai

def MatrizDeRigidezIntermediaria(area,E,L,angle):
    C = math.cos(angle)
    S = math.sin(angle)
    K = []
    K.append([C**2,C*S,-C**2,-C*S])
    K.append([C*S,S**2,-C*S,-S**2])
    K.append([-C**2,-C*S,C**2,C*S])
    K.append([-C*S,-S**2,C*S,S**2])
    K = np.array(K)
    print("k: \n",K,"\n")
    print("Area: ",area)
    print("E: ",E)
    print("L: ",L)
    print("angle: ",angle)
    return np.array((area*E/L)*K)

def MatrizZerada(n):
    K = identidade(n)
    for i in range(len(K)):
        for j in range(len(K)):
            if i == j:
                K[i][j] = 0
    return K
    
def MatrizDeRigidez(i,j,k,K):
    print("i: ",i," j: ",j)
    if (2*i) < 24 and (2*j) < 24:
        K[2*i-1, 2*i-1]+=k[0][0]
        K[2*i, 2*i-1]+=k[1][0]
        K[2*j-1, 2*i-1]+=k[2][0]
        K[2*j, 2*i-1]+=k[3][0]
        K[2*i-1, 2*i]+=k[0][1]
        K[2*i, 2*i]+=k[1][1]
        K[2*j-1, 2*i]+=k[2][1]
        K[2*j, 2*i]+=k[3][1]
        K[2*i-1, 2*j-1]+=k[0][2]
        K[2*i, 2*j-1]+=k[1][2]
        K[2*j-1, 2*j-1]+=k[2][2]
        K[2*j, 2*j-1]+=k[3][2]
        K[2*i-1, 2*j]+=k[0][3]
        K[2*i, 2*j]+=k[1][3]
        K[2*j-1, 2*j]+=k[2][3]
        K[2*j, 2*j]+=k[3][3]
    return np.array(K)

def ArrumaZerrosVetor(vetor):
    newList = np.array([])
    for i in range(len(vetor)):
        if abs(vetor[i]) < 1e-6:
            newList = np.append(newList,0)
        else:
            newList = np.append(newList,vetor[i])
    return newList

def Verifica(A,autovetores,autovalores):
    newList = np.array([])
    for i in range(len(autovalores)):
        for j in range(len(autovalores)):
            newList = np.append(newList,autovetores[j][i]) 
        print("Av: \n",ArrumaZerrosVetor(A@newList),"\n")
        print("lambdav: \n",autovalores[i]*newList,"\n")
        newList = np.array([])
    
    
# ------------------------------------------------------------------------- #
#                  ]                 Tarefa 
# ------------------------------------------------------------------------- #

def tarefa1(escolha):
    if escolha == 1:
        A = LerArquivo("input-a")
        print("Matriz de entrada: \n",A,"\n")
        H, HT = TransformacaoHouseholder(A)
        iteracoes, autovalores, autovetores, matrizDeAutovalores = QR(H,HT,'s')
        autovalores = matrizAutovalores(autovalores)
        autovetores = ordemAutovetores(autovetores)
        Verifica(A,autovetores,autovalores)
        print("autovetores: \n",np.round(autovetores,7),"\n")
        print("Autovetores normalizados: \n",np.round(normalizacao(autovetores),7),"\n")
        #print("Verificação A*v = lambda*v")
        print("T = H A H.T: \n",autovetores@matrizDeAutovalores@autovetores.T,"\n")
        print("Matriz inicial: \n",A,"\n")
        #print("Av: \n",arrumaZeros(A@autovetores),"\n")
        #autovalores = np.array([[7,0,0,0],[0,2,0,0],[0,0,-1,0],[0,0,0,-2]])
        #print("lambdav: \n",arrumaZeros((autovalores@autovetores).T),"\n")
    elif escolha == 2:
        A = LerArquivo("input-b")
        print("Matriz de entrada: \n",A,"\n")
        H, HT = TransformacaoHouseholder(A)
        iteracoes, autovalores, autovetores, matrizDeAutovalores = QR(H,HT,'s')
        autovalores = matrizAutovalores(autovalores)
        print("Autovalores analíticos: \n",AutovaloresAnaliticos(),"\n")
        autovetores = ordemAutovetores(autovetores)
        #EscreverArquivo(np.round(autovetores,7)) # Caso queira gerar um arquivo com a matriz de autovetores, para melhor exibição e analise dos resultados, retire o primeiro hashtag dessa linha
        print("Verificação A*v = lambda*v")
        print("T = H A H.T: \n",autovetores@matrizDeAutovalores@autovetores.T,"\n")
        #EscreverArquivo(np.round(H,2)) # Caso queira gerar um arquivo com a matriz H (tridiagonal simérica), para melhor exibição e analise dos resultados, retire o primeiro hashtag dessa linha
        #EscreverArquivo(np.round(T,0)) # Caso queira gerar um arquivo com a matriz T, para melhor exibição e analise dos resultados, retire o primeiro hashtag dessa linha
        
        
def main():
    print("1) Testes")
    print("2) Aplicações: Treliças Planas")
    escolha = int(input("Escolha qual item do exercício programa: "))
    if escolha == 1:
        print("1) Teste a")
        print("2) Teste b")
        testes = int(input("Escolha qual item dos testes: "))
        print("\nresultados obtidos foram: \n")
        tarefa1(testes)
    elif escolha == 2:
        print("\nresultados obtidos foram: \n")
        LerArquivoC('input-c')                
main()
