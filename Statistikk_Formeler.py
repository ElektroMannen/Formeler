import math
from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import norm, t
#-----------------------------
#Grunnleggende formeler
#Gjennomsnitt
def gjennomsnitt(dataArray):
    sum = 0
    for i in range (0, len(dataArray)):
        sum += dataArray[i]
    return sum/len(dataArray)

def median(dataSett):
    data = sorted(dataSett)
    if ((len(dataSett) % 2) == 1):
        return data[len(dataSett)//2]
    else: 
        return (data[len(dataSett) // 2 - 1] + data[len(dataSett) // 2]) / 2
    
def variasjonsBredde(dataSett):
    return max(dataSett)- min(dataSett)

def sum_XiminSnittoppi2(dataArray):
    snitt = gjennomsnitt(dataArray)
    sum = 0
    for i in range (0,len(dataArray)):
       sum += (dataArray[i]-snitt)**2
    return sum

def sumXminSnittogYminSnitt(dataSettX, dataSettY):
    sum = 0
    for i in range (0,len(dataSettX)):
        sum += (dataSettX[i] - gjennomsnitt(dataSettX))*(dataSettY[i] - gjennomsnitt(dataSettY))
    return sum
#--------------------
#Emperiske formeler

def emperiskVarians(dataSett):
    return sum_XiminSnittoppi2(dataSett)/(len(dataSett)-1)

def emperiskStandardaavik(dataSett):
    return math.sqrt(emperiskVarians(dataSett))

def emperiskKorrelasjon(dataSettX, dataSettY):
    teller = 0
    nevner = math.sqrt(sum_XiminSnittoppi2(dataSettX))*math.sqrt(sum_XiminSnittoppi2(dataSettY))
    return teller/nevner

#Gir ut y vedi hvis det anngies x
def minsteKvadratsumsRetteLinje(dataSettX,dataSettY, x = None):
    b = (sumXminSnittogYminSnitt(dataSettX,dataSettY))/(sum_XiminSnittoppi2(dataSettX))
    a = gjennomsnitt(dataSettY) - (b*gjennomsnitt(dataSettX))
    if(x != None):
        y = a + b*x
        return y,a,b
    else:
        return None,a,b

#---------------
#Heldelser
def P_AnB(A,B):
    return A * B

def P_AuB(A,B,AnB = None):
    if(AnB != None):
        return A + B - AnB 
    else:
        return A + B - P_AnB(A,B)


def P_AIB(A,B,AnB=None):
    if(AnB != None):
        return (AnB)/B
    else:
        return P_AuB(A,B)/B

def P_BIA(A,B):
    return (P_AIB(A,B)*B)/A
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#-----Alt under her må testes og skrives på nytt--------
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#Stokastiske formeler

#Biomisk fordeling hvor k er medlemer, n er antall forsøk, p er sjansen for at det 
def biomisk_fordelt_stokastisk_variabel(k, n, p):
    biomiske_koeffisienten = 0
    teller = math.factorial(n)
    nevner = math.factorial(k) * math.factorial(n-k)
    biomiske_koeffisienten = teller/(nevner)
    return biomiske_koeffisienten*(p**k)*((1-p)**(n-k))

def stokastisk_variabel_P_mindre_eller_lik(k,n,p):
    sum = 0
    for i in range(0,k+1):
        sum += biomisk_fordelt_stokastisk_variabel(i,n,p)
        print(f"X = {i} er {sum} ")
    return sum

def stokastisk_variabel_P_større_eller_lik(k,n,p):
    sum = 0
    for i in range(0,k+1):
        sum += biomisk_fordelt_stokastisk_variabel(i,n,p)
        print(f"X = {i} er {sum} ")
    return 1-sum

#Stokastisk variaberl hvor P(x < k) og P(k < y)
def stokastisk_variabel_P_imellom(StørstK,MinstK,n,p):
    sumMin = 0
    sumMax = 0
    for i in range(1,MinstK+1):
        sumMin += biomisk_fordelt_stokastisk_variabel(i,n,p)
    for i in range(1,StørstK+1):
        sumMax += biomisk_fordelt_stokastisk_variabel(i,n,p)
    return sumMax-sumMin
    
#--------------------------
#Div formeler

#Varians for en Diskret Sannsynlighetsfordeling
def stokastisk_varians(dataSett):
    var = 0
    Ex = 0
    #Regner ut forventningsverdien
    for i in range (0,len(dataSett)):
        Ex += i * dataSett[i]
    #Finner variansen
    var = sum_XiminSnittoppi2(dataSett)
    return var 

#trekninger uten Trilbakelegging
def trekkning_Uten_Tilbakelegging(antallValg, trekkninger):
    teller = 1
    nevner1 = 1
    nevner2 = 1
    for i in range (1,antallValg+1):
        teller *= (i)
    for i in range (1,trekkninger+1):
        nevner1 *= i
    for i in range (1,antallValg + 1-trekkninger):
        nevner2 *= i
    print(f"Svaret er {teller/(nevner1*nevner2)}")
    return teller/(nevner1*nevner2)

#Possison fordelingens formel
def poissonfordeling(k,lamb,t):
    teller= 0
    nevner = 0
    teller = ((lamb*t)**k) * math.e**(-lamb*t)
    nevner = math.factorial(k)
    return teller/nevner

#Kumolativ posissonfordeling
def kumolativPoisson(k,lamb,t):
    var = 0
    for i in range(0,k+1):
        print(f"X = {i} er {poissonfordeling(i,lamb,t)} ")
        var += poissonfordeling(i,lamb,t)
    return var

#Kumulativ fordeling
def KumulativFordeling(k,lamb):
    return 1-math.e**(-lamb*k)

#Standard avvik
def standardAvik(s):
    return math.sqrt(s)

#Standard normalfordeling
def standardNormalFordeling(z):
    return (1/(math.sqrt(2*math.pi)))*(math.e**((-z**2)/(2)))


#Generell normalfordeling
def GenerellNormalFordeling(x,dataArray):
    førstehalvdel = 1/(standardAvik(emperiskVarians(dataArray))*math.sqrt(2*math.pi))
    andrehalvdel = math.e**(-((x-gjennomsnitt(dataArray))**2)/(2*emperiskVarians(dataArray)))
    return førstehalvdel*andrehalvdel

#Generell normalfordeling
def ManuellGenerellNormalFordeling(x,Savik,gjennomsnitt):

    return (x-gjennomsnitt)/(Savik)

#Korvanse
def Korvanse(dataX,dataY,n):
    teller = 0
    nevner = 0
    for i in range(0,len(dataX)):
        teller += (dataX[i]-gjennomsnitt(dataX))*(dataY[i]-gjennomsnitt(dataY))
    
    nevner = len(dataX)-1
    return teller/nevner

#Konfidensienl intervall
#-----------------------
def konfidensiel_Intervall_uten_standardavvik(datasett,prosent):
    alpha = 1 - prosent/100 
    midt = t.ppf(1-(alpha/2),len(datasett)-1)
    midt1 = 2.262
    bakerst = standardAvik(emperiskVarians(datasett))/(math.sqrt(len(datasett)))
    return gjennomsnitt(datasett) + (midt1 * bakerst), gjennomsnitt(datasett) - (midt1 * bakerst), midt


#----------------------------
# Rapport utskriving av datasett
def rapport(siffer,dataSettX=None, dataSettY = None):
    if (dataSettX != None): 
        print("-------------------------------------------\nRapoort for datasett")
        print("         X verdier")
        print(f"Gjennommsnitt X: {round(gjennomsnitt(dataSettX),siffer)}")
        print(f"Median X: {round(median(dataSettX),siffer)}")
        print(f"Variasjonsbredde X: {round(variasjonsBredde(dataSettX),siffer)}")
        print(f"Emperisk varians X: {round(emperiskVarians(dataSettX),siffer)}")
        print(f"Emperisk standardavvik X: {round(emperiskStandardaavik(dataSettX),siffer)}")
        print(f"VAR(X): ")
        print(f"SD(X): {round(standardAvik(emperiskVarians(dataSettX)),siffer)}")
        print("-------------------------------------------")
    if (dataSettY != None):
        print("         Y verdier")
        print(f"Gjennommsnitt Y: {round(gjennomsnitt(dataSettY),siffer)}")
        print(f"Median X: {round(median(dataSettY),siffer)}")
        print(f"Variasjonsbredde Y: {round(variasjonsBredde(dataSettY),siffer)}")
        print(f"Emperisk varians Y: {round(emperiskVarians(dataSettY),siffer)}")
        print(f"Emperisk standardavvik XY {round(emperiskStandardaavik(dataSettY),siffer)}")
        print(f"VAR(Y): ")
        print(f"SD(Y): {round(standardAvik(emperiskVarians(dataSettY)),siffer)}")
        print("-------------------------------------------")
    if((dataSettX != None) and (dataSettY != None)):
        print("         Felles verdier")
        print(f"Emperisk korrelasjon er: {round(emperiskKorrelasjon(dataSettX,dataSettY),siffer)}")
        print(f"Minste kvadratsums rette linje (y = a+bx): y = {round(minsteKvadratsumsRetteLinje(dataSettX,dataSettY)[1],siffer)} + {round(minsteKvadratsumsRetteLinje(dataSettX,dataSettY)[2],siffer)}x")
        print("-------------------------------------------\n")

    return "Rapport printet"