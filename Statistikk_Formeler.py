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

def sum_XiminSnittoppi2(dataArray):
    snitt = gjennomsnitt(dataArray)
    sum = 0
    for i in range (0,len(dataArray)):
       sum += (dataArray[i]-snitt)**2
    return sum
#--------------------
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

#Emperisk varians
def emperiskVarians(dataSett):
    return sum_XiminSnittoppi2(dataSett)/(len(dataSett)-1)

#Uavhengig
def uavhengig_Hendelse(Pa,Pb):
    return Pa*Pb

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


#---------------------------
# Rapport utskriving av datasett
def rapport(dataSett):
    print("-------------------------------------------\nRapoort for datasett")
    print(f"\nGjennommsnitt: {round(gjennomsnitt(dataSett),3)}")
    print(f"Emperisk varians: {round(emperiskVarians(dataSett),3)}")
    print(f"VAR(X): ")
    print(f"SD(X): {round(standardAvik(emperiskVarians(dataSett)),3)}")
    print("-------------------------------------------\n")

    return "Rapport printet"