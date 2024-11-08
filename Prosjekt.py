import Formeler as fm

Diameter = [24.3,
22.1,
22.2,
22.5,
23.3,
23.2,
24,
23.5,
23.9]

#|-----------------|
#|Oppgave 1: Del 1 |
#|-----------------|

#Ta målinger og finn gjennomsnitt og standardavvik
print(f"\n|-----------------|\n|Oppagve 1: del 1 |\n|-----------------|")
print(f"Gjennomsnitt: {round(fm.gjennomsnitt(Diameter),4)}")
print(f"Emperisk standardavvik: {round(fm.standardAvik(fm.sampleVarians(Diameter)),4)}")

#|-----------------|
#|Oppgave 1: Del 2 |
#|-----------------|
print(f"\n|-----------------|\n|Oppagve 1: del 2 |\n|-----------------|")
print(f"Standard Usikkerhet: {round(fm.sampleVarians(Diameter),4)}")
print(f"Forventningsverdi 1: {round(fm.konfidensiel_Intervall_uten_standardavvik(Diameter,95)[0],4)}")
print(f"Forventningsverdi 2: {round(fm.konfidensiel_Intervall_uten_standardavvik(Diameter,95)[1],4)}")
print(f"Med en t_a/2: {round(fm.konfidensiel_Intervall_uten_standardavvik(Diameter,95)[2],4)} 95% burde være {2.262}")

#|-----------------|
#|Oppgave 1: Del 3 |
#|-----------------|
print(f"\n|-----------------|\n|Oppagve 1: del 3 |\n|-----------------|")