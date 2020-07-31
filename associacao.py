from trocadores import *
import numpy as np

'''
@author: Igor Bastos
@version: 1.0
'''

class Associa:
    """Recebe um argumento 'tipo', que pode ser 'contra', 'paralelo' ou 'misto'. Tanto 'contra' quanto 'paralelo' se referem à associação em série; 'misto' é a associação combinada em série e paralelo"""
    def __init__(self,tipo):
        self.tipo = tipo
    def P1(self,P1b=1.0,R1=1.0,n=2.0):
        if self.tipo == 'contra':
            """Retorna P1. Exige como inputs R1, P1b e n"""
            A = (1-R1*P1b)**(n)
            B = (1-P1b)**(n)
            return((A-B)/(A-R1*B))
        if self.tipo == 'paralelo':
            """Retorna P1. Exige como inputs R1, P1b e n"""
            return(1/(1+R1)*(1-(1-(1+R1)*P1b)**(n)))
        if self.tipo == 'misto':
            """Retorna P1. Exige como inputs P1b e n"""
            return(1-(1-P1b)**(n))
    def NUT1(self,NUT1b,n=2.0):
        """Retorna NUT1. Exige como inputs NUT1b e n"""
        if self.tipo == 'contra':
            return n*NUT1b
        if self.tipo == 'paralelo':
            return n*NUT1b
        if self.tipo == 'misto':
            return n*NUT1b
    def P1b(self,P1=1.0,R1=1.0,n=2.0):
        if self.tipo == 'contra':
            """Retorna P1b. Exige como inputs R1, P1 e n"""
            Fa = ((1-R1*P1)/(1-P1))**(1/n)
            return ((Fa-1)/(Fa-R1))
        if self.tipo == 'paralelo':
            """Retorna P1b. Exige como inputs R1, P1 e n"""
            A = 1/(1+R1)
            B = (1+R1)*P1
            C = (1-B)**(1/n)
            D = 1-C
            return(A*D)
        if self.tipo == 'misto':
            """Retorna P1b. Exige como inputs P1 e n"""
            return(1-(1-P1)**(1/n))
    def NUT1b(self,NUT1,n=2.0):
        return NUT1/n
    def n(self,P1=1.0,R1=1.0,P1b=1.0):
        if self.tipo == 'contra':
            """Retorna n. Exige como inputs R1, P1 e P1b"""
            A = (1-R1*P1)
            B = (1-P1)
            C = (1-R1*P1b)
            D = (1-P1b)
            return((np.log((A)/(B)))/(np.log((C)/(D))))
        if self.tipo == 'paralelo':
            """Retorna n. Exige como inputs R1, P1 e P1b"""
            A = (1+R1)*P1
            B = (1+R1)*P1b
            return((np.log(1-A))/(np.log(1-B)))
        if self.tipo == 'misto':
            """Retorna n. Exige como inputs P1 e P1b"""
            return((np.log(1+P1))/(1+P1b))
