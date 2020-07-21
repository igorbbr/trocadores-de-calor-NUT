import numpy as np
from scipy.optimize import fsolve
import math

'''
@author: Igor Bastos
@version: 1.1
'''

def dados(a,b):
    '''Retorna uma array 2x1 na forma [a,b], sendo a = P1 ou NUT1 (a depender do caso), e b = R1.'''
    return np.array([[a],[b]])

def coth(x):
    return(math.cosh(x)/math.sinh(x))

excep1 = '{} não é uma variável possível de ser calculada pelo método eps-NUT'
excep2 = '{} não é uma variável possível de ser calculada pelo método P-NUT'

class Trocador:
    '''Cria um objeto Trocador, que recebe um parâmetro tipo = "eps", "nut", "p" ou "pnut", a depender do que o usuário queira calcular.'''
    def __init__(self, tipo):
        self.tipo = tipo
    def TCparalelo(self, Y):
        '''eps-NUT, correntes paralelas'''
        Cr = Y[1]
        if self.tipo == 'nut':
            eps = Y[0]
            return -math.log(1-eps*(1+Cr))/(1+Cr)[0]
        if self.tipo == 'eps':
            NUT = Y[0]
            return (1-math.exp(-NUT*(1+Cr)))/(1+Cr)[0]
        else:
            raise Exception(excep1.format(self.tipo))
    def TCcontra(self,Y):
        '''eps-NUT, correntes contrárias'''
        Cr = Y[1]
        if self.tipo == 'nut':
            eps = Y[0]
            if Cr < 1:
                return ((1)/(1-Cr))*(math.log((1-Cr*eps)/(1-eps)))[0]
            if Cr == 1:
                return (eps/(1-eps))[0]
            else:
                raise Exception('Cr não pode ser maior do que 1.0. Valor atual = {}'.format(Cr))
        if self.tipo == 'eps':
            NUT = Y[0]
            if Cr < 1:
                return ((1-math.exp(-NUT*(1-Cr)))/(1-Cr*math.exp(-NUT*(1-Cr))))[0]
            if Cr == 1:
                return (NUT/(1+NUT))[0]
            else:
                raise Exception('Cr não pode ser maior do que 1.0. Valor atual = {}'.format(Cr))
        else:
            raise Exception(excep1.format(self.tipo))
    def TC1051(self,Y):
        '''P-NUT, Correntes contrárias'''
        R1 = Y[1]
        if self.tipo == 'pnut':
            P1 = Y[0]
            return 1/(1-R1)*math.log((1-R1*P1)/(1-P1))[0]
        if self.tipo == 'p':
            NUT1 = Y[0]
            return((1-math.exp(-NUT1*(1-R1)))/(1-R1*math.exp(-NUT1*(1-R1))))[0]
        else:
            raise Exception(excep2.format(self.tipo))
    def TC1052(self,Y):
        '''P-NUT, Correntes paralelas'''
        R1 = Y[1]
        if self.tipo == 'pnut':
            P1 = Y[0]
            return ((1)/(1+R1)*math.log((1)/(1-P1*(1+R1))))[0]
        if self.tipo == 'p':
            NUT1 = Y[0]
            return((1-math.exp(-NUT1*(1+R1)))/(1+R1))[0]
        else:
            raise Exception(excep2.format(self.tipo))
    def TEMAE2misturado(self,Y):
        '''P-NUT, TEMA E, 2 passes, misturado'''
        R1 = Y[1]
        E = (1+R1**(2))**(1/2)
        if self.tipo == 'pnut':
            P1 = Y[0]
            return(1/E * math.log((2-P1*(1+R1-E))/(2-P1*(1+R1+E))))[0]
        if self.tipo == 'p':
            NUT1 = Y[0]
            return((2)/(1+R1+E*coth(E*NUT1/2)))[0]
        else:
            raise Exception(excep2.format(self.tipo))
    def TEMAE2divsim(self,Y,chute_inicial=1):
        '''P-NUT, TEMA E, 2 passes, dividido simetricamente. Permite especificar um valor para o chute inicial na determinação de PNUT1 = f(P1, R1). Padrão = 1'''
        R1 = Y[1]
        def TEMAE2divsim_P1(NUT1):
            E = np.exp(NUT1)
            B = np.exp(-NUT1*R1/2)
            return(R1**(-1)*(1-((2-R1)*(2*E+R1*B))/((2+R1)*(2*E-R1/B))))
        if self.tipo == 'p':
            NUT1 = Y[0]
            return TEMAE2divsim_P1(NUT1)[0]
        if self.tipo == 'pnut':
            P1 = Y[0]
            mod = lambda NUT: TEMAE2divsim_P1(NUT) - P1
            return fsolve(mod,chute_inicial)[0]
        else:
            raise Exception(excep2.format(self.tipo))
    def TEMAE3(self,Y,chute_inicial=1):
        '''P-NUT, TEMA E, 3 passes. Permite especificar um valor para o chute inicial na determinação de PNUT1 = f(P1, R1). Padrão = 1'''
        R1 = Y[1]
        def TEMAE3_P1(NUT1):
            def lamma(i):
                if i == 1:
                    return(-(3/2)+(9/4+R1*(R1-1))**(1/2))
                if i == 2:
                    return(-(3/2)-(9/4+R1*(R1-1))**(1/2))
                if i == 3:
                    return(R1)
            def X(i):
                delta = lamma(1)-lamma(2)
                return(np.exp(lamma(i)*NUT1/3)/(2*delta))
            delta = lamma(1)-lamma(2)
            A = ((X(1)*(R1+lamma(1))*(R1-lamma(2)))/(2*lamma(1)))-X(3)*delta-((X(2)*(R1+lamma(2))*(R1-lamma(1)))/(2*lamma(2)))+((1)/(1-R1))
            B = X(1)*(R1-lamma(2))-X(2)*(R1-lamma(1))+X(3)*delta
            C = X(2)*(3*R1+lamma(1))-X(1)*(3*R1+lamma(2))+X(3)*delta
            return(R1**(-1)*(1-C/(A*C+B**2)))
        if self.tipo == 'p':
            NUT1 = Y[0]
            return TEMAE3_P1(NUT1)[0]
        if self.tipo == 'pnut':
            P1 = Y[0]
            mod = lambda NUT: TEMAE3_P1(NUT) - P1
            return fsolve(mod,chute_inicial)[0]
        else:
            raise Exception(excep2.format(self.tipo))
    def TEMAE4(self,Y,chute_inicial=1):
        '''P-NUT, TEMA E, 4 passes. Permite especificar um valor para o chute inicial na determinação de PNUT1 = f(P1, R1). Padrão = 1'''
        R1 = Y[1]
        def TEMAE4_P1(NUT1):
            D = (4+R1**2)**(1/2)
            B = math.tanh(R1*NUT1/4)
            A = coth(D*NUT1/4)
            return (4*(2*(1+R1)+D*A+R1*B)**(-1))
        if self.tipo == 'p':
            NUT1 = Y[0]
            return TEMAE4_P1(NUT1)[0]
        if self.tipo == 'pnut':
            P1 = Y[0]
            mod = lambda NUT: TEMAE4_P1(NUT) - P1
            return fsolve(mod,chute_inicial)[0]
        else:
            raise Exception(excep2.format(self.tipo))
    def TEMAG1(self,Y,chute_inicial=1):
        '''P-NUT, TEMA G, 1 passe. Permite especificar um valor para o chute inicial na determinação de PNUT1 = f(P1, R1). Padrão = 1'''
        R1 = Y[1]
        def TEMAG1_P1(NUT1):
            A = (1/(1+R1))*(1-math.exp((-NUT1*(1+R1))/2))
            D = math.exp(-NUT1*(1-R1)/2)
            B = (1-D)/(1-R1*D)
            return (A + B - A*B*(1+R1)+R1*(A*B**2))
        if self.tipo == 'p':
            NUT1 = Y[0]
            return TEMAG1_P1(NUT1)[0]
        if self.tipo == 'pnut':
            P1 = Y[0]
            mod = lambda NUT: TEMAG1_P1(NUT) - P1
            return fsolve(mod,chute_inicial)[0]
        else:
            raise Exception(excep2.format(self.tipo))
    def TEMAG2(self,Y,chute_inicial=1):
        '''P-NUT, TEMA G, 2 passes. Permite especificar um valor para o chute inicial na determinação de PNUT1 = f(P1, R1). Padrão = 1'''
        R1 = Y[1]
        def TEMAG2_P1(NUT1):
            alfa = math.exp(-NUT1*(2+R1)/4)
            beta = math.exp(-NUT1*(2-R1)/2)
            A = (-2*R1*(1-alfa)**2)/(2+R1)
            B = (4-beta*(2+R1))/(2-R1)
            return ((B-alfa**2)/(A+2+R1*B))
        if self.tipo == 'p':
            NUT1 = Y[0]
            return TEMAG2_P1(NUT1)[0]
        if self.tipo == 'pnut':
            P1 = Y[0]
            mod = lambda NUT: TEMAG2_P1(NUT) - P1
            return fsolve(mod,chute_inicial)[0]
        else:
            raise Exception(excep2.format(self.tipo))
    def TEMAH1(self,Y,chute_inicial=1):
        '''P-NUT, TEMA H, 1 passe. Permite especificar um valor para o chute inicial na determinação de PNUT1 = f(P1, R1). Padrão = 1'''
        R1 = Y[1]
        def TEMAH1_P1(NUT1):
            D = math.exp(-NUT1*(1-R1/2)/2)
            A = ((1)/(1+R1/2))*(1-math.exp(-NUT1*(1+R1/2)/2))
            B = (1-D)/(1-R1*D/2)
            E = (A+B-A*B*R1/2)/2
            return (E*(1+(1-B*R1/2)*(1-A*R1/2+A*B*R1))-A*B*(1-B*R1/2))
        if self.tipo == 'p':
            NUT1 = Y[0]
            return TEMAH1_P1(NUT1)[0]
        if self.tipo == 'pnut':
            P1 = Y[0]
            mod = lambda NUT: TEMAH1_P1(NUT) - P1
            return fsolve(mod,chute_inicial)[0]
        else:
            raise Exception(excep2.format(self.tipo))
    def TEMAH2(self,Y,chute_inicial=1):
        '''P-NUT, TEMA H, 2 passes. Permite especificar um valor para o chute inicial na determinação de PNUT1 = f(P1, R1). Padrão = 1'''
        R1 = Y[1]
        def TEMAH2_P1(NUT1):
            beta = (NUT1*(4-R1))/8
            alfa = (NUT1*(4+R1))/8
            D = (1-math.exp(-alfa))/((4/R1)+1)
            E = (1-math.exp(-beta))/((4/R1)-1)
            H = (1-math.exp(-2*beta))/((4/R1)-1)
            G = (1-D)**(2)*(D**(2)+E**(2))+D**(2)*(1+E)**2
            B = (1+H)*(1+E)**2
            return ((1/R1)*(1-(1-D)**4/(B-4*G/R1)))
        if self.tipo == 'p':
            NUT1 = Y[0]
            return TEMAH2_P1(NUT1)[0]
        if self.tipo == 'pnut':
            P1 = Y[0]
            mod = lambda NUT: TEMAH2_P1(NUT) - P1
            return fsolve(mod,chute_inicial)[0]
        else:
            raise Exception(excep2.format(self.tipo))
    def TEMAJ1(self,Y,chute_inicial=1):
        '''P-NUT, TEMA J, 1 passe. Permite especificar um valor para o chute inicial na determinação de PNUT1 = f(P1, R1). Padrão = 1'''
        R1 = Y[1]
        def TEMAJ1_P1(NUT1):
            A = np.exp(NUT1)
            B = np.exp(-NUT1*R1/2)
            return(R1**(-1)*(1-((2-R1)*(2*A+R1*B))/((2+R1)*(2*A-R1/B))))
        if self.tipo == 'p':
            NUT1 = Y[0]
            return TEMAJ1_P1(NUT1)[0]
        if self.tipo == 'pnut':
            P1 = Y[0]
            mod = lambda NUT: TEMAJ1_P1(NUT) - P1
            return fsolve(mod,chute_inicial)[0]
        else:
            raise Exception(excep2.format(self.tipo))
    def TEMAJ2(self,Y,chute_inicial=1):
        '''P-NUT, TEMA J, 2 passes. Permite especificar um valor para o chute inicial na determinação de PNUT1 = f(P1, R1). Padrão = 1'''
        R1 = Y[1]
        def TEMAJ2_P1(NUT1):
            A = np.exp(NUT1)
            lamma = (1+R1**(2)/4)**(1/2)
            B = (pow(A,lamma)+1)/(pow(A,lamma)-1)
            C = (pow(A,(1+lamma)/2))/(lamma-1+(1+lamma)*pow(A,lamma))
            D = 1+(lamma*pow(A,(1-lamma)/2))/(pow(A,lamma)+1)
            return (1+R1/2+lamma*B-2*lamma*C*D)**(-1)
        if self.tipo == 'p':
            NUT1 = Y[0]
            return TEMAJ2_P1(NUT1)
        if self.tipo == 'pnut':
            P1 = Y[0]
            mod = lambda NUT: TEMAJ2_P1(NUT) - P1
            return fsolve(mod,chute_inicial)[0]
        else:
            raise Exception(excep2.format(self.tipo))
    def TEMAJ4(self,Y,chute_inicial=1):
        '''P-NUT, TEMA J, 4 passes. Permite especificar um valor para o chute inicial na determinação de PNUT1 = f(P1, R1). Padrão = 1'''
        R1 = Y[1]
        def TEMAJ4_P1(NUT1):
            A = np.exp(NUT1)
            lamma = (1+R1**(2)/16)**(1/2)
            B = (pow(A,lamma)+1)/(pow(A,lamma)-1)
            C = (pow(A,(1+lamma)/2))/(lamma-1+(1+lamma)*pow(A,lamma))
            D = 1+(lamma*pow(A,(1-lamma)/2))/(pow(A,lamma)+1)
            E = np.exp(R1*NUT1/2)
            return (1+R1/4*((1+3*E)/(1+E))+lamma*B-2*lamma*C*D)**(-1)
        if self.tipo == 'p':
            NUT1 = Y[0]
            return TEMAJ4_P1(NUT1)[0]
        if self.tipo == 'pnut':
            P1 = Y[0]
            mod = lambda NUT: TEMAJ4_P1(NUT) - P1
            return fsolve(mod,chute_inicial)[0]
        else:
            raise Exception(excep2.format(self.tipo))