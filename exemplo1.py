from trocadores import *
help(Trocador)
class Fluido:
    def __init__(self, fluid):
        self.fluid = fluid

oleo = Fluido('PNF2')
agua = Fluido('H2O')

# Exemplo 1 (Cap. 10) - Método P-NUT

# Dados do enunciado

oleo.Tent = 120.0 + 273.15
oleo.C = 4200.0
oleo.cp = 2.35e3

agua.Tent = 20.0 + 273.15
agua.vazao = 0.8
agua.cp = 4.18e3

U = 600.0
A = 10.0
razao = 1.4

# Resolução

Ch = oleo.C
Cc = agua.vazao * agua.cp

# eps-NUT
print('--- [ Método eps-NUT ] ---\n')

Cmin = min(Ch, Cc)
Cmax = max(Ch, Cc)
Cr = Cmin/Cmax
NUT = U*A/Cmin

tc0 = Trocador('eps')
valores0 = dados(NUT, Cr)

epsp = tc0.TCparalelo(valores0)
print('eps (paralelo) = \n', epsp, '\n')
qmaxp = Cmin*(oleo.Tent - agua.Tent)
qtcp = epsp*qmaxp
print('qtc (paralelo) = \n', qtcp, 'W\n')

epsc = tc0.TCcontra(valores0)
print('eps (contra) = \n', epsc, '\n')
qmaxc = Cmin*(oleo.Tent - agua.Tent)
qtcc = epsc*qmaxc
print('qtc (contra) = \n', qtcc, 'W\n')

# P-NUT
print('--- [ Método P-NUT ] ---\n')

R1 = Cc/Ch
NUT1 = U*A/Cc

valores1 = dados(NUT1, R1)

tc1 = Trocador('p')
P1p = tc1.TC1052(valores1)
print('P1 (paralelo) = \n', P1p, '\n')

P1c = tc1.TC1051(valores1)
print('P1 (contra) = \n', P1c, '\n')




