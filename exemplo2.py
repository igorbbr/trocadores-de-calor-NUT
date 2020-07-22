from trocadores import *

# Exemplo 2 (Cap. 10)

class Fluido:
    def __init__(self,fluido):
        self.fluido = fluido    
g = Fluido('gas')
ag = Fluido('water')

def tprint(texto):
    print('--- [ ',texto,' ] ---\n')

def pprint(varname,var,unit=' '):
    print(varname,'= \n',var,unit,'\n')

# Dados do enunciado

g.cp = 1009.0
g.Te = 450.0 + 273.15 # K
g.Ts = 130.0 + 273.15 # K

ag.cp = 4190.0 # J/kgK
ag.m = 2.0 # kg/s
ag.Te = 30.0 + 273.15 # K
ag.Ts = 125.0 + 273.15 # K

Uq = 200.0 # W/m²K

# Balanço de energia (água = f2)
tprint('Balanço de energia')

qtc = ag.m*ag.cp*(ag.Ts-ag.Te)
pprint('qtc', qtc, 'W')

g.m = qtc/(g.cp*(g.Te-g.Ts))
pprint('mg',g.m,'kg/s')

# Capacidade térmica (gás = quente, água = frio)
tprint('Capacidades térmicas')

Ch = g.m*g.cp
pprint('Ch',Ch,'W/K')
Cc = ag.m*ag.cp
pprint('Cc',Cc,'W/K')

# Método eps-NUT
tprint('Método eps-NUT')

Cmin=min(Cc,Ch)
Cmax=max(Cc,Ch)
Cr = Cmin/Cmax
pprint('Cr',Cr)

qmax = Cmin*(g.Te-ag.Te)
pprint('qmax',qmax,'W')

eps = qtc/qmax
pprint('eps',eps)

tc = Trocador('pnut')
data = dados(eps,Cr)
nut = tc.TC1053(data)
pprint('NUT',nut)

Aq = nut*Cmin/Uq
pprint('Aq',Aq,'m²')
