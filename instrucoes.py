# Instruções de uso

# Antes de tudo, importe a biblioteca 'trocadores', a seguir
from trocadores import *

# É necessário criar um vetor [a,b] para guardar na memória os valores utilizados nos cálculos
# Para isso, utilizaremos a função dados(a,b)
help(dados)
# "Retorna uma array 2x1 na forma [a,b], sendo a = P1 ou NUT1 (a depender do caso), e b = R1."

# Se eu quero calcular P1, então a = NUT1; se eu quero calcular NUT1, então a = P1
# Vamos salvar os dados em uma variável "A"
A = dados(0.5, 0.4)
# Para este exemplo, temos um trocador de calor tipo TEMA J 4 passes, e queremos calcular NUT1
# Primeiramente, devemos criar um objeto da classe Trocador, salvo na memória como "tc1"

# Devemos inserir como argumento o tipo de cálculo que pretendemos fazer
# No nosso caso, queremos calcular NUT1 pelo método P-NUT. Logo, inserimos 'pnut' no argumento
tc1 = Trocador('pnut')
# Caso quiséssemos calcular P1, colocaríamos 'p' no lugar de 'pnut'

# Para fazer os cálculos de NUT1, precisamos especificar o tipo de trocador de calor
# Além disso, precisamos inserir os valores de P1 e R1, que já foram armazenados na variável "A" anteriormente
# Armazenaremos o valor de NUT1 em uma variável chamada "NUT1"
NUT1 = tc1.TEMAJ4(A)
print('NUT1', NUT1)

# O cálculo é um processo iterativo que parte de um chute inicial para NUT1
# O chute inicial é definido por padrão como NUT1 = 1.0, como é possível ver aqui:
help(Trocador.TEMAJ4)
# "Permite especificar um valor para o chute inicial na determinação de PNUT1 = f(P1, R1). Padrão = 1"
# Caso queiramos alterar o valor do chute inicial, basta especificá-lo como um segundo argumento:
NUT1_alternativo = tc1.TEMAJ4(A,25.0)
print('NUT1 com chute 25.0 = ', NUT1_alternativo)

# A depender do valor do chute inicial, o processo para o cálculo de NUT1 pode convergir em outro lugar
# Recomendo não mudar o valor do chute inicial, a não ser que você saiba o que está fazendo
# Para ver a lista de trocadores de calor incluídos nesta biblioteca, basta usar a função help():
help(Trocador)