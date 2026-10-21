# PP2 — Cálculo Numérico

Requer Python 3.10 ou superior, sem bibliotecas externas.

Na pasta do projeto, execute:

```bash
python main.py
```

Se estiver na pasta M2, use `python PP2/main.py`.

- `main.py`: dicionários com os dados e exemplos dos 11 métodos.
- `metodos.py`: funções dos métodos de resolução.
- `auxiliares.py`: validações, cálculos de apoio, tabelas e relatórios.
- `resultado.md`: passo a passo da última execução válida.

O `main.py` executa todos os exemplos em sequência. Para resolver somente
uma questão, mantenha a chamada desejada dentro de `if __name__ == "__main__":`.
Cada chamada exibe o resultado no terminal e sobrescreve `resultado.md`;
este README permanece intacto. Importar os módulos não executa exemplos.

## Sistemas lineares

Declare um dicionário com `matriz` e `termos`, como no exemplo do `main.py`:

```python
SISTEMA_2 = {
    "matriz": [
        [10, 2, 1],
        [1, 5, 1],
        [2, 3, 10],
    ],
    "termos": [7, -8, 6],
}
```

Cada linha de `matriz` representa uma equação, e cada coluna corresponde
à mesma incógnita em todas as linhas. Inclua zero para termos ausentes.
A matriz deve ser quadrada, com um termo independente por equação.
O exemplo acima tem solução `(1, -2, 1)`.

Importe as funções utilizadas de `metodos`, por exemplo:

```python
from metodos import gauss, gauss_seidel
```

Passe o dicionário diretamente, sem `**`. As opções do método vêm depois dele.

### Gauss

Transforma o sistema em triangular e resolve por substituição retroativa.

```python
gauss(SISTEMA_2)
gauss(SISTEMA_2, pivotamento="parcial")
```

`pivotamento` aceita `"nenhum"`, `"parcial"` ou `"completo"`; o padrão é
`"nenhum"`, com troca de linhas quando o pivô é numericamente nulo.

### Gauss-Jordan

Reduz a matriz à identidade para obter diretamente a solução.

```python
jordan(SISTEMA_2)
```

Use o mesmo dicionário com `matriz` e `termos`.

### Decomposição LU

Decompõe a matriz em duas matrizes triangulares para resolver o sistema.

```python
decomposicao_lu(SISTEMA_2)
```

Use o mesmo formato de sistema. L é triangular inferior e U é triangular
superior; a implementação usa pivotação parcial e calcula `P*A = L*U`.

### Cramer

Calcula cada incógnita pela razão entre determinantes.

```python
cramer(SISTEMA_2)
```

Use o mesmo formato de sistema, com determinante não nulo.

### Gauss-Seidel

Aproxima a solução por iterações, utilizando imediatamente os valores atualizados.

```python
gauss_seidel(SISTEMA_2, tolerancia=1e-8, max_iteracoes=100)
```

Recebe o mesmo dicionário de sistema. `inicial=[0, 0, 0]` permite escolher
os valores iniciais; se omitido, começa com zeros. `tolerancia` controla
a diferença máxima entre iterações, e `max_iteracoes` limita as repetições.
Os padrões são `0.01` e `10`, respectivamente. Use `casas_decimais=4` apenas
quando quiser arredondar cada atualização; o padrão é não arredondar.
Confira `convergiu` e o resíduo no retorno: esse critério de parada sozinho
não garante o erro real da solução.

Quando o exercício pedir distância relativa, acrescente `erro_relativo=True`:

```python
gauss_seidel(SISTEMA_2, tolerancia=0.001, erro_relativo=True)
```

Nesse modo, a parada usa `max|x(k)-x(k-1)| / max|x(k)|`. O histórico conserva
`diferenca` absoluta e também informa `diferenca_relativa`. O padrão continua
sendo a distância absoluta.

### Substituição direta

Resolve sistemas triangulares inferiores, da primeira equação à última.

```python
SISTEMA_TRIANGULAR_INFERIOR = {
    "matriz": [[2, 0, 0], [1, 3, 0], [4, -1, 2]],
    "termos": [2, 7, 8],
}

substituicao_direta(SISTEMA_TRIANGULAR_INFERIOR)
```

Os coeficientes acima da diagonal devem ser zero, e a diagonal deve ser
não nula. A solução deste exemplo é `(1, 2, 3)`.

### Substituição retroativa

Resolve sistemas triangulares superiores, da última equação à primeira.

```python
SISTEMA_TRIANGULAR_SUPERIOR = {
    "matriz": [[2, -1, 3], [0, 4, 2], [0, 0, 5]],
    "termos": [9, 14, 15],
}

substituicao_retroativa(SISTEMA_TRIANGULAR_SUPERIOR)
```

Os coeficientes abaixo da diagonal devem ser zero, e a diagonal deve ser
não nula. A solução deste exemplo também é `(1, 2, 3)`.

## Integração numérica

Trapézios e Simpson aceitam o mesmo formato de dicionário. Para uma função:

```python
INTEGRAL = {
    "funcao": lambda x: x**2,
    "intervalo": (0, 2),
    "subintervalos": 4,
}
```

`funcao` calcula o valor de y para cada x, `intervalo` contém os limites
inicial e final, e `subintervalos` define em quantas partes dividir o intervalo.
O limite final deve ser maior que o inicial. Funções de `math` podem ser
usadas após a importação; as trigonométricas trabalham em radianos.

### Trapézios

Aproxima uma integral pela soma das áreas de trapézios.

```python
trapezios(INTEGRAL)
```

Use pelo menos um subintervalo.

### Simpson

Aproxima uma integral usando parábolas a cada dois subintervalos.

```python
simpson(INTEGRAL)
```

Use um número par de subintervalos, com pelo menos dois.

### Trapézios com precisão

Refina os trapézios até a diferença entre aproximações atender à tolerância.

```python
trapezios_com_precisao(
    INTEGRAL["funcao"], INTEGRAL["intervalo"], casas_decimais=4
)
```

Essa função recebe a função e o intervalo separadamente, aproveitando
os valores do dicionário. Ela escolhe a quantidade de subintervalos;
o campo `subintervalos` não é utilizado. `casas_decimais` define o erro
solicitado como `0.5 * 10**(-casas_decimais)`. O refinamento dobra a malha
e usa a estimativa de Richardson, com margem de segurança de 2. Se a última
malha atingir o limite sem dobrar, a estimativa considera a razão efetiva.

### Simpson com precisão

Refina Simpson até a diferença entre aproximações atender à tolerância.

```python
simpson_com_precisao(
    INTEGRAL["funcao"], INTEGRAL["intervalo"], casas_decimais=4
)
```

Os parâmetros funcionam como no refinamento dos trapézios. Ambos aceitam
`max_subintervalos=1000` para limitar o refinamento e `unidade="u.a."` para
rotular o resultado; esses são os valores padrão. Confira `convergiu` e
`erro_estimado` no retorno. Sem uma cota de derivada, a estimativa é empírica
e `erro_garantido` é `False`: não prova o erro real para uma função arbitrária.

### Erro solicitado e menor número de subintervalos pela cota

Use `erro_maximo` quando o enunciado fornecer epsilon diretamente. Para
usar a cota teórica, informe também `limite_derivada`: um limite superior
de `|f''(x)|` em trapézios ou `|f''''(x)|` em Simpson, válido em todo o intervalo.

```python
from math import exp

trapezios_com_precisao(
    lambda x: exp(-x), (1, 2),
    erro_maximo=1e-5, limite_derivada=exp(-1),
)
simpson_com_precisao(
    lambda x: exp(-x), (1, 2),
    erro_maximo=1e-5, limite_derivada=exp(-1),
)
```

Os resultados usam, respectivamente, **56 e 4 subintervalos**. São os menores
n que satisfazem as cotas `M2*(b-a)**3/(12*n**2)` e
`M4*(b-a)**5/(180*n**4)`, com n par em Simpson. Isso não significa o menor
n pelo erro real, que pode ser menor que a cota suficiente.

O retorno inclui `erro_estimado`, `erro_maximo` e `erro_garantido`.
A garantia é da cota de truncamento e depende da regularidade da função
e do limite de derivada correto; não inclui arredondamento de `float`.
Se `max_subintervalos` impedir atingir a cota, `convergiu` e `erro_garantido`
serão `False`. Um limite de derivada zero permite n=1 em trapézios ou n=2
em Simpson, quando a fórmula é exata para o polinômio.

### Integrais impróprias em um extremo

Declare `impropria="inicio"` ou `impropria="fim"` para não avaliar a função
no ponto singular. As mesmas funções de trapézios e Simpson integram faixas
sucessivamente menores e estimam a contribuição que falta junto ao extremo.

```python
from math import sqrt, tan, pi

IMPROPRIA = {
    "funcao": lambda x: 6 / sqrt(3 * x),
    "intervalo": (0, 9),
    "impropria": "inicio",
    "erro_maximo": 1e-6,
    "subintervalos": 10,
}

simpson(IMPROPRIA)
trapezios(IMPROPRIA)
```

O resultado aproxima `12*sqrt(3)`, cerca de **20.784609690827**.
`subintervalos` é a malha inicial de cada faixa; ela pode aumentar para
controlar a quadratura. Não é a aplicação da fórmula fechada original
com somente 10 pontos de divisão: aquela fórmula não está definida no polo.
Os padrões são 16 subintervalos iniciais, `max_subintervalos=16384` por faixa
e `max_refinamentos=60` faixas. Não misture esse modo com tabelas ou perfis.

```python
DIVERGENTE = {
    "funcao": tan,
    "intervalo": (1, pi / 2),
    "impropria": "fim",
}

simpson(DIVERGENTE)
```

Essa integral diverge. O teste numérico identifica parcelas que não decaem
e retorna `integral=None`, `convergiu=False` e `possivel_divergencia=True`.
`aproximacao` guarda somente a soma truncada, que não deve ser usada como
resposta. Se apenas faltar precisão, o resultado também fica inconclusivo;
isso não constitui uma prova de divergência.

A estimativa de cauda é empírica; `erro_garantido` permanece `False` neste
modo. Uma função arbitrária pode ter comportamento que as amostras não
revelam. O modo regular também recusa polos suspeitos nos extremos, em vez
de aceitar silenciosamente números enormes causados por arredondamento.
Essa sondagem é uma proteção numérica, não uma análise completa do domínio.

## Outros dicionários para trapézios e Simpson

### Tabela de pontos e valores

```python
TABELA = {
    "pontos": [0, 1, 2, 3, 4],
    "valores": [0, 1, 4, 9, 16],
    "unidade": "u.a.",
}

trapezios(TABELA)
simpson(TABELA)
```

Os pontos devem ser crescentes e igualmente espaçados, com um valor por
ponto. Para n subintervalos, informe n+1 pontos. Simpson continua exigindo
n par. `unidade` é opcional e apenas rotula o resultado.

### Perfis de alturas, como na lagoa e no terreno

```python
AREA = {
    "area_base": 100,
    "perfis": [
        {"alturas": [0, 3, 4, 3, 0], "passo": 2, "sinal": -1},
    ],
    "unidade": "km²",
}

simpson(AREA)
```

`alturas` contém as medidas, incluindo os extremos; `passo` é a distância
positiva entre elas. Cada perfil pode ter um passo diferente.
`sinal=1` soma a contribuição e `sinal=-1` desconta; o padrão é 1.
`area_base` é opcional e vale zero por padrão. O total é a área base mais
as contribuições dos perfis. A mesma entrada funciona em `trapezios(AREA)`.

### Vários trechos de tabelas ou funções

```python
DADOS = {
    "trechos": [
        {"pontos": [0, 1, 2], "valores": [0, 1, 4]},
        {"funcao": lambda x: x**2, "intervalo": (2, 4), "subintervalos": 4},
    ],
}

simpson(DADOS)
```

Cada trecho deve respeitar as condições do método. Os resultados são
somados; `sinal=-1` também permite descontar trechos. Não misture `trechos`
e `perfis` no mesmo dicionário, nem dados de tabela ou função fora deles.
Mantenha as unidades coerentes: tempos em minutos devem ser divididos
por 60 antes de integrar velocidades em km/h para obter distâncias em km.

## Retornos e funções de apoio

Todos os métodos retornam um dicionário com `relatorio`. Para guardar o
resultado e acessar os números:

```python
resultado = gauss(SISTEMA_2)
solucao = resultado["solucao"]
residuo = resultado["residuo"]

resultado = simpson(INTEGRAL)
integral = resultado["integral"]
```

Nas substituições triangulares, o retorno contém `solucao` e `relatorio`.
Gauss-Seidel também informa `convergiu`, `iteracoes` e `historico`;
os métodos com precisão informam `convergiu`, `subintervalos`, `historico`
e os campos de erro descritos acima.
Com perfis ou trechos, os cálculos individuais ficam em `perfis` ou `trechos`.

Também é possível importar de `auxiliares` estas funções de apoio:

```python
from auxiliares import calcular_residuo, classificar_sistema, criterios_convergencia

classificar_sistema(SISTEMA_2)
criterios_convergencia(SISTEMA_2)
calcular_residuo(SISTEMA_2, solucao=[1, -2, 1])
```

Elas classificam o sistema, verificam condições suficientes de convergência
e calculam `r = b - A*x`, respectivamente.

## Verificação

Para executar a suíte de testes a partir da pasta do projeto:

```bash
python -B -m unittest discover -s tests -v
```

Os testes conferem resultados anteriores, cotas de erro, refinamento,
distância relativa, integrais impróprias e preservação do relatório em erros.
