from math import ceil, isclose, isfinite
from pathlib import Path


ARQUIVO_RESULTADO = Path(__file__).resolve().with_name("resultado.md")
TOLERANCIA_PIVO = 1e-12
TOLERANCIA_ITERACAO = 0.01
MAX_ITERACOES = 10
CASAS_EXIBICAO = 10
CAMPOS_TABELA = {"pontos", "valores"}
CAMPOS_FUNCAO = {"funcao", "intervalo", "subintervalos"}
CAMPOS_PERFIL = {"alturas", "passo"}
PASSO_SONDAGEM = 1e-6
RAZAO_SUSPEITA_POLO = 1e8
RAZAO_CRESCIMENTO_POLO = 1.5
MIN_PARCELAS_DIVERGENCIA = 12
RAZAO_SEM_DECAIMENTO = 1 - 1e-10


def _numero(valor):
    return format(valor, f".{CASAS_EXIBICAO}g")


def _tabela(cabecalho, linhas):
    linhas = [cabecalho, ["---"] * len(cabecalho), *linhas]
    return "\n".join(
        "| " + " | ".join(
            _numero(valor) if type(valor) in (int, float) else str(valor)
            for valor in linha
        ) + " |" for linha in linhas
    )


def _tabela_matriz(matriz, termos=None, variaveis=None):
    variaveis = range(len(matriz[0])) if variaveis is None else variaveis
    cabecalho = [f"x{coluna + 1}" for coluna in variaveis]
    linhas = [list(linha) for linha in matriz]
    if termos is not None:
        cabecalho.append("b")
        linhas = [linha + [termo] for linha, termo in zip(linhas, termos)]
    return _tabela(cabecalho, linhas)


def _publicar(resultado, partes):
    relatorio = "\n\n".join(partes) + "\n"
    ARQUIVO_RESULTADO.write_text(relatorio, encoding="utf-8")
    print(relatorio)
    resultado["relatorio"] = relatorio
    return resultado


def _preparar_sistema(matriz, termos):
    coeficientes = [list(map(float, linha)) for linha in matriz]
    independentes = list(map(float, termos))
    ordem = len(coeficientes)
    if not ordem or any(len(linha) != ordem for linha in coeficientes):
        raise ValueError("A matriz deve ser quadrada e não vazia.")
    if len(independentes) != ordem:
        raise ValueError("O vetor de termos independentes deve ter um valor por linha.")
    valores = [valor for linha in coeficientes for valor in linha] + independentes
    if not all(map(isfinite, valores)):
        raise ValueError("Os coeficientes e termos devem ser números finitos.")
    return coeficientes, independentes


def _calcular_residuo(matriz, termos, solucao):
    return [
        termo - sum(coeficiente * valor for coeficiente, valor in zip(linha, solucao))
        for linha, termo in zip(matriz, termos)
    ]


def calcular_residuo(dados, solucao):
    matriz, termos = _preparar_sistema(dados["matriz"], dados["termos"])
    solucao = list(map(float, solucao))
    if len(solucao) != len(termos) or not all(map(isfinite, solucao)):
        raise ValueError("A solução deve conter um número finito por incógnita.")
    residuo = _calcular_residuo(matriz, termos, solucao)
    return _publicar(
        {"residuo": residuo, "norma_maxima": max(map(abs, residuo))},
        ["# Resíduo", "r = b - A*x", _tabela(["Linha", "r"], [
            [indice + 1, valor] for indice, valor in enumerate(residuo)
        ]), f"Norma máxima: {_numero(max(map(abs, residuo)))}."],
    )


def _substituir(matriz, termos, inferior=False, variaveis=None):
    ordem = len(termos)
    variaveis = range(ordem) if variaveis is None else variaveis
    solucao = [0.0] * ordem
    indices = range(ordem) if inferior else range(ordem - 1, -1, -1)
    etapas = []
    for linha in indices:
        colunas = range(linha) if inferior else range(linha + 1, ordem)
        soma = sum(matriz[linha][coluna] * solucao[coluna] for coluna in colunas)
        pivo = matriz[linha][linha]
        if pivo == 0:
            raise ValueError(f"Pivô nulo na linha {linha + 1}; não há solução única.")
        solucao[linha] = (termos[linha] - soma) / pivo
        etapas.append(
            f"x{variaveis[linha] + 1} = "
            f"({_numero(termos[linha])} - ({_numero(soma)})) / "
            f"{_numero(pivo)} = {_numero(solucao[linha])}"
        )
    return solucao, etapas


def _eliminar(matriz, termos, pivotamento):
    if pivotamento not in ("nenhum", "parcial", "completo"):
        raise ValueError("Pivotamento deve ser 'nenhum', 'parcial' ou 'completo'.")
    superior = [linha.copy() for linha in matriz]
    transformados = termos.copy()
    ordem = len(termos)
    inferior = [[float(i == j) for j in range(ordem)] for i in range(ordem)]
    permutacao = list(range(ordem))
    variaveis = list(range(ordem))
    etapas = ["## Matriz aumentada inicial", _tabela_matriz(matriz, termos)]
    for coluna in range(ordem):
        linha_pivo, coluna_pivo = coluna, coluna
        if pivotamento == "completo":
            linha_pivo, coluna_pivo = max(
                ((i, j) for i in range(coluna, ordem) for j in range(coluna, ordem)),
                key=lambda posicao: abs(superior[posicao[0]][posicao[1]]),
            )
        elif pivotamento == "parcial" or abs(superior[coluna][coluna]) <= (
            TOLERANCIA_PIVO * max(abs(linha[variaveis[coluna]]) for linha in matriz)
        ):
            linha_pivo = max(
                range(coluna, ordem), key=lambda i: abs(superior[i][coluna])
            )
        escala_coluna = max(abs(linha[variaveis[coluna_pivo]]) for linha in matriz)
        if abs(superior[linha_pivo][coluna_pivo]) <= TOLERANCIA_PIVO * escala_coluna:
            raise ValueError(
                "Pivô numericamente nulo: não foi possível obter solução única. "
                "Use classificar_sistema()."
            )
        if linha_pivo != coluna:
            for estrutura in (superior, transformados, permutacao):
                estrutura[coluna], estrutura[linha_pivo] = (
                    estrutura[linha_pivo], estrutura[coluna]
                )
            origem, destino = inferior[coluna], inferior[linha_pivo]
            origem[:coluna], destino[:coluna] = destino[:coluna], origem[:coluna]
            etapas.extend([
                f"Troca: L{coluna + 1} <-> L{linha_pivo + 1}.",
                _tabela_matriz(superior, transformados, variaveis),
            ])
        if coluna_pivo != coluna:
            for linha in superior:
                linha[coluna], linha[coluna_pivo] = linha[coluna_pivo], linha[coluna]
            variaveis[coluna], variaveis[coluna_pivo] = (
                variaveis[coluna_pivo], variaveis[coluna]
            )
            etapas.extend([
                f"Troca de colunas: C{coluna + 1} <-> C{coluna_pivo + 1}. "
                f"Ordem das incógnitas: {[f'x{i + 1}' for i in variaveis]}.",
                _tabela_matriz(superior, transformados, variaveis),
            ])
        pivo = superior[coluna][coluna]
        for linha in range(coluna + 1, ordem):
            if superior[linha][coluna] == 0:
                continue
            multiplicador = -superior[linha][coluna] / pivo
            inferior[linha][coluna] = -multiplicador
            for indice in range(coluna + 1, ordem):
                superior[linha][indice] += multiplicador * superior[coluna][indice]
            superior[linha][coluna] = 0.0
            transformados[linha] += multiplicador * transformados[coluna]
            etapas.extend([
                f"m{linha + 1},{coluna + 1} = -a{linha + 1},{coluna + 1} / "
                f"a{coluna + 1},{coluna + 1} = {_numero(multiplicador)}; "
                f"L{linha + 1} <- L{linha + 1} + "
                f"({_numero(multiplicador)}) * L{coluna + 1}.",
                _tabela_matriz(superior, transformados, variaveis),
            ])
    return superior, transformados, inferior, permutacao, variaveis, etapas


def _publicar_sistema(matriz, termos, solucao, etapas, rodape=(), **resultado):
    residuo = _calcular_residuo(matriz, termos, solucao)
    if not all(map(isfinite, solucao + residuo)):
        raise ArithmeticError("A solução ou o resíduo excedeu a faixa numérica.")
    partes = [
        *etapas,
        "## Solução", _tabela(["Incógnita", "Valor"], [
            [f"x{i + 1}", valor] for i, valor in enumerate(solucao)
        ]), "## Resíduo: r = b - A*x",
        _tabela(["Linha", "b", "A*x", "r"], [
            [i + 1, termo, termo - residuo[i], residuo[i]]
            for i, termo in enumerate(termos)
        ]), f"Norma máxima do resíduo: {_numero(max(map(abs, residuo)))}.",
        *rodape,
    ]
    return _publicar({"solucao": solucao, "residuo": residuo, **resultado}, partes)


def classificar_sistema(dados):
    matriz, termos = _preparar_sistema(dados["matriz"], dados["termos"])
    ordem = len(termos)
    aumentada = []
    for linha, termo in zip(matriz, termos):
        escala = max(map(abs, linha)) or abs(termo) or 1.0
        aumentada.append([valor / escala for valor in linha] + [termo / escala])
    pivos = []
    for coluna in range(ordem):
        linha_atual = len(pivos)
        linha_pivo = max(
            range(linha_atual, ordem), key=lambda i: abs(aumentada[i][coluna])
        )
        if abs(aumentada[linha_pivo][coluna]) <= TOLERANCIA_PIVO:
            continue
        aumentada[linha_atual], aumentada[linha_pivo] = (
            aumentada[linha_pivo], aumentada[linha_atual]
        )
        pivo = aumentada[linha_atual][coluna]
        aumentada[linha_atual] = [valor / pivo for valor in aumentada[linha_atual]]
        base = aumentada[linha_atual]
        for indice, linha in enumerate(aumentada):
            if indice != linha_atual:
                fator = linha[coluna]
                aumentada[indice] = [a - fator * b for a, b in zip(linha, base)]
        pivos.append(coluna)
    impossivel = any(
        all(abs(valor) <= TOLERANCIA_PIVO for valor in linha[:-1])
        and abs(linha[-1]) > TOLERANCIA_PIVO for linha in aumentada
    )
    classificacao = "impossível" if impossivel else (
        "possível e determinado" if len(pivos) == ordem else "possível e indeterminado"
    )
    particular, base = None, []
    if not impossivel:
        particular = [0.0] * ordem
        for linha, coluna in enumerate(pivos):
            particular[coluna] = aumentada[linha][-1]
        for livre in (coluna for coluna in range(ordem) if coluna not in pivos):
            vetor = [0.0] * ordem
            vetor[livre] = 1.0
            for linha, coluna in enumerate(pivos):
                vetor[coluna] = -aumentada[linha][livre]
            base.append(vetor)
    partes = [
        "# Classificação do sistema", _tabela_matriz(matriz, termos),
        f"Sistema **{classificacao}**; posto(A) = {len(pivos)}; "
        f"posto(A|b) = {len(pivos) + int(impossivel)}.",
        "## Forma reduzida",
        _tabela_matriz(
            [linha[:-1] for linha in aumentada],
            [linha[-1] for linha in aumentada],
        ),
    ]
    if particular is not None:
        partes.append(f"Solução particular: {list(map(_numero, particular))}.")
    if base:
        partes.append(
            "Solução geral: x = particular + t1*v1 + t2*v2 + ...; "
            "parâmetros reais livres."
        )
        partes.extend(
            f"v{i + 1} = {list(map(_numero, vetor))}"
            for i, vetor in enumerate(base)
        )
    return _publicar({"classificacao": classificacao, "solucao_particular": particular,
                      "base_nucleo": base, "posto": len(pivos)}, partes)


def _determinante(matriz):
    coeficientes = [linha.copy() for linha in matriz]
    determinante = 1.0
    for coluna in range(len(matriz)):
        pivo = max(
            range(coluna, len(matriz)),
            key=lambda i: abs(coeficientes[i][coluna]),
        )
        escala_coluna = max(abs(linha[coluna]) for linha in matriz)
        if abs(coeficientes[pivo][coluna]) <= TOLERANCIA_PIVO * escala_coluna:
            return 0.0
        if pivo != coluna:
            coeficientes[coluna], coeficientes[pivo] = (
                coeficientes[pivo], coeficientes[coluna]
            )
            determinante = -determinante
        diagonal = coeficientes[coluna][coluna]
        determinante *= diagonal
        for linha in range(coluna + 1, len(matriz)):
            fator = coeficientes[linha][coluna] / diagonal
            for indice in range(coluna + 1, len(matriz)):
                coeficientes[linha][indice] -= fator * coeficientes[coluna][indice]
    return determinante


def _criterios_convergencia(matriz):
    betas, linhas = [], []
    for i, linha in enumerate(matriz):
        diagonal = abs(linha[i])
        if diagonal == 0:
            raise ValueError(f"Diagonal nula na linha {i + 1}; reordene as equações.")
        soma = sum(abs(valor) for j, valor in enumerate(linha) if j != i)
        linhas.append({
            "diagonal": diagonal, "soma_fora_diagonal": soma,
            "atende": diagonal >= soma,
        })
        beta = (sum(abs(linha[j]) * betas[j] for j in range(i))
                + sum(abs(linha[j]) for j in range(i + 1, len(matriz)))) / diagonal
        betas.append(beta)
    return {
        "linhas": linhas, "betas": betas,
        "dominancia_estrita": all(
            linha["diagonal"] > linha["soma_fora_diagonal"] for linha in linhas
        ),
        "sassenfeld": max(betas) < 1,
    }


def _relatar_criterios(criterios):
    return ["## Critérios de convergência", _tabela(
        ["Linha", "|a_ii|", "Soma fora da diagonal", "Critério >=", "Beta"],
        [[i + 1, linha["diagonal"], linha["soma_fora_diagonal"],
          "Atende" if linha["atende"] else "Não atende", criterios["betas"][i]]
         for i, linha in enumerate(criterios["linhas"])],
    ), "A aula usa >= no critério das linhas. Dominância estrita em todas as linhas "
        "ou max(beta) < 1 são condições suficientes de convergência. "
        "Igualdade sozinha não garante convergência; "
        "falhar nesses testes não prova divergência."]


def criterios_convergencia(dados):
    matriz = dados["matriz"]
    matriz, _ = _preparar_sistema(matriz, [0] * len(matriz))
    criterios = _criterios_convergencia(matriz)
    return _publicar(criterios, [
        "# Análise de convergência", *_relatar_criterios(criterios)
    ])


def _integrar_dados(pontos, valores, metodo, unidade, titulo):
    pontos, valores = list(map(float, pontos)), list(map(float, valores))
    if len(pontos) != len(valores) or len(pontos) < 2:
        raise ValueError("Forneça ao menos dois pontos, cada um com seu valor.")
    if not all(map(isfinite, pontos + valores)):
        raise ValueError("Os pontos e valores devem ser números finitos.")
    intervalos = len(pontos) - 1
    passo = pontos[1] - pontos[0]
    if passo <= 0 or any(pontos[i + 1] <= pontos[i] for i in range(intervalos)):
        raise ValueError("Os pontos devem estar em ordem estritamente crescente.")
    if any(
        not isclose(
            pontos[i + 1] - pontos[i], passo,
            rel_tol=1e-9, abs_tol=abs(passo) * 1e-9,
        ) for i in range(intervalos)
    ):
        raise ValueError(
            "As fórmulas compostas das aulas exigem pontos igualmente espaçados."
        )
    if metodo == "simpson" and intervalos % 2:
        raise ValueError(
            "Simpson 1/3 exige um número par de subintervalos (ao menos dois)."
        )
    if metodo == "simpson":
        pesos = [
            1 if i in (0, intervalos) else (4 if i % 2 else 2)
            for i in range(intervalos + 1)
        ]
    else:
        pesos = [1] + [2] * (intervalos - 1) + [1]
    divisor = 3 if metodo == "simpson" else 2
    produtos = [peso * valor for peso, valor in zip(pesos, valores)]
    soma_ponderada = sum(produtos)
    integral = passo / divisor * soma_ponderada
    if not isfinite(integral):
        raise ArithmeticError(
            "A integral excedeu a faixa numérica; verifique a escala dos dados."
        )
    linhas = zip(range(intervalos + 1), pontos, valores, pesos, produtos)
    partes = [f"# {titulo}", f"n = {intervalos}; a = {_numero(pontos[0])}; "
              f"b = {_numero(pontos[-1])}; h = (b - a)/n = {_numero(passo)}.",
              _tabela(["i", "x", "y", "cte", "y*cte"], linhas),
              f"Soma ponderada = {_numero(soma_ponderada)}.",
              f"I = h/{divisor} * soma(y*cte) = {_numero(passo)}/{divisor} * "
              f"{_numero(soma_ponderada)} = **{_numero(integral)} {unidade}**."]
    return {
        "integral": integral, "passo": passo, "pesos": pesos,
        "produtos": produtos, "soma_ponderada": soma_ponderada,
        "pontos": pontos, "valores": valores,
    }, partes


def _amostrar_funcao(funcao, inicio, fim, subintervalos):
    if not isinstance(subintervalos, int) or subintervalos < 1:
        raise ValueError("O número de subintervalos deve ser um inteiro positivo.")
    if not isfinite(inicio) or not isfinite(fim) or fim <= inicio:
        raise ValueError(
            "Os limites devem ser finitos e o limite final deve ser "
            "maior que o inicial."
        )
    passo = (fim - inicio) / subintervalos
    pontos = [inicio + i * passo for i in range(subintervalos)] + [fim]
    try:
        valores = [funcao(ponto) for ponto in pontos]
        distancia = (fim - inicio) * PASSO_SONDAGEM
        for indice, extremo, sentido in ((0, inicio, 1), (-1, fim, -1)):
            adjacente = valores[1] if indice == 0 else valores[-2]
            if abs(valores[indice]) <= RAZAO_SUSPEITA_POLO * abs(adjacente):
                continue
            vizinhos = [abs(funcao(extremo + sentido * distancia * i))
                        for i in (1, 2)]
            if (abs(valores[indice]) > RAZAO_SUSPEITA_POLO * max(vizinhos)
                    and vizinhos[0] > RAZAO_CRESCIMENTO_POLO * vizinhos[1]):
                raise ValueError("Possível polo ou descontinuidade no extremo.")
    except (ValueError, ZeroDivisionError, OverflowError) as erro:
        raise ValueError(
            "A função não pôde ser amostrada como integral regular. "
            "Verifique o domínio ou informe impropria='inicio' ou 'fim'."
        ) from erro
    return pontos, valores


def _dados_integracao(dados):
    extras = dados.keys() - (CAMPOS_TABELA | CAMPOS_FUNCAO | CAMPOS_PERFIL)
    if extras:
        raise TypeError(f"Campos de integração desconhecidos: {sorted(extras)}.")
    campos = {chave for chave, valor in dados.items() if valor is not None}
    if campos == CAMPOS_PERFIL:
        passo, alturas = float(dados["passo"]), list(dados["alturas"])
        if not isfinite(passo) or passo <= 0:
            raise ValueError("O passo deve ser positivo e finito.")
        return [indice * passo for indice in range(len(alturas))], alturas
    if campos == CAMPOS_FUNCAO:
        if len(dados["intervalo"]) != 2:
            raise ValueError("Informe intervalo=(inicio, fim) e subintervalos.")
        return _amostrar_funcao(
            dados["funcao"], *dados["intervalo"], dados["subintervalos"]
        )
    if campos == CAMPOS_TABELA:
        return dados["pontos"], dados["valores"]
    raise ValueError(
        "Use pontos/valores, funcao/intervalo/subintervalos ou alturas/passo."
    )


def _integrar(metodo, dados, titulo):
    if "impropria" in dados:
        resultado, partes = _integrar_impropria(metodo, dados, titulo)
        return _publicar(resultado, partes)
    dados = dict(dados)
    unidade = dados.pop("unidade", "u.a.")
    area_base = float(dados.pop("area_base", 0))
    perfis, trechos = dados.pop("perfis", None), dados.pop("trechos", None)
    extras = dados.keys() - (CAMPOS_TABELA | CAMPOS_FUNCAO)
    if extras:
        raise TypeError(f"Campos de integração desconhecidos: {sorted(extras)}.")
    if not isfinite(area_base):
        raise ValueError("A área base deve ser finita.")
    if perfis is not None and trechos is not None:
        raise ValueError("Escolha perfis ou trechos, sem misturar os formatos.")
    chave_resultado = "perfis" if perfis is not None else "trechos"
    nome_parte = "Perfil" if perfis is not None else "Trecho"
    trechos = perfis if perfis is not None else trechos
    if trechos is None:
        if area_base != 0:
            raise ValueError("Informe perfis ou trechos para usar uma área base.")
        pontos, valores = _dados_integracao(dados)
        resultado, partes = _integrar_dados(pontos, valores, metodo, unidade, titulo)
        return _publicar(resultado, partes)
    if any(valor is not None for valor in dados.values()):
        raise ValueError("Ao usar perfis ou trechos, coloque os dados dentro deles.")
    trechos = list(trechos)
    if not trechos:
        raise ValueError("Informe pelo menos um perfil ou trecho.")
    resultados, partes = [], [
        f"# {titulo}", f"Área base = {_numero(area_base)} {unidade}.",
    ]
    for indice, trecho in enumerate(trechos, start=1):
        dados_trecho = dict(trecho)
        sinal = dados_trecho.pop("sinal", 1)
        if sinal not in (-1, 1):
            raise ValueError("O sinal deve ser 1 para somar ou -1 para descontar.")
        pontos, valores = _dados_integracao(dados_trecho)
        resultado, detalhes = _integrar_dados(
            pontos, valores, metodo, unidade, f"{nome_parte} {indice}"
        )
        detalhes[0] = f"## {nome_parte} {indice}"
        contribuicao = sinal * resultado["integral"]
        resultado.update({"sinal": sinal, "contribuicao": contribuicao})
        detalhes.append(
            f"Contribuição = {sinal} * {_numero(resultado['integral'])} "
            f"= {_numero(contribuicao)} {unidade}."
        )
        resultados.append(resultado)
        partes.extend(detalhes)
    integral = area_base + sum(r["contribuicao"] for r in resultados)
    if not isfinite(integral):
        raise ArithmeticError("A soma das integrais excedeu a faixa numérica.")
    parcelas = " + ".join(f"({_numero(r['contribuicao'])})" for r in resultados)
    partes.append(
        f"Integral total = {_numero(area_base)} + {parcelas} "
        f"= **{_numero(integral)} {unidade}**."
    )
    return _publicar({
        "integral": integral, "area_base": area_base,
        chave_resultado: resultados,
    }, partes)


def _integrar_com_precisao(
    funcao, inicio, fim, casas_decimais=4, metodo="simpson",
    max_subintervalos=1000, unidade="u.a.",
    erro_maximo=None, limite_derivada=None,
    subintervalos_iniciais=None, publicar=True,
):
    if not isinstance(casas_decimais, int) or casas_decimais < 0:
        raise ValueError("As casas decimais devem ser um inteiro não negativo.")
    incremento = 2 if metodo == "simpson" else 1
    minimo = incremento if limite_derivada is not None else 2 * incremento
    if not isinstance(max_subintervalos, int) or max_subintervalos < minimo:
        raise ValueError(
            f"O limite de subintervalos deve ser um inteiro >= {minimo}."
        )
    tolerancia = 0.5 * 10 ** (-casas_decimais) if erro_maximo is None else erro_maximo
    if not isfinite(tolerancia) or tolerancia <= 0:
        raise ValueError("O erro máximo deve ser positivo e finito.")
    if not isfinite(inicio) or not isfinite(fim) or fim <= inicio:
        raise ValueError("Informe limites finitos, com o final maior que o inicial.")
    if limite_derivada is not None and (
        not isfinite(limite_derivada) or limite_derivada < 0
    ):
        raise ValueError("O limite da derivada deve ser finito e não negativo.")
    ordem = 4 if metodo == "simpson" else 2
    maximo = max_subintervalos - max_subintervalos % incremento
    quantidade = (
        incremento if subintervalos_iniciais is None else subintervalos_iniciais
    )
    if (not isinstance(quantidade, int) or quantidade < incremento
            or quantidade % incremento or quantidade > maximo):
        raise ValueError("Quantidade inicial de subintervalos inválida para o método.")
    coeficiente = None
    if limite_derivada is not None:
        divisor = 180 if metodo == "simpson" else 12
        coeficiente = 0.0 if limite_derivada == 0 else (
            (fim - inicio) ** (ordem + 1) * limite_derivada / divisor
        )
        raiz = (coeficiente / tolerancia) ** (1 / ordem)
        quantidade = maximo if not isfinite(raiz) else min(
            maximo, max(incremento, ceil(raiz))
        )
        quantidade += -quantidade % incremento
        while quantidade > incremento and (
            coeficiente / (quantidade - incremento) ** ordem <= tolerancia
        ):
            quantidade -= incremento
        while coeficiente / quantidade ** ordem > tolerancia and quantidade < maximo:
            quantidade += incremento
        quantidade = min(quantidade, maximo)
    historico, anterior, convergiu = [], None, False
    quantidade_anterior = None
    while True:
        pontos, valores = _amostrar_funcao(funcao, inicio, fim, quantidade)
        resultado, detalhes = _integrar_dados(
            pontos, valores, metodo, unidade, "Última tabela"
        )
        integral = resultado["integral"]
        diferenca = None if anterior is None else abs(integral - anterior)
        erro_estimado = None
        if coeficiente is not None:
            erro_estimado = coeficiente / quantidade ** ordem
        elif anterior is not None:
            erro_estimado = diferenca / (
                (quantidade / quantidade_anterior) ** ordem - 1
            )
        historico.append([
            quantidade, _numero(integral),
            "-" if diferenca is None else _numero(diferenca),
        ])
        alvo = tolerancia if coeficiente is not None else tolerancia / 2
        if erro_estimado is not None and erro_estimado <= alvo:
            convergiu = True
            break
        if coeficiente is not None or quantidade == maximo:
            break
        anterior = integral
        quantidade_anterior = quantidade
        quantidade = min(2 * quantidade, maximo)
    partes = [
        "# Integração por comparação de resultados", f"Método: {metodo}; "
        f"erro solicitado <= {_numero(tolerancia)}.",
        _tabela(["n", "Integral", "Diferença absoluta"], historico),
        "Critério satisfeito." if convergiu else
        "**Limite atingido sem satisfazer o critério.**",
        (f"Cota de erro = {_numero(erro_estimado)}; usa o limite informado "
         f"para |f^({ordem})|. A garantia depende dessa cota ser válida "
         "em todo o intervalo e da regularidade da função. "
         "Quando o critério é satisfeito, n é o menor permitido pela cota. "
         "A cota de truncamento não inclui o arredondamento de float.")
        if coeficiente is not None else
        "Estimativa de Richardson entre malhas dobradas, com margem de 2, "
        "limitada por max_subintervalos. "
        "É um critério empírico; não garante o erro real nem o n mínimo "
        "para uma função arbitrária. Informe limite_derivada para usar "
        "a cota teórica.", *detalhes,
    ]
    resultado.update({
        "subintervalos": quantidade, "convergiu": convergiu,
        "historico": historico,
        "erro_estimado": erro_estimado, "erro_maximo": tolerancia,
        "erro_garantido": coeficiente is not None and convergiu,
    })
    return _publicar(resultado, partes) if publicar else (resultado, partes)


def _integrar_impropria(metodo, dados, titulo):
    dados = dict(dados)
    extremidade = dados.pop("impropria")
    erro_maximo = dados.pop("erro_maximo", 1e-6)
    max_refinamentos = dados.pop("max_refinamentos", 60)
    max_subintervalos = dados.pop("max_subintervalos", 16384)
    subintervalos = dados.pop("subintervalos", 16)
    unidade = dados.pop("unidade", "u.a.")
    funcao, intervalo = dados.pop("funcao"), dados.pop("intervalo")
    if dados:
        raise TypeError(f"Campos de integração desconhecidos: {sorted(dados)}.")
    if extremidade not in ("inicio", "fim"):
        raise ValueError("impropria deve ser 'inicio' ou 'fim'.")
    if len(intervalo) != 2:
        raise ValueError("Informe intervalo=(inicio, fim).")
    inicio, fim = intervalo
    if not isfinite(inicio) or not isfinite(fim) or fim <= inicio:
        raise ValueError("Informe limites finitos, com o final maior que o inicial.")
    if not isfinite(erro_maximo) or erro_maximo <= 0:
        raise ValueError("O erro máximo deve ser positivo e finito.")
    if not isinstance(max_refinamentos, int) or max_refinamentos < 1:
        raise ValueError("max_refinamentos deve ser um inteiro positivo.")
    total, erro_quadratura = 0.0, 0.0
    historico, parcelas = [], []
    convergiu, possivel_divergencia = False, False
    cauda = None
    for indice in range(1, max_refinamentos + 1):
        distancia = (fim - inicio) / 2 ** indice
        esquerda, direita = (
            (inicio + distancia, inicio + 2 * distancia)
            if extremidade == "inicio" else
            (fim - 2 * distancia, fim - distancia)
        )
        if esquerda == direita or esquerda == inicio or direita == fim:
            if indice > 1:
                break
        largura = direita - esquerda
        normalizada = lambda t: largura * funcao(esquerda + largura * t)
        parcela, _ = _integrar_com_precisao(
            normalizada, 0, 1, metodo=metodo,
            erro_maximo=erro_maximo / (8 * indice ** 2),
            max_subintervalos=max_subintervalos,
            subintervalos_iniciais=subintervalos, publicar=False,
        )
        if not parcela["convergiu"]:
            break
        valor = parcela["integral"]
        total += valor
        erro_quadratura += parcela["erro_estimado"]
        parcelas.append(valor)
        if len(parcelas) >= 4:
            recentes = parcelas[-4:]
            if all(a * b > 0 for a, b in zip(recentes, recentes[1:])):
                razoes = [abs(b / a) for a, b in zip(recentes, recentes[1:])]
                razao = max(razoes)
                cauda = abs(valor) * razao / (1 - razao) if razao < 1 else None
                possivel_divergencia = (
                    indice >= MIN_PARCELAS_DIVERGENCIA
                    and min(razoes) >= RAZAO_SEM_DECAIMENTO
                )
            else:
                cauda = None
        historico.append([indice, distancia, valor, total, cauda])
        if cauda is not None and cauda + erro_quadratura <= erro_maximo / 2:
            convergiu = True
            break
        if possivel_divergencia:
            break
    resultado = {
        "integral": total if convergiu else None, "aproximacao": total,
        "convergiu": convergiu, "erro_garantido": False,
        "possivel_divergencia": possivel_divergencia, "historico": historico,
        "erro_estimado": None if cauda is None else cauda + erro_quadratura,
    }
    partes = [
        f"# {titulo}: integral imprópria", f"Extremidade singular: {extremidade}.",
        _tabela(["Parcela", "Distância", "Integral", "Soma", "Cauda estimada"], [
            [*linha[:-1], "-" if linha[-1] is None else linha[-1]]
            for linha in historico
        ]),
        f"Integral aproximada = **{_numero(total)} {unidade}**." if convergiu else
        "**Não foi validada uma integral finita; a soma truncada não é a resposta.**",
        "Há indícios numéricos de divergência: as parcelas não decaem." if
        possivel_divergencia else "As parcelas foram integradas sem avaliar o polo.",
        "A cauda foi estimada por razões entre parcelas. É um teste empírico, "
        "sem prova de convergência, divergência ou garantia universal do erro.",
    ]
    return resultado, partes
