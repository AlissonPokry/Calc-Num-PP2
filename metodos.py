from math import isfinite

from auxiliares import (
    MAX_ITERACOES,
    TOLERANCIA_ITERACAO,
    _criterios_convergencia,
    _determinante,
    _eliminar,
    _integrar,
    _integrar_com_precisao,
    _numero,
    _preparar_sistema,
    _publicar,
    _publicar_sistema,
    _relatar_criterios,
    _substituir,
    _tabela,
    _tabela_matriz
)


def substituicao_retroativa(dados):
    matriz, termos = _preparar_sistema(dados["matriz"], dados["termos"])
    if any(matriz[i][j] != 0 for i in range(len(matriz)) for j in range(i)):
        raise ValueError("A matriz deve ser triangular superior.")
    solucao, etapas = _substituir(matriz, termos)
    return _publicar({"solucao": solucao}, [
        "# Substituição retroativa", _tabela_matriz(matriz, termos), *etapas
    ])


def substituicao_direta(dados):
    matriz, termos = _preparar_sistema(dados["matriz"], dados["termos"])
    if any(matriz[i][j] != 0 for i in range(len(matriz))
           for j in range(i + 1, len(matriz))):
        raise ValueError("A matriz deve ser triangular inferior.")
    solucao, etapas = _substituir(matriz, termos, inferior=True)
    return _publicar({"solucao": solucao}, [
        "# Substituição direta", _tabela_matriz(matriz, termos), *etapas
    ])


def gauss(dados, pivotamento="nenhum", titulo="Eliminação de Gauss"):
    matriz, termos = _preparar_sistema(dados["matriz"], dados["termos"])
    superior, transformados, _, _, variaveis, etapas = _eliminar(
        matriz, termos, pivotamento
    )
    reordenada, substituicoes = _substituir(
        superior, transformados, variaveis=variaveis
    )
    solucao = [0.0] * len(termos)
    for coluna, variavel in enumerate(variaveis):
        solucao[variavel] = reordenada[coluna]
    return _publicar_sistema(
        matriz, termos, solucao,
        [f"# {titulo}", *etapas, "## Substituição retroativa", *substituicoes],
        triangular=superior, termos_transformados=transformados,
        ordem_variaveis=variaveis,
    )


def jordan(dados, titulo="Método de Gauss-Jordan"):
    matriz, termos = _preparar_sistema(dados["matriz"], dados["termos"])
    reduzida, solucao, _, _, _, etapas = _eliminar(matriz, termos, "parcial")
    for coluna in range(len(termos) - 1, -1, -1):
        pivo = reduzida[coluna][coluna]
        base = [valor / pivo for valor in reduzida[coluna]]
        reduzida[coluna] = base
        solucao[coluna] /= pivo
        etapas.extend([f"L{coluna + 1} <- L{coluna + 1} / {_numero(pivo)}.",
                       _tabela_matriz(reduzida, solucao)])
        for linha in range(coluna):
            fator = -reduzida[linha][coluna]
            if fator == 0:
                continue
            reduzida[linha] = [a + fator * b for a, b in zip(reduzida[linha], base)]
            reduzida[linha][coluna] = 0.0
            solucao[linha] += fator * solucao[coluna]
            etapas.extend([
                f"L{linha + 1} <- L{linha + 1} + "
                f"({_numero(fator)}) * L{coluna + 1}.",
                _tabela_matriz(reduzida, solucao),
            ])
    return _publicar_sistema(
        matriz, termos, solucao, [f"# {titulo}", *etapas], reduzida=reduzida
    )


def decomposicao_lu(dados, titulo="Decomposição LU"):
    matriz, termos = _preparar_sistema(dados["matriz"], dados["termos"])
    superior, _, inferior, indices, _, etapas = _eliminar(matriz, termos, "parcial")
    permutacao = [
        [float(j == indice) for j in range(len(termos))] for indice in indices
    ]
    termos_permutados = [termos[indice] for indice in indices]
    intermediario, diretas = _substituir(inferior, termos_permutados, inferior=True)
    solucao, retroativas = _substituir(superior, intermediario)
    return _publicar_sistema(matriz, termos, solucao, [
        f"# {titulo}", *etapas, "## Fatoração P*A = L*U",
        "P:", _tabela_matriz(permutacao), "L:", _tabela_matriz(inferior),
        "U:", _tabela_matriz(superior), "## Substituição direta: L*y = P*b",
        *[etapa.replace("x", "y") for etapa in diretas],
        "## Substituição retroativa: U*x = y", *retroativas,
    ], inferior=inferior, superior=superior, permutacao=permutacao,
    intermediario=intermediario)


def cramer(dados, titulo="Regra de Cramer"):
    matriz, termos = _preparar_sistema(dados["matriz"], dados["termos"])
    determinante = _determinante(matriz)
    if determinante == 0:
        raise ValueError("det(A) é nulo numericamente. Use classificar_sistema().")
    if not isfinite(determinante):
        raise ArithmeticError("O determinante excedeu a faixa numérica; prefira Gauss.")
    solucao, determinantes = [], []
    partes = [
        f"# {titulo}", _tabela_matriz(matriz, termos),
        f"det(A) = {_numero(determinante)}.",
    ]
    for coluna in range(len(termos)):
        substituida = [linha.copy() for linha in matriz]
        for linha, termo in zip(substituida, termos):
            linha[coluna] = termo
        determinante_coluna = _determinante(substituida)
        valor = determinante_coluna / determinante
        if not isfinite(valor):
            raise ArithmeticError(
                "A razão de determinantes não é finita; prefira Gauss."
            )
        determinantes.append(determinante_coluna)
        solucao.append(valor)
        partes.extend([f"## D{coluna + 1}", _tabela_matriz(substituida),
                       f"x{coluna + 1} = det(D{coluna + 1}) / det(A) = "
                       f"{_numero(determinante_coluna)} / "
                       f"{_numero(determinante)} = {_numero(valor)}."])
    return _publicar_sistema(
        matriz, termos, solucao, partes, determinante=determinante,
        determinantes_colunas=determinantes,
    )


def gauss_seidel(
    dados, inicial=None, tolerancia=TOLERANCIA_ITERACAO,
    max_iteracoes=MAX_ITERACOES, casas_decimais=None, titulo="Gauss-Seidel",
    *, erro_relativo=False,
):
    matriz, termos = _preparar_sistema(dados["matriz"], dados["termos"])
    if not isfinite(tolerancia) or tolerancia <= 0:
        raise ValueError("A tolerância deve ser positiva e finita.")
    if not isinstance(max_iteracoes, int) or max_iteracoes < 1:
        raise ValueError("O número máximo de iterações deve ser um inteiro positivo.")
    if not isinstance(erro_relativo, bool):
        raise ValueError("erro_relativo deve ser True ou False.")
    if casas_decimais is not None and (
        not isinstance(casas_decimais, int) or casas_decimais < 0
    ):
        raise ValueError("As casas decimais devem ser um inteiro não negativo ou None.")
    solucao = [0.0] * len(termos) if inicial is None else list(map(float, inicial))
    if len(solucao) != len(termos) or not all(map(isfinite, solucao)):
        raise ValueError("O vetor inicial deve conter um número finito por incógnita.")
    criterios = _criterios_convergencia(matriz)
    partes = [
        f"# {titulo}", _tabela_matriz(matriz, termos),
        *_relatar_criterios(criterios), "## Equações de iteração",
    ]
    for i, linha in enumerate(matriz):
        parcelas = " + ".join(f"({_numero(valor)})*x{j + 1}({'k+1' if j < i else 'k'})"
                              for j, valor in enumerate(linha) if j != i) or "0"
        partes.append(
            f"x{i + 1}(k+1) = ({_numero(termos[i])} - ({parcelas})) "
            f"/ {_numero(linha[i])}"
        )
    arredondamento = "sem arredondar" if casas_decimais is None else casas_decimais
    partes.append(
        f"Tolerância = {_numero(tolerancia)}; "
        f"máximo de iterações = {max_iteracoes}; "
        f"arredondamento durante as iterações = {arredondamento}."
    )
    if erro_relativo:
        partes.append("Distância relativa = max|x(k)-x(k-1)| / max|x(k)|.")
    historico = [{"iteracao": 0, "solucao": solucao.copy(), "diferenca": None}]
    if erro_relativo:
        historico[0]["diferenca_relativa"] = None
    linhas_iteracoes = [[0, *solucao, "-", "Continua"]]
    convergiu = False
    for iteracao in range(1, max_iteracoes + 1):
        anterior = solucao.copy()
        for i, linha in enumerate(matriz):
            soma = sum(valor * solucao[j] for j, valor in enumerate(linha) if j != i)
            solucao[i] = (termos[i] - soma) / linha[i]
            if casas_decimais is not None:
                solucao[i] = round(solucao[i], casas_decimais)
        if not all(map(isfinite, solucao)):
            raise ArithmeticError(
                "A iteração produziu valores não finitos; "
                "verifique convergência e escala."
            )
        diferenca = max(abs(novo - antigo) for novo, antigo in zip(solucao, anterior))
        historico.append({
            "iteracao": iteracao, "solucao": solucao.copy(),
            "diferenca": diferenca,
        })
        diferenca_parada = diferenca
        if erro_relativo:
            norma = max(map(abs, solucao))
            relativa = diferenca / norma if norma else (
                float("inf") if diferenca else 0
            )
            historico[-1]["diferenca_relativa"] = relativa
            diferenca_parada = relativa
        convergiu = diferenca_parada < tolerancia
        linhas_iteracoes.append([
            iteracao, *solucao, diferenca_parada, "Para" if convergiu else "Continua",
        ])
        if convergiu:
            break
    partes.extend([
        "## Iterações", _tabela(
            ["k", *[f"x{i + 1}" for i in range(len(termos))],
             "max_dif_rel" if erro_relativo else "max_dif", "Parada"],
            linhas_iteracoes,
        ),
        "Critério de parada satisfeito." if convergiu else
        "**Limite de iterações atingido sem satisfazer a tolerância. "
        "A aproximação não foi aceita.**",
    ])
    return _publicar_sistema(
        matriz, termos, solucao, partes, convergiu=convergiu,
        iteracoes=iteracao, historico=historico, criterios=criterios,
        rodape=[
            "O teste entre iterações não limita sozinho o erro da solução; "
            "confira também o resíduo."
        ],
    )


def trapezios(dados, titulo="Regra dos trapézios composta"):
    return _integrar("trapezios", dados, titulo)


def simpson(dados, titulo="Regra de Simpson 1/3 composta"):
    return _integrar("simpson", dados, titulo)


def trapezios_com_precisao(
    funcao, intervalo, casas_decimais=4, max_subintervalos=1000, unidade="u.a.",
    *, erro_maximo=None, limite_derivada=None,
):
    return _integrar_com_precisao(
        funcao, *intervalo, casas_decimais=casas_decimais,
        metodo="trapezios", max_subintervalos=max_subintervalos, unidade=unidade,
        erro_maximo=erro_maximo, limite_derivada=limite_derivada,
    )


def simpson_com_precisao(
    funcao, intervalo, casas_decimais=4, max_subintervalos=1000, unidade="u.a.",
    *, erro_maximo=None, limite_derivada=None,
):
    return _integrar_com_precisao(
        funcao, *intervalo, casas_decimais=casas_decimais,
        metodo="simpson", max_subintervalos=max_subintervalos, unidade=unidade,
        erro_maximo=erro_maximo, limite_derivada=limite_derivada,
    )
