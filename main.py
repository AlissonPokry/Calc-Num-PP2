from math import *
from metodos import *

SISTEMA_1 = {
    "matriz": [ [3, -2, 5], 
                [6, -9, 12], 
                [-5, 0, 2]],
    "termos": [20, 51, 1],
}

SISTEMA_2 = {
    "matriz": [ [10, 2, 1], 
                [1, 5, 1], 
                [2, 3, 10]],
    "termos": [7, -8, 6],
}

SISTEMA_3 = {
    "matriz": [ [5, 1, 1], 
                [3, 4, 1], 
                [3, 3, 6]],
    "termos": [5, 6, 0],
}

SISTEMA_4 = {
    "matriz": [
        [2, -1, -1, 1],
        [3, 4, -2, -2],
        [3, -2, 4, -1],
        [2, -4, 0, 3],
    ],
    "termos": [3, 3, -4, -4],
}

SISTEMA_TRIANGULAR_INFERIOR = {
    "matriz": [
        [2, 0, 0],
        [1, 3, 0],
        [4, -1, 2],
    ],
    "termos": [2, 7, 8],
}

SISTEMA_TRIANGULAR_SUPERIOR = {
    "matriz": [
        [2, -1, 3],
        [0, 4, 2],
        [0, 0, 5],
    ],
    "termos": [9, 14, 15],
}

LAGOA = {
    "perfis": [
        {
            "alturas": [0, 3, 6, 9, 10, 9, 8, 6, 0],
            "passo": 6,
        },
        {
            "alturas": [0, 4, 5, 9, 8, 7, 0],
            "passo": 8,
        },
    ],
    "unidade": "km²",
}

ALTURA_TERRENO = 65
LARGURA_ESQUERDA = 40
BASE_MAIOR_DIREITA = 70
BASE_MENOR_DIREITA = 15
TERRENO = {
    "area_base": (
        LARGURA_ESQUERDA * ALTURA_TERRENO
        + (BASE_MAIOR_DIREITA + BASE_MENOR_DIREITA) * ALTURA_TERRENO / 2
    ),
    "perfis": [
        {
            "alturas": [0, 20, 12.5, 9, 0],
            "passo": 10,
            "sinal": -1,
        },
        {
            "alturas": [0, 21, 22, 20.5, 0],
            "passo": 10,
            "sinal": -1,
        },
    ],
    "unidade": "km²",
}

INTEGRAL = {
    "funcao": lambda x: exp(-cos(x)) / sqrt(2 * x + 4),
    "intervalo": (0, 2),
    "subintervalos": 5,
}

if __name__ == "__main__":
    gauss(SISTEMA_1)
    gauss(SISTEMA_2)
    gauss(SISTEMA_3)
    gauss(SISTEMA_4)
    jordan(SISTEMA_2)
    decomposicao_lu(SISTEMA_2)
    cramer(SISTEMA_2)
    gauss_seidel(SISTEMA_2, tolerancia=1e-8, max_iteracoes=100)
    substituicao_direta(SISTEMA_TRIANGULAR_INFERIOR)
    substituicao_retroativa(SISTEMA_TRIANGULAR_SUPERIOR)
    trapezios_com_precisao(
        INTEGRAL["funcao"], INTEGRAL["intervalo"], casas_decimais=4
    )
    simpson_com_precisao(
        INTEGRAL["funcao"], INTEGRAL["intervalo"], casas_decimais=4
    )
    trapezios(INTEGRAL)
    trapezios(LAGOA)
    simpson(LAGOA)
    simpson(TERRENO)
