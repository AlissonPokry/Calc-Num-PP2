import contextlib
import io
from math import *
from pathlib import Path
import sys
import unittest
from tempfile import TemporaryDirectory

PASTA = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PASTA))
import auxiliares
import main
import metodos


class CorrecaoNumerica(unittest.TestCase):
    def setUp(self):
        self.temporaria = TemporaryDirectory(dir=Path(__file__).resolve().parent)
        self.arquivo_anterior = auxiliares.ARQUIVO_RESULTADO
        auxiliares.ARQUIVO_RESULTADO = Path(self.temporaria.name) / 'resultado.md'

    def tearDown(self):
        auxiliares.ARQUIVO_RESULTADO = self.arquivo_anterior
        self.temporaria.cleanup()

    def chamar(self, metodo, *args, **kwargs):
        terminal = io.StringIO()
        with contextlib.redirect_stdout(terminal):
            resultado = metodo(*args, **kwargs)
        salvo = auxiliares.ARQUIVO_RESULTADO.read_text(encoding='utf-8')
        self.assertEqual(salvo, resultado['relatorio'])
        self.assertEqual(terminal.getvalue(), salvo+'\n')
        return resultado

    def test_resultados_anteriores_dos_exercicios(self):
        for dados, esperado in (
            (main.SISTEMA_1,[1,-1,3]),(main.SISTEMA_2,[1,-2,1]),
            (main.SISTEMA_3,[1,1,-1]),
            (main.SISTEMA_4,[-13/14,-285/28,-251/28,-100/7]),
        ):
            for metodo in (metodos.gauss,metodos.jordan,metodos.decomposicao_lu,metodos.cramer):
                with self.subTest(metodo=metodo.__name__, esperado=esperado):
                    r=self.chamar(metodo,dados)
                    for valor, exato in zip(r['solucao'],esperado):
                        self.assertAlmostEqual(valor,exato,places=12)
        for metodo,dados,esperado in (
            (metodos.trapezios,main.LAGOA,570),
            (metodos.simpson,main.LAGOA,1784/3),
            (metodos.simpson,main.TERRENO,4192.5),
        ):
            self.assertAlmostEqual(self.chamar(metodo,dados)['integral'],esperado,places=12)

    def test_erro_preserva_relatorio_anterior(self):
        self.chamar(metodos.gauss,main.SISTEMA_1)
        anterior=auxiliares.ARQUIVO_RESULTADO.read_bytes()
        with self.assertRaises(ValueError):
            self.chamar(metodos.simpson,{'funcao':tan,'intervalo':(1,pi/2),'subintervalos':10})
        self.assertEqual(auxiliares.ARQUIVO_RESULTADO.read_bytes(),anterior)

    def test_refinamento_atende_erro_das_quatro_integrais(self):
        for funcao, intervalo, exata in (
            (lambda x: exp(-x), (1, 2), exp(-1)-exp(-2)),
            (sqrt, (1, 4), 14/3),
            (lambda x: 1/sqrt(x), (2, 5), 2*(sqrt(5)-sqrt(2))),
            (lambda x: 1/(1+x), (0, .6), log(1.6)),
            (lambda x: sqrt(sin(x)), (1, pi/2), .555162600078254),
        ):
            for metodo in (metodos.trapezios_com_precisao, metodos.simpson_com_precisao):
                with self.subTest(metodo=metodo.__name__, intervalo=intervalo):
                    r = self.chamar(metodo, funcao, intervalo, casas_decimais=5)
                    self.assertTrue(r['convergiu'])
                    self.assertLess(abs(r['integral']-exata), 5e-6)

    def test_minimo_pela_cota_de_derivadas(self):
        for funcao, intervalo, m2, m4, nt, ns in (
            (lambda x: exp(-x), (1,2), exp(-1), exp(-1), 56, 4),
            (sqrt, (1,4), 1/4, 15/16, 238, 20),
            (lambda x: 1/sqrt(x), (2,5), 3/(4*2**2.5), 105/(16*2**4.5), 173, 16),
            (lambda x: 1/(1+x), (0,.6), 2, 24, 60, 6),
        ):
            for metodo, limite, esperado in (
                (metodos.trapezios_com_precisao,m2,nt),
                (metodos.simpson_com_precisao,m4,ns),
            ):
                with self.subTest(metodo=metodo.__name__, intervalo=intervalo):
                    r = self.chamar(metodo,funcao,intervalo,erro_maximo=1e-5,limite_derivada=limite)
                    self.assertEqual(r['subintervalos'], esperado)
                    self.assertTrue(r['convergiu'])
                    self.assertTrue(r['erro_garantido'])
                    self.assertLessEqual(r['erro_estimado'], 1e-5)

    def test_cota_zero_aceita_polinomio_exato(self):
        r = self.chamar(metodos.trapezios_com_precisao,lambda x:2*x+1,(0,2),limite_derivada=0)
        self.assertEqual(r['subintervalos'],1)
        self.assertEqual(r['integral'],6)
        self.assertTrue(r['erro_garantido'])

    def test_limite_nao_confirma_precisao(self):
        r = self.chamar(metodos.trapezios_com_precisao,exp,(0,1),erro_maximo=1e-12,limite_derivada=e,max_subintervalos=8)
        self.assertFalse(r['convergiu'])
        self.assertFalse(r['erro_garantido'])

    def test_erro_relativo_seidel(self):
        dados={'matriz':[[3,-2,5,1],[-6,4,-8,1],[9,-6,19,1],[6,-4,-6,15]],'termos':[7,-9,23,11]}
        r=self.chamar(metodos.gauss_seidel,dados,max_iteracoes=2,tolerancia=.001,erro_relativo=True)
        self.assertAlmostEqual(r['historico'][1]['diferenca_relativa'],1)
        self.assertAlmostEqual(r['historico'][2]['diferenca_relativa'],.3375)
        self.assertFalse(r['convergiu'])

    def test_erro_relativo_zero(self):
        r=self.chamar(metodos.gauss_seidel,{'matriz':[[1]],'termos':[0]},erro_relativo=True)
        self.assertTrue(r['convergiu'])
        self.assertEqual(r['historico'][1]['diferenca_relativa'],0)
        r=self.chamar(metodos.gauss_seidel,{'matriz':[[1]],'termos':[0]},inicial=[1],erro_relativo=True)
        self.assertEqual(r['historico'][1]['diferenca_relativa'],inf)
        self.assertEqual(r['iteracoes'],2)

    def test_validacoes_precisao(self):
        for opcoes in ({'erro_maximo':0},{'erro_maximo':float('nan')},{'limite_derivada':-1},{'limite_derivada':float('inf')}):
            with self.subTest(opcoes=opcoes):
                with self.assertRaises(ValueError):
                    self.chamar(metodos.simpson_com_precisao,exp,(0,1),**opcoes)

    def test_integral_impropria_convergente(self):
        for metodo in (metodos.trapezios, metodos.simpson):
            for n in (10, 20):
                with self.subTest(metodo=metodo.__name__, n=n):
                    r=self.chamar(metodo,{
                        'funcao':lambda x:6/sqrt(3*x),'intervalo':(0,9),
                        'subintervalos':n,'impropria':'inicio','erro_maximo':1e-6,
                    })
                    self.assertTrue(r['convergiu'])
                    self.assertLess(abs(r['integral']-12*sqrt(3)),1e-6)
                    self.assertFalse(r['erro_garantido'])

    def test_polo_nao_produz_integral_finita(self):
        for metodo in (metodos.trapezios,metodos.simpson):
            with self.subTest(metodo=metodo.__name__):
                with self.assertRaises(ValueError):
                    self.chamar(metodo,{'funcao':tan,'intervalo':(1,pi/2),'subintervalos':10})
                r=self.chamar(metodo,{
                    'funcao':tan,'intervalo':(1,pi/2),'impropria':'fim',
                })
                self.assertIsNone(r['integral'])
                self.assertFalse(r['convergiu'])
                self.assertTrue(r['possivel_divergencia'])

    def test_modo_regular_exige_funcao_definida(self):
        with self.assertRaises(ValueError):
            self.chamar(metodos.simpson,{'funcao':lambda x:1/sqrt(x),'intervalo':(0,1),'subintervalos':10})

    def test_refinamento_limite_nao_dobrado(self):
        r=self.chamar(metodos.trapezios_com_precisao,lambda x:exp(-x),(1,2),casas_decimais=9,max_subintervalos=513)
        if r['convergiu']:
            self.assertLess(abs(r['integral']-(exp(-1)-exp(-2))),5e-10)

    def test_quadratica_regular_com_raizes_nas_sondagens(self):
        a,b=1-1e-6,1-2e-6
        r=self.chamar(metodos.simpson,{'funcao':lambda x:(x-a)*(x-b),'intervalo':(0,1),'subintervalos':10})
        self.assertAlmostEqual(r['integral'],1/3-(a+b)/2+a*b,places=14)

    def test_cascas_nulas_nao_confirmam_cauda(self):
        for expoente in (-.5,-1):
            with self.subTest(expoente=expoente):
                r=self.chamar(metodos.simpson,{'funcao':lambda x:x**expoente if x<.01 else 0,'intervalo':(0,1),'impropria':'inicio'})
                self.assertFalse(r['convergiu'] and r['integral']==0)

    def test_convergencia_lenta_nao_e_divergencia(self):
        r=self.chamar(metodos.simpson,{'funcao':lambda x:x**(-.999),'intervalo':(0,1),'impropria':'inicio','max_refinamentos':24})
        self.assertFalse(r['convergiu'])
        self.assertFalse(r['possivel_divergencia'])

    def test_validacoes_impropria(self):
        for opcoes in ({'impropria':'outro'},{'erro_maximo':0},{'max_refinamentos':0},{'subintervalos':3}):
            with self.subTest(opcoes=opcoes):
                dados={'funcao':lambda x:1/sqrt(x),'intervalo':(0,1),'impropria':'inicio'}
                dados.update(opcoes)
                with self.assertRaises(ValueError):
                    self.chamar(metodos.simpson,dados)


if __name__ == '__main__':
    unittest.main(verbosity=2)
