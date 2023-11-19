import utils
import euler

def sinAmortiguar():

    'Funcion sin amortiguacion'
    def f(x):
        return (utils.k/utils.m)*(utils.c-x)

    'Solucion analitica sin amortiguacion'
    def analitica(t):
        return utils.c - utils.c * utils.np.cos(((k/m)**0.5)*t)

    'Matriz inversa para la solucion de euler implicito sin amortiguacion'
    divisor = (utils.h**2)*utils.k+utils.h*utils.lambda_+utils.m
    A_menos1 = utils.np.array([[utils.m/divisor , -utils.h*utils.k/divisor],
                        [utils.h*utils.m/divisor , (utils.h*utils.lambda_+utils.m)/divisor]])

    termino_indep = utils.np.array([utils.h*((utils.k*utils.c/utils.m)+(utils.lambda_*c_prim/utils.m)),0])

    t=utils.np.arange(0,utils.intervalo+utils.h,utils.h)


    'Metodo de euler explicito sin amortiguar'
    euler.imprimirEulerExplicito(f,utils.intervalo,utils.h,analitica,utils.t)
    euler.imprimirEulerImplicito()
