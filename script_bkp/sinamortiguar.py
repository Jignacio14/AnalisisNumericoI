import utils
import euler
import rungekutta2

def sinAmortiguar():

    'Funcion sin amortiguacion'
    def f(x):
        return (utils.k/utils.m)*(utils.c-x)

    'Matriz inversa para la solucion de euler implicito sin amortiguacion'
    divisor = (utils.h**2)*utils.k+utils.h*utils.lambda_+utils.m
    A_invertida = utils.np.array([[utils.m/divisor , -utils.h*utils.k/divisor],
                        [utils.h*utils.m/divisor , (utils.h*utils.lambda_+utils.m)/divisor]])

    termino_indep = utils.np.array([utils.h*((utils.k*utils.c/utils.m)+(utils.lambda_*utils.c_prim/utils.m)),0])

    t=utils.np.arange(0,utils.intervalo+utils.h,utils.h)


    'Metodo de euler explicito sin amortiguar'
    euler.imprimirEulerExplicito(f,utils.intervalo,utils.h,utils.t)
    euler.imprimirEulerImplicito(A_invertida, termino_indep,utils.intervalo,utils.h,t)
    
    'Metodo de Runge Kutta sin amortiguar'
    t=utils.np.arange(0,utils.intervalo+utils.h,utils.h)
    rungekutta2.imprimirRungeKutta2(f,utils.intervalo,utils.h,t)
