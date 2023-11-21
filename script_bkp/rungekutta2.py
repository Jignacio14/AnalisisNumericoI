import matplotlib.pyplot as plt
import numpy as np
import utils

def runge_kutta_segundo_orden(f, intervalo, paso):
    # Calculamos la cantidad de puntos en el intervalo
    cantidad = int(intervalo/paso) + 1

    # Creamos arrays para almacenar los resultados de u y v
    u = np.zeros(cantidad, dtype=float)
    v = np.zeros(cantidad, dtype=float)

    # Condiciones iniciales
    u[0] = 0
    v[0] = 0

    n = 0  # Inicializamos el índice del bucle
    i = 0  # Inicializamos el tiempo

    while n < cantidad-1:
        # Predictor
        u_actual_predictor = u[n] + paso * v[n]
        v_actual_predictor = v[n] + paso * f(u[n])

        # Corrector
        u[n+1] = u[n] + (paso/2) * (v[n] + v_actual_predictor)
        v[n+1] = v[n] + (paso/2) * (f(u[n]) + f(u_actual_predictor))

        n += 1
        i += paso

    return u

def armarGraficoRK2OK(aproximacion, t, h, titulo, nombre):
    plt.plot(t, aproximacion)
    plt.xlabel('t')
    plt.ylabel('y')
    plt.title(titulo + str(h))
    plt.plot(t, utils.analitica(t), 'r--')
    name = nombre + str(h) + '.png'
    plt.savefig(name)
    plt.show()
    return

def armarGraficoRK2Error(t, aproximacion, h, titulo, nombre):
    plt.plot(t, aproximacion)
    plt.xlabel('t')
    plt.ylabel('error')
    plt.title(titulo + str(h))
    name = nombre + str(h) + '.png'
    plt.savefig(name)
    plt.show()
    return

def iterarRK2(aproximacion, intervalo, h):
    paso = 0
    for j in range(int((intervalo/h) + 1)):
        aproximacion[j] = (utils.analitica(paso) - aproximacion[j])
        paso += h
    return

def imprimirRungeKutta2(f, intervalo, h, t):
    aproximacion = imprimirRungeKutta2Extras(f, intervalo, h)
    armarGraficoRK2OK(aproximacion, t, h, utils.y_t_RK2, utils.y_nombre_RK2)
    iterarRK2(aproximacion, intervalo, h)
    armarGraficoRK2Error(aproximacion, t, h, utils.e_t_RK2, utils.e_nombre_RK2)
    return

def rungeKuttaO2ConExtras(intervalo, paso, extras):
    cantidad = intervalo/paso + 1 
    u = np.zeros(int(cantidad), dtype=float)
    v = np.zeros(int(cantidad), dtype=float)
    u[0] = 0
    v[0] = 0
    i = 0
    n = 0
    while n < cantidad-1:
        u_actual_predictor = (u[(n)] + paso * v[n])
        v_actual_predictor = (v[n] + paso * f2(u[n], extras[0], i, v[n], extras[1]))
        u[(n+1)] = (u[n] + (paso/2) * (v[n] + v_actual_predictor))
        v[n+1] = (v[n] + (paso/2) * (f2(u[n], extras[0], i, v[n], extras[1]) + f2(u_actual_predictor, extras[0], i+paso, v_actual_predictor, extras[1])))
        n += 1
        i += paso
    return u, v


def f2(y,c,t,y_prim,c_prim):
    lambda_ = 750
    res =  (utils.k/utils.m)*(c(t)-y) + (lambda_/utils.m) * (c_prim(t) - y_prim)
    return res

def imprimirAmortiguado(t,aproximada,k,h):
    plt.plot(t, aproximada)
    plt.xlabel('t')
    plt.ylabel('y')
    plt.title('y(t) RK2 k = ' + str(k))
    name = 'RK2Amortiguado' + str(h) + '.png'
    plt.savefig(name)
    plt.show()
    return

def imprimirCompresion(t,aproximada,h):
    plt.plot(t, aproximada)
    plt.xlabel('t')
    plt.ylabel('compresion')
    plt.title('compresion(t) del amortiguador ')
    name = 'Compresion' + str(h) + '.png'
    plt.savefig(name)
    plt.show()
    return

def iterarCompresion(intervalo,h,aproximada,k,lam,c):
    paso = 0
    minimo = 0
    for j in range(int((intervalo/h) + 1)):
        aproximada[j] = (aproximada[j] - float(c(paso)))
        if aproximada[j] < minimo:
            minimo = aproximada[j]
        paso += h
    print('maxima compresion k = ' + str(k) + ' lambda = ' + str(lam) + '  :  ' + str(minimo))
def imprimirRungeKutta2Extras(vfun, intervalo, h, extras, t, c, k, lam):
    x = np.arange(0, intervalo+h, h)
    y = vfun(x)
    plt.plot(x, y, 'r')
    aproximada, v = rungeKuttaO2ConExtras(intervalo, h, extras)
    u = np.copy(aproximada)
    imprimirAmortiguado(t,aproximada,k,h)
    iterarCompresion(intervalo,h,aproximada,k,lam,c)
    imprimirCompresion(t,aproximada,h)
    return u, v
