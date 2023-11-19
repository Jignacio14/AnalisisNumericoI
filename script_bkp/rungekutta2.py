import matplotlib.pyplot as plt
import numpy as np
import utils

# Método de Runge-Kutta de segundo orden
def runge_kutta_segundo_orden(y0, v0, m, k, lambda_, c, tiempo_total, paso_tiempo):
    cantidad_pasos = int(tiempo_total / paso_tiempo) + 1

    y = np.zeros(cantidad_pasos)
    v = np.zeros(cantidad_pasos)

    y[0] = y0
    v[0] = v0

    for n in range(1, cantidad_pasos):
        # Primer paso RK2
        k1_y = paso_tiempo * v[n-1]
        k1_v = paso_tiempo * ((-k / m) * (y[n-1] - c) - (lambda_ / m) * v[n-1])

        # Segundo paso RK2
        k2_y = paso_tiempo * (v[n-1] + k1_v)
        k2_v = paso_tiempo * ((-k / m) * (y[n-1] + k1_y - c) - (lambda_ / m) * (v[n-1] + k1_v))

        # Actualizar soluciones utilizando RK2
        y[n] = y[n-1] + 0.5 * (k1_y + k2_y)
        v[n] = v[n-1] + 0.5 * (k1_v + k2_v)

    return y, v

def runge_kuta_2(f,intervalo, paso):
    cantidad = intervalo/paso +1 
    u=np.zeros(int(cantidad),dtype=float)
    v=np.zeros(int(cantidad),dtype=float)
    u[0]=0
    v[0]=0
    i=0
    n=0
    while n<cantidad-1:
        u_actual_predictor = (u[(n)] + paso * v[n])
        v_actual_predictor = (v[n] + paso * f(u[n]))
        u[(n+1)]= (u[n] + (paso/2) *(v[n] + v_actual_predictor) )
        v[n+1] = (v[n] + (paso/2) * (f(u[n]) + f(u_actual_predictor)) )
        n+=1
        i+=paso
    return u

def armarGraficoRK2OK(aproximacion,t,h,titulo,nombre):
    plt.plot(t,aproximacion)
    plt.xlabel('t')
    plt.ylabel('y')
    plt.title(titulo + str(h))
    plt.plot(t,utils.analitica(t) , 'r--')
    name = nombre+str(h)+'.png'
    plt.savefig(name)
    plt.show()
    return

def armarGraficoRK2Error(t,aproximacion,h,titulo,nombre):
    plt.plot(t,aproximacion)
    plt.xlabel('t')
    plt.ylabel('error')
    plt.title(titulo + str(h))
    name = nombre+str(h)+'.png'
    plt.savefig(name)
    plt.show()
    return

def iterarRK2(aproximacion,intervalo,h):
    paso=0
    for j in range (int((intervalo/h)+1)):
        aproximacion[j] = (utils.analitica(paso) - aproximacion[j])
        paso += h
    return

def imprimirRungeKutta2(f,intervalo,h,t):
    aproximacion = runge_kuta_2(f,intervalo, h)
    armarGraficoRK2OK(aproximacion,t,h,utils.y_t_RK2,utils.y_nombre_RK2)
    iterarRK2(aproximacion,intervalo,h)
    armarGraficoRK2Error(aproximacion,t,h,utils.e_t_RK2,utils.e_nombre_RK2)