import utils
import matplotlib.pyplot as plt
import numpy as np

def euler_explicito(f,intervalo, paso):
    cantidad = intervalo/paso +1 
    u=np.zeros(int(cantidad),dtype=float)
    v=np.zeros(int(cantidad),dtype=float)
    u[0]=0
    v[0]=0
    i=0
    n=0
    while n<cantidad-1:
        u_actual = (u[(n)] + paso * v[n])
        u[(n+1)]= u_actual
        v_actual = v[n] + paso * f(u[n])
        v[n+1] = v_actual
        n+=1
        i+=paso
    return u

def euler_implicito(A_invertida,termino_indep,intervalo, paso):
    cantidad = intervalo/paso +1
    res = np.zeros(int(cantidad*2))
    res.shape=(2,int(cantidad))
    i=0
    n=0
    while n<cantidad-1:
        aux = np.array(A_invertida @ res[:,n] + termino_indep)
        res[:,(n+1)]= aux
        n+=1
        i+=paso
    return res[1]


def armarGraficoEulerOK(aproximacion,t,h,titulo,nombre):
    plt.plot(t,aproximacion)
    plt.xlabel('t')
    plt.ylabel('y')
    plt.title(titulo + str(h))
    plt.plot(t,utils.analitica(t) , 'r--')
    name = nombre+str(h)+'.png'
    plt.savefig(name)
    plt.show()
    return

def armarGraficoEulerError(t,aproximacion,h,titulo,nombre):
    plt.plot(t,aproximacion)
    plt.xlabel('t')
    plt.ylabel('error')
    plt.title(titulo + str(h))
    name = nombre+str(h)+'.png'
    plt.savefig(name)
    plt.show()
    return

def iterarEuler(aproximacion,intervalo,h):
    paso = 0
    for i in range(int((intervalo/h)+1)):
        aproximacion[i] = (utils.analitica(paso) - aproximacion[i])
        paso+=h

def imprimirEulerExplicito(f,intervalo,h,t):
    aproximacion = euler_explicito(f,intervalo,h)
    armarGraficoEulerOK(aproximacion,t,h,utils.y_t_EulerExplicito,utils.y_nombre_EulerExplicito)
    iterarEuler(aproximacion,intervalo,h)
    armarGraficoEulerError(t,aproximacion,h,utils.e_t_EulerExplicito,utils.e_nombre_EulerExplicito)
    return

def imprimirEulerImplicito(A_invertida, termino_indep,intervalo,h,t):
    aproximacion = euler_implicito(A_invertida,termino_indep,intervalo,h)
    armarGraficoEulerOK(aproximacion,t,h,utils.y_t_EulerImplicito,utils.y_nombre_EulerImplicito)
    iterarEuler(aproximacion,intervalo,h)
    armarGraficoEulerError(t,aproximacion,h,utils.e_t_EulerImplicito,utils.e_nombre_EulerImplicito)
    return