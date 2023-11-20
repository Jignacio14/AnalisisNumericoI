import matplotlib.pyplot as plt
import rungekutta2
import utils
import numpy as np

def c(t):
    if 1 <= t < 1.1:
        return t - 1.0
    elif 1.1 <= t < 1.3:
        return 0.1
    elif 1.3 <= t < 1.4:
        return 1.4 - t
    else:
        return 0.0

def c_prim(t):
    if t < 1:
        return 0
    elif 1 <= t < 1.1:
        return 1
    elif 1.1 <= t < 1.3:
        return 0
    elif 1.3 <= t < 1.4:
        return -1
    else:
        return 0

def aceleracion(c, t, c_prim, u, v):
    return (utils.k / utils.m) * (c(t) - u) + (utils.lambda_ / utils.m) * (c_prim(t) - v)

def maximaCompresion(intervalo, h, extras, c):
    paso = 0
    aproximada, v = rungekutta2.rungeKuttaO2ConExtras(utils.f2, intervalo, h, extras)
    minimo = 0
    for j in range(int((intervalo / h) + 1)):
        aproximada[j] = (aproximada[j] - float(c(paso)))
        if aproximada[j] < minimo:
            minimo = aproximada[j]
        paso += h
    return minimo

def f1(y, c, t, y_prim, c_prim):
    lambda_ = 750
    res = (utils.k / utils.m) * (c(t) - y) + (lambda_ / utils.m) * (c_prim(t) - y_prim)
    return res

def amortiguar():
    v_n = np.vectorize(c)
    x = np.arange(0, utils.intervalo + utils.h, utils.h)
    y = v_n(x)
    plt.plot(x, y, 'r')
    
    t = np.arange(0, utils.intervalo + utils.h, utils.h)
    extras = [c, c_prim]

    rungekutta2.imprimirRungeKutta2Extras(v_n, utils.intervalo, utils.h, extras, t, c, 15000, 500)
    
    # AGREGAR MÁS CASOS

    for i in range(int(utils.lambdaV.shape[0])-1):
        for j in range(int(utils.kV.shape[0])-1):
            act = 0
            act = maximaCompresion(utils.intervalo, utils.h, extras, c)
            if act >= utils.maxiCompresion:
                ponderacionActual = utils.lambdaV[i] / 750 + utils.kV[j] / utils.k
                if ponderacionActual < minimaPonderacion:
                    minimaPonderacion = ponderacionActual
                    utils.lambdaAmortiguado = utils.lambdaV[i]
                    utils.kAmortiguado = utils.kV[j]
                    utils.compresionFinal = act
                    
    print('K elegido = ' + str(utils.kAmortiguado) + ' y lambda electo = ' + str(utils.lambdaAmortiguado) +
          'con compresion = ' + str(utils.compresionFinal))
    
    utils.kAmortiguado = 80500
    utils.lambdaAmortiguado = 4550
    t = np.arange(0, utils.intervalo + utils.h, utils.h)
    u, v = rungekutta2.imprimirRungeKutta2Extras(v_n, utils.intervalo, utils.h, extras, t, c, utils.kAmortiguado, utils.lambdaAmortiguado)

    afun = np.vectorize(aceleracion)
    x = np.arange(0, utils.intervalo + utils.h, utils.h)
    y = afun(c, x, c_prim, u, v)
    plt.xlabel('t')
    plt.ylabel('y\'\'(t)')
    plt.title('Aceleracion')
    plt.plot(x, y, 'r')
    plt.show()
    return
