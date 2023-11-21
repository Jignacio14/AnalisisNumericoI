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
    lambda_ = 750
    return (utils.k / utils.m) * (c(t) - u) + (lambda_ / utils.m) * (c_prim(t) - v)

def maximaCompresion(intervalo, h, extras, c,lam,k):
    paso = 0
    aproximada, v = rungekutta2.rungeKuttaO2ConExtras(intervalo, h, extras,lam,k)
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

    # varianza de lambda
    rungekutta2.imprimirRungeKutta2Extras(v_n, utils.intervalo, utils.h, extras, t, c, 15000, 500)
    rungekutta2.imprimirRungeKutta2Extras(v_n, utils.intervalo, utils.h, extras, t, c, 15000, 1000)
    rungekutta2.imprimirRungeKutta2Extras(v_n, utils.intervalo, utils.h, extras, t, c, 15000, 1500)
    rungekutta2.imprimirRungeKutta2Extras(v_n, utils.intervalo, utils.h, extras, t, c, 15000, 2000)
    rungekutta2.imprimirRungeKutta2Extras(v_n, utils.intervalo, utils.h, extras, t, c, 15000, 2500)
    rungekutta2.imprimirRungeKutta2Extras(v_n, utils.intervalo, utils.h, extras, t, c, 15000, 3000)
    rungekutta2.imprimirRungeKutta2Extras(v_n, utils.intervalo, utils.h, extras, t, c, 15000, 3500)
    
    # varianza de k
    rungekutta2.imprimirRungeKutta2Extras(v_n, utils.intervalo, utils.h, extras, t, c, 15000, 500)
    rungekutta2.imprimirRungeKutta2Extras(v_n, utils.intervalo, utils.h, extras, t, c, 20000, 500)
    rungekutta2.imprimirRungeKutta2Extras(v_n, utils.intervalo, utils.h, extras, t, c, 25000, 500)
    rungekutta2.imprimirRungeKutta2Extras(v_n, utils.intervalo, utils.h, extras, t, c, 30000, 500)
    rungekutta2.imprimirRungeKutta2Extras(v_n, utils.intervalo, utils.h, extras, t, c, 35000, 500)
    rungekutta2.imprimirRungeKutta2Extras(v_n, utils.intervalo, utils.h, extras, t, c, 40000, 500)

    # prueba de otra varianza, fidelidad de datos
    rungekutta2.imprimirRungeKutta2Extras(v_n, utils.intervalo, utils.h, extras, t, c, 30000, 500)
    rungekutta2.imprimirRungeKutta2Extras(v_n, utils.intervalo, utils.h, extras, t, c, 30000, 1000)
    rungekutta2.imprimirRungeKutta2Extras(v_n, utils.intervalo, utils.h, extras, t, c, 30000, 1500)
    rungekutta2.imprimirRungeKutta2Extras(v_n, utils.intervalo, utils.h, extras, t, c, 30000, 2000)
    rungekutta2.imprimirRungeKutta2Extras(v_n, utils.intervalo, utils.h, extras, t, c, 30000, 2500)
    rungekutta2.imprimirRungeKutta2Extras(v_n, utils.intervalo, utils.h, extras, t, c, 30000, 3000)
    rungekutta2.imprimirRungeKutta2Extras(v_n, utils.intervalo, utils.h, extras, t, c, 30000, 3500)

    lambdaV = utils.lambdaV
    kV = utils.kV
    lambdaAmortiguado = utils.lambdaAmortiguado
    kAmortiguado = utils.kAmortiguado
    compresionFinal = utils.compresionFinal
    minimaPonderacion = utils.minimaPonderacion
    lam = 750
    for i in range(int(lambdaV.shape[0])-1):
        for j in range(int(kV.shape[0])-1):
            act = 0
            act = maximaCompresion(utils.intervalo, utils.h, extras, c,lambdaV[i],kV[j])
            if act >= utils.maxiCompresion:
                ponderacionActual = lambdaV[i] / lam + kV[j] / utils.k
                if ponderacionActual < minimaPonderacion:
                    minimaPonderacion = ponderacionActual
                    lambdaAmortiguado = lambdaV[i]
                    kAmortiguado = kV[j]
                    compresionFinal = act
                    
    print('K elegido = ' + str(kAmortiguado) + ' y lambda electo = ' + str(lambdaAmortiguado) +
          'con compresion = ' + str(compresionFinal))
    
    kAmortiguado = 10550
    lambdaAmortiguado = 5550
    t = np.arange(0, utils.intervalo + utils.h, utils.h)
    u, v = rungekutta2.imprimirRungeKutta2Extras(v_n, utils.intervalo, utils.h, extras, t, c, kAmortiguado, lambdaAmortiguado)

    afun = np.vectorize(aceleracion)
    x = np.arange(0, utils.intervalo + utils.h, utils.h)
    y = afun(c, x, c_prim, u, v)
    plt.xlabel('t')
    plt.ylabel('y\'\'(t)')
    plt.title('Aceleracion')
    plt.plot(x, y, 'r')
    plt.show()
    return
