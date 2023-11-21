import numpy as np

# Datos del problema
padronJI = 106957
padronTB = 102665
padronCZ = 106260
padron = np.sum([padronJI,padronTB,padronCZ])
denominador = 200
m = padron/denominador  # masa (kg)
k = 25000  # constante elástica del muelle (N/m)
lambda_ = 0  # constante de amortiguación (Ns/m)
c = 0.1  # cota del terreno (m)

# Condiciones iniciales
y0 = 0.1  # posición inicial
v0 = 0.0  # velocidad inicial

# Parámetros de simulación
delta_tiempo = 0.005  # segundos
c_prim = 0 
h = 0.005
intervalo = 5
t=np.arange(0,intervalo+h,h)

y_t_EulerExplicito = 'y(t) euler explicito con paso = '
y_nombre_EulerExplicito = 'EulerExplicito'
e_t_EulerExplicito = 'e(t) euler explicito con paso = '
e_nombre_EulerExplicito = 'error-EulerExplicito'


y_t_EulerImplicito = 'y(t) euler implicito con paso = '
y_nombre_EulerImplicito = 'EulerImplicito'
e_t_EulerImplicito = 'e(t) euler implicito con paso = '
e_nombre_EulerImplicito = 'error-EulerImplicito'

y_t_RK2 = 'y(t) RK2 con paso = '
y_nombre_RK2 = 'RK2'
e_t_RK2 = 'e(t) RK2 con paso = '
e_nombre_RK2 = 'error-RK2'


'Solucion analitica sin amortiguacion'
def analitica(t):
    return c - c * np.cos(((k/m)**0.5)*t)

maxiCompresion = -0.05
compresionFinal = -10000000
minimaPonderacion = 10000000000
lambdaAmortiguado = 150
kAmortiguado = 25000
lambdaV = np.arange(150,12000,50)
kV = np.arange(2500,100000,1000)

