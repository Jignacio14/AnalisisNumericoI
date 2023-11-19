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

h = 0.005
intervalo = 5
t=np.arange(0,intervalo+h,h)

y_t_EulerExplicito = 'y(t) euler explicito con paso = '
y_nombre_EulerExplicito = 'EulerExplicito'
e_t_EulerExplicito = 'e(t) euler explicito con paso = '
e_nombre_EulerExplicito = 'error-EulerExplicito'


y_t_EulerImplicito = 'y(t) euler implicito con paso = '
y_nombre_EulerImplicito = 'EulerImplicito'
e_t_EulerImplicito = 'e(t) euler Implicito con paso = '
e_nombre_EulerImplicito = 'error-EulerImplicito'