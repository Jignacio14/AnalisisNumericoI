import numpy as np

# Método de Euler explícito para ecuación diferencial de segundo orden
def euler_explicito(y0, v0, m, k, lambda_, c, tiempo_total, paso_tiempo):
    cantidad_pasos = int(tiempo_total / paso_tiempo) + 1

    # Inicializar arrays para almacenar las soluciones
    y = np.zeros(cantidad_pasos)
    v = np.zeros(cantidad_pasos)

    # Condiciones iniciales
    y[0] = y0
    v[0] = v0

    for n in range(1, cantidad_pasos):
        # Método de Euler explícito
        # u_{n+1} = u_n + h * f(u_n, t_n) 
        # f_1(u_n,v_n,t_n) = v_n
        # u_{n+1} = u_n + h * v_n 
        # v_{n+1} = v_n + h f_2(u_n,v_n,t_n)
        y[n] = y[n-1] + paso_tiempo * v[n-1]
        v[n] = v[n-1] + paso_tiempo * ((-k / m) * (y[n-1] - c) - (lambda_ / m) * v[n-1])

    return y, v

# Método de Euler implícito para ecuación diferencial de segundo orden
def euler_implicito(y0, v0, m, k, lambda_, c, tiempo_total, paso_tiempo):
    cantidad_pasos = int(tiempo_total / paso_tiempo) + 1

    # Inicializar arrays para almacenar las soluciones
    y = np.zeros(cantidad_pasos)
    v = np.zeros(cantidad_pasos)

    y[0] = y0
    v[0] = v0

    for n in range(1, cantidad_pasos):
        # Método de Euler implícito
        y[n] = (y[n-1] + paso_tiempo * v[n-1]) / (1 + (paso_tiempo * lambda_ / m))
        v[n] = v[n-1] + paso_tiempo * ((-k / m) * (y[n] - c) - (lambda_ / m) * v[n-1])

    return y, v