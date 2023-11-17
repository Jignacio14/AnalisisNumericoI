import numpy as np

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