import numpy as np

# Método de Euler explícito para ecuación diferencial de segundo orden
def euler_explicito(y0, y1_0, f, intervalo, paso):
    cantidad_pasos = int(intervalo / paso) + 1

    # Inicializar arrays para almacenar las soluciones
    y0_solucion = np.zeros(cantidad_pasos)
    y1_solucion = np.zeros(cantidad_pasos)

    # Condiciones iniciales
    y0_solucion[0] = y0
    y1_solucion[0] = y1_0

    for n in range(1, cantidad_pasos):
        # Método de Euler explícito
        # u_{n+1} = u_n + h * f(u_n, t_n) 
        # f_1(u_n,v_n,t_n) = v_n
        # u_{n+1} = u_n + h * v_n 
        # v_{n+1} = v_n + h f_2(u_n,v_n,t_n)
        y0_solucion[n] = y0_solucion[n-1] + paso * y1_solucion[n-1]
        y1_solucion[n] = y1_solucion[n-1] + paso * f(y0_solucion[n-1], y1_solucion[n-1])

    return y0_solucion, y1_solucion

# Método de Euler implícito para ecuación diferencial de segundo orden
def euler_implicito(y0, y1_0, f, intervalo, paso):
    cantidad_pasos = int(intervalo / paso) + 1

    # Inicializar arrays para almacenar las soluciones
    y0_solucion = np.zeros(cantidad_pasos)
    y1_solucion = np.zeros(cantidad_pasos)

    # Condiciones iniciales
    y0_solucion[0] = y0
    y1_solucion[0] = y1_0

    for n in range(1, cantidad_pasos):
        # Método de Euler implícito
        y0_siguiente = y0_solucion[n-1] + paso * y1_solucion[n-1]
        y1_siguiente = y1_solucion[n-1] + paso * f(y0_siguiente, y1_solucion[n-1])

        # Actualizar soluciones utilizando el resultado implícito
        y0_solucion[n] = y0_siguiente
        y1_solucion[n] = y1_siguiente

    return y0_solucion, y1_solucion

# Ecuación diferencial de segundo orden
def f(y0, y1):
    k = 1.0
    m = 1.0
    c = 2.0
    lambda_ = 0.1
    c0_prime = 1.0

    return (k / m) * (c - y0) + (lambda_ / m) * (c0_prime - y1)

def main():
    # Condiciones iniciales
    y0_0 = 1.0
    y1_0 = 0.0

    # Intervalo de tiempo y tamaño del paso
    intervalo_tiempo = 5.0
    paso = 0.1

    # Aplicar el método de Euler explícito
    solucion_y0_exp, solucion_y1_exp = euler_explicito(y0_0, y1_0, f, intervalo_tiempo, paso)

    # Aplicar el método de Euler implícito
    solucion_y0_imp, solucion_y1_imp = euler_implicito(y0_0, y1_0, f, intervalo_tiempo, paso)

    # Imprimir o visualizar los resultados según sea necesario
    print("Resultados de Euler Explícito:")
    print("y0:", solucion_y0_exp)
    print("y1:", solucion_y1_exp)

    print("\nResultados de Euler Implícito:")
    print("y0:", solucion_y0_imp)
    print("y1:", solucion_y1_imp)

if __name__ == "__main__":
    main()
