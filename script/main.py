import euler
import rungekutta2

def main():
    # Datos del problema
    padronJI = 0
    padronTB = 0
    padronCZ = 106260
    denominador = 200
    m1 = padronJI/denominador  # masa (kg)
    m2 = padronTB/denominador
    m3 = padronCZ/denominador
    k = 25000  # constante elástica del muelle (N/m)
    lambda_ = 0  # constante de amortiguación (Ns/m)
    c = 0.1  # cota del terreno (m)

    # Condiciones iniciales
    y0 = 0.1  # posición inicial
    v0 = 0.0  # velocidad inicial

    # Parámetros de simulación
    delta_tiempo = 0.005  # segundos


    return

if __name__ == "__main__":
    main()