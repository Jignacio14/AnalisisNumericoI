import numpy
import euler
import rungekutta2

def main():
    # Datos del problema
    padronJI = 106957
    padronTB = 102665
    padronCZ = 106260
    padron = numpy.sum([padronJI,padronTB,padronCZ])
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

    euler.euler_explicito(y0,v0,m,k,lambda_,c,0,delta_tiempo)
    euler.euler_implicito(y0,v0,m,k,lambda_,c,0,delta_tiempo)

    rungekutta2.runge_kutta_segundo_orden(y0,v0,m,k,lambda_,c,0,delta_tiempo)


    return

if __name__ == "__main__":
    main()