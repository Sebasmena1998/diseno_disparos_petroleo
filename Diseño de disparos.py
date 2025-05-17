print("*** Diseño de disparos ***\n")
# Import
import math
import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate
from scipy.interpolate import interp1d

## Funciones
def pausar():
    input("Presiona Enter para continuar con el siguiente paso...")

def redondeo(valor):
    decimal = valor - int(valor)
    if decimal < 0.5:
        return int(valor)
    else:
        return int(valor) + 1

def estimar_tiros_por_pie(penetracion):
    datos = {
        1:12,
        3: 8,
        4: 6,
        6: 4,
        10: 2,
        20: 1
    }

    puntos = sorted(datos.items())

    if penetracion <= puntos[0][0]:
        return datos[puntos[0][0]]
    if penetracion >= puntos[-1][0]:
        return datos[puntos[-1][0]]

    for i in range(len(puntos) - 1):
        x0, y0 = puntos[i]
        x1, y1 = puntos[i + 1]
        if x0 <= penetracion <= x1:
            # Interpolación lineal
            tpf = y0 + (y1 - y0) * (penetracion - x0) / (x1 - x0)
            return redondeo(tpf)

def temperatura_limite(m, b, tiempo):
    return m * math.log10(tiempo) + b

def elegir_explosivo(tiempo, temperatura):
    explosivos = {
        "PYX": (-22, 555),
        "HNS": (-21, 480),
        "HMX": (-20, 390),
        "RDX": (-18, 310)
    }

    adecuados = []
    for nombre_explosivo, (m, b) in explosivos.items():
        temp_max = temperatura_limite(m, b, tiempo)
        if temperatura <= temp_max:
            adecuados.append(nombre_explosivo)

    if adecuados:
        return adecuados[-1]  # Retorna el menos resistente que aún es válido
    else:
        return "Ninguno (la temperatura es demasiado alta para estos explosivos)"

def interpolar_constantes(fase_usuario):
    tabla_constantes = {
        0:    [0.250, -2.091, 0.0453, 5.1313, 1.8672, 0.1600, 2.675],
        180:  [0.500, -2.025, 0.0943, 3.0373, 1.8111, 0.0620, 4.532],
        120:  [0.648, -2.018, 0.0634, 1.6136, 1.7770, 0.0650, 5.320],
        90:   [0.726, -1.985, 0.1038, 1.5674, 1.6935, 0.0910, 6.155],
        60:   [0.813, -1.898, 0.1023, 1.3654, 1.6490, 0.0030, 7.509],
        45:   [0.860, -1.788, 0.2398, 1.1915, 1.6392, 0.000046, 8.791]
    }

    fases_ordenadas = sorted(tabla_constantes.keys())

    if fase_usuario in tabla_constantes:
        return tabla_constantes[fase_usuario]

    for i in range(len(fases_ordenadas) - 1):
        f1, f2 = fases_ordenadas[i], fases_ordenadas[i + 1]
        if f1 < fase_usuario < f2:
            c1, c2 = tabla_constantes[f1], tabla_constantes[f2]
            interpoladas = [
                c1[j] + (c2[j] - c1[j]) * (fase_usuario - f1) / (f2 - f1)
                for j in range(len(c1))
            ]
            return interpoladas
    return None

def obtener_RP(beta_0):
    log_beta_calc = np.log10(beta_0)
    log_RP = interp_func(log_beta_calc)
    return 10 ** log_RP

# Función para calcular beta_0
def calcular_beta0(Lpf, tpf_est, dperf, anisotropia):
    return Lpf * (tpf_est ** 1.5) * (dperf ** 0.5) * (anisotropia ** (-5 / 8))


#--------------- Paso 1: Determinar el tipo de Formación (Consolidada o No consolidada) ---------------------------
print("*** Favor de ingresar la información necesaria ***\n")
nombre_pozo = input("Ingrese el nombre del pozo: ")
# Preguntar si se tienen datos de registros sónicos o de densidad
while True:
    respuesta_datos_registro = input("¿Tiene datos de registros sónicos o de densidad? (si/no): ").lower().strip()

    if respuesta_datos_registro == "si":
        while True:
            sonico_o_densidad = input(
                "¿Tiene los datos del registro sónico en el intervalo de interés? (si/no): ").lower().strip()
            if sonico_o_densidad == "si":
                try:
                    valor_sonico = float(input("Introduce el valor del registro sónico (μs/pie): "))
                    formacion = "Consolidada" if valor_sonico <= 100 else "No consolidada"
                    break
                except ValueError:
                    print("\nError: Por favor ingrese un valor numérico válido.")
            elif sonico_o_densidad == "no":
                try:
                    valor_densidad = float(input("\nIntroduce el valor del registro de densidad (grs/cc): "))
                    formacion = "Consolidada" if valor_densidad >= 2.4 else "No consolidada"
                    break
                except ValueError:
                    print("\nError: Por favor ingrese un valor numérico válido.")
            else:
                print("\nOpción no válida. Por favor responda 'si' o 'no'.")
        break

    elif respuesta_datos_registro == "no":
        while True:
            seleccion_manual = input("\n¿La formación es consolidada o no consolidada? (c/nc): ").lower().strip()
            if seleccion_manual == "c":
                formacion = "Consolidada"
                break
            elif seleccion_manual == "nc":
                formacion = "No consolidada"
                break
            else:
                print("\nOpción no válida. Por favor ingrese 'c' o 'nc'.")
        break

    else:
        print("\nOpción no válida. Por favor responda 'si' o 'no'.")
#--------------- Paso 2: Determinar la permeabilidad de la formación ---------------------------
while True:
    # Preguntar si se tienen datos directos de permeabilidad
    respuesta_permeabilidad = input("¿Tiene datos medidos de permeabilidad? (si/no): ").lower().strip()

    if respuesta_permeabilidad == "si":
        # Obtener permeabilidad directamente
        while True:
            try:
                permeabilidad = float(input("Ingrese el valor de permeabilidad en milidarcies (mD): "))
                if permeabilidad >= 0:  # Validar que no sea negativo
                    break
                else:
                    print("Error: La permeabilidad no puede ser negativa. Intente nuevamente.")
            except ValueError:
                print("Error: Debe ingresar un valor numérico válido.")
        break

    elif respuesta_permeabilidad == "no":
        # Calcular permeabilidad usando porosidad y saturación de agua
        print("\nCalcularemos la permeabilidad usando la porosidad (ϕ) y saturación de agua (Sw)")

        # Obtener porosidad con validación
        while True:
            try:
                porosidad = float(input("Ingrese la porosidad en porcentaje (%): "))
                if 0 <= porosidad <= 100:
                    break
                else:
                    print("Error: La porosidad debe estar entre 0% y 100%.")
            except ValueError:
                print("Error: Ingrese un valor numérico válido.")

        # Obtener saturación de agua con validación
        while True:
            try:
                saturacion_agua = float(input("Ingrese la saturación de agua en porcentaje (%): "))
                if 0 <= saturacion_agua <= 100:
                    break
                else:
                    print("Error: La saturación debe estar entre 0% y 100%.")
            except ValueError:
                print("Error: Ingrese un valor numérico válido.")

        # Calcular permeabilidad usando la ecuación de Timur
        permeabilidad = 0.136 * ((porosidad ** 4.4) / (saturacion_agua ** 2))
        break

    else:
        print("Opción no válida. Por favor responda 'si' o 'no'.\n")
#--------------- Paso 3: Determinar el mínimo desbalance para superar el factor de daño ---------------------------
while True:
    tipo_fluido = input("¿Tipo de fluido en el intervalo, Gas o Aceite? (g/a): ").lower().strip()

    if tipo_fluido == "g":
        Umin = 2500 / (permeabilidad ** 0.17)
        break
    elif tipo_fluido == "a":
        Umin = 3500 / (permeabilidad ** 0.37)
        break
    else:
        print("Error: Opción no válida. Ingrese 'g' o 'a'.\n")

##--------------- Paso 4: Determinar el maximo desbalance seguro ---------------------------
if respuesta_datos_registro == "no" or formacion == "Consolidada":
    print('Debido a que la formación es "Consolidada" o no cuenta con los datos de los registros.')
    print("Podemos considerar UM como el limite de presion del casing, tuberia o accesorio ")
    try:
        UMax = float(input("Ingrese el limite de presion en PSI: "))
        UMax = UMax*0.80
    except ValueError:
        print("\nError: Por favor ingrese un valor numérico válido.")
else:
    if tipo_fluido == "a":
        if sonico_o_densidad == "si":
            UMax = 3600 - 20*valor_sonico
        else:
            UMax = 2340*valor_densidad-4000
    else:
        if sonico_o_densidad == "si":
            UMax = 4750 - 25*valor_sonico
        else:
            UMax = 2900*valor_densidad-4700

##--------------- Paso 5: Calcular el promedio entre Umin y UMax ---------------------------
Uprom = (Umin+UMax)/2

## Considerando como estandar que la invasion sea somera
U= (Umin+Uprom)/2
## Calculamos la PH
while True:
    try:
        PRESION_YACIMIENTO = float(input("Ingrese la presión de yacimiento en PSI: "))
        DENSIDAD_FLUIDO_TERMINACION = float(input("Ingrese la densidad del fluido de terminación en gr/cc: "))
        FONDO_POZO = float(input("Ingrese la profundidad de fondo de pozo en metros (m): "))
        presion_hidrostatica = PRESION_YACIMIENTO - U
        altura_columna = presion_hidrostatica/(1.4228*DENSIDAD_FLUIDO_TERMINACION)
        NIVEL_FLUIDO = FONDO_POZO-altura_columna
        break
    except ValueError:
        print("\nError: Por favor ingrese un valor numérico válido.")
##--------------- Paso 6: Determinar la longitud estimada de penetración ---------------------------
while True:
    try:
        if respuesta_permeabilidad == "si":
            porosidad = float(input("Ingrese la porosidad en porcentaje (%): "))
            break
        respuesta_compresibilidad_formacion = input("¿Tiene el dato de la compresibilidad de la formación? (si/no): ").lower().strip()
        if respuesta_compresibilidad_formacion == "si":
            compresibilidad_formacion = float(input("Ingrese la compresibilidad de la formación en psi: "))
            break
        else:
            compresibilidad_formacion = 28510 - (1023.3*porosidad)
            break
    except ValueError:
        print("\nError: Por favor ingrese un valor numérico válido.")

##Manejando como estandar una   Arenisca berea con una Compresibilidad de 7000 pies y una penetración de 9.21 oulgadas
PENETRACION = 9.21
RESISTENCIA_COMPRESIVA = 7000

lnLpf = math.log(PENETRACION) + 0.086 * (RESISTENCIA_COMPRESIVA-compresibilidad_formacion) * (1*10**-3)
Lpf = math.e**lnLpf

##--------------- Paso 7: Determinar la densidad de disparo ---------------------------

tpf_estimado = estimar_tiros_por_pie(Lpf)

if formacion == "Consolidada":
    while True:
        tipo_terminacion = input("Introduzca el tipo de terminación: Natural (n), Estimulada (e): ").lower().strip()
        if tipo_terminacion == "n":
            nombre_terminacion = "Natural"
            tpf_estimado *= 3
            tipo_carga = "DP"
            fase = 360 / tpf_estimado
            break
        elif tipo_terminacion == "e":
            nombre_terminacion = "Estimulada"
            tpf_estimado *= 3
            tipo_carga = "GH"
            fase = 180 / tpf_estimado
            break
        else:
            print("\nOpción no válida. Por favor ingrese 'n' o 'e'.")
else:
    tpf_estimado *= 5
    tipo_carga = "BH"
    fase = 360 / tpf_estimado
    nombre_terminacion = "Control de arenas"




##--------------- Paso 8: Determinar el tipo de explosivo ---------------------------
while True:
    profundidad_base= float(input("Inserte la profundidad de la base del intervalo deseado en metros: "))
    profundidad_cima= float(input("Inserte la profundidad de la cima del intervalo deseado en metros: "))
    respuesta_tiempo_exposicion = input("¿Conoce el tiempo de exposición al que va a estar sometido la pistola? (si/no): ").strip().lower()
    if respuesta_tiempo_exposicion == "si":
        tiempo_exposicion = float(input("Ingrese el tiempo de exposición en hrs: "))
        tiempo_exposicion_seguro= tiempo_exposicion*1.5
        break
    elif respuesta_tiempo_exposicion =="no":
        #Considerando un promedio de 10 minutos por tuberia y sabiendo que una tuberia mide 9 metros
        LONGITUD_TRAMO =9
        tiempo_exposicion= ((profundidad_base / LONGITUD_TRAMO) * 10) / 60
        #aplicando un margen de seguridad de 50%
        tiempo_exposicion_seguro = tiempo_exposicion*1.5
        break
    else:
        print("\nOpción no válida. Por favor responda 'si' o 'no'.")

temperatura_intervalo_centigrados = float(input("Ingrese la temperatura en el intervado en Centrigrados (°C) "))
temperatura_intervalo_f = (temperatura_intervalo_centigrados*(9/5))+32
temperatura_intervalo_f_seguro = temperatura_intervalo_f*1.1
temperatura_intervalo_centigrados_seguro = temperatura_intervalo_centigrados*1.1

explosivo = elegir_explosivo(tiempo_exposicion_seguro, temperatura_intervalo_f_seguro)
##--------------- Paso 9: Calcular el daño por efecto del disparo ---------------------------

# --- Entrada de datos requeridos ---
try:
    diametro_pozo_pulgadas = float(input("Inserte el diámetro del pozo en pulgadas: "))
except ValueError:
    print("\nError: Por favor ingrese valores numéricos válidos.")
    exit()

rw_pulgadas = diametro_pozo_pulgadas/2
# --- Obtener constantes interpoladas según fase ---
constantes = interpolar_constantes(fase)
if constantes is None:
    print("Fase fuera del rango. No se puede calcular.")
    exit()

# Asignación de constantes
A0, A1, A2, B1, B2, C1, C2 = constantes

# --- Calculo r'w(0) ---
rw_pies = (rw_pulgadas)/12  # pasar

if fase not in [0, 360]:
    r_w0 = A0 * (rw_pies + (Lpf / 12))
else:
    r_w0 = (Lpf / 12) / 4

# --- Calcular SH ---
sh = math.log((rw_pies/r_w0))
# --- Valores para rperf según tipo de carga ---
if tipo_carga == "DP":
    rperf = 0.365
elif tipo_carga == "GH":
    rperf = 0.5725
else:
    rperf = 0.78

dperf = rperf*2

# --- Parámetros geométricos y de anisotropía ---
ANISOTROPIA = 10
hperf = 1 / tpf_estimado
rD = ((rperf / 2) * (1 / 12)) / (2 * hperf) * (1 + math.sqrt(1 / ANISOTROPIA))
hD = (hperf / (Lpf / 12)) * math.sqrt(ANISOTROPIA)
# --- Calculo de a y b ---
a = (A1 * math.log(rD,10))+ A2
b = B1 * rD + B2
# --- Calcular SV ---
SV = (10 ** a) * (hD ** (b - 1)) * (rD ** b)
# --- Calcular Swb ---
rwD = rw_pies / ((Lpf / 12) + rw_pies)
Swb = C1 * math.exp(C2 * rwD)
Sp= sh+SV+Swb

##--------------- Paso 10: Determinar el tipo de disparo ---------------------------
if profundidad_base >1500:
    tipo_disparo = "TCP"
else:
    angulo_maximo = float(input("Ingrese la maximo angulo de inclinación del pozo: "))
    if angulo_maximo >65:
        tipo_disparo = "Wireline"
    else:
        tipo_disparo = "TCP"


#-------------- Paso 11: Calculos de RP -----------------------------------------
beta_0_points = np.array([0.1, 1, 10, 100, 1000])  # Valores de beta_0 (log)
RP_N4_points = np.array([0.15, 0.35, 0.6, 0.8, 0.95])  # RP para N=4 SPF

# Crear función de interpolación (log-log)
log_beta = np.log10(beta_0_points)
log_RP_N4 = np.log10(RP_N4_points)
interp_func = interp1d(log_beta, log_RP_N4, kind='linear', fill_value='extrapolate')


# Llamamos a la funcion Calcular beta_0 y Obtener Relación de Productividad (RP)
beta_0 = calcular_beta0(Lpf, tpf_estimado, dperf, ANISOTROPIA)
RP_grafica = obtener_RP(beta_0)

#RP_ecuación
re_pulgadas=float(input("Inserte el radio de drene del pozo (re) en pulgadas: "))
re_pies= re_pulgadas*(1/12)
RP_ecuacion = math.log(re_pies/rw_pies)/(math.log(re_pies/rw_pies)+ Sp)





## ------------ Resumen de resultados -------------------------------


print(f'''\n*** Los Resultados son los siguientes:

El pozo {nombre_pozo}, busca dispara en el intervalo {profundidad_cima}m - {profundidad_base}m,
un intervalo de un espesor de {profundidad_base-profundidad_cima}m, al ser una formación {formacion}
con una permeabilidad de {permeabilidad:.4f}mD, se busca una presión bajo-balance por sus beneficios al momento
de disparar como la limpieza de los canales, es por eso que con una Presión Hidrostatica (PH) de {presion_hidrostatica:.4f} psi
el diferencial de balance calculado es de {U:.4f} psi, entonces la columna de fluido de control para obtener este
diferencial es de {altura_columna:.2f} metros, dandonos un nivel de fluido de {NIVEL_FLUIDO:.2f}.
Gracias a los estudios de Thopmson podemos determinar la longitud de penetración estimada que es de {Lpf:.2f} pulgadas,
y al tener una terminación {nombre_terminacion}, Podemos determinar que 
la desindad de tiro adecuada es de {tpf_estimado} spf con una fase de {fase}° con un tipo de carga {tipo_carga}
lo cual ayudara a las carecteristicas de la formación a tener un mejor Indice de productividad.
Con todos los valores obtenidos determinacmos que este disparo puede generar un daño por disparos (Sp) de {Sp:.4f}.
Por los datos de la profundidad o del angulo maximo de inclinación se recomienda que sean disparos {tipo_disparo},
buscando elegir la pistola de mayor tamaño disponible.
''')

pausar()

# Gráfica de referencia
plt.figure(figsize=(8, 5))
plt.loglog(beta_0_points, RP_N4_points, 'ro-',)
plt.loglog(beta_0, RP_grafica, 'bs', label=f'RP={RP_grafica:.2f}')
plt.xlabel('Factor βo (log)')
plt.ylabel('Relación de Productividad (RP)')
plt.title('Relación entre βo y RP')
plt.grid(True, which="both", ls="--")
plt.legend()
plt.show()

print(f'''Con todos los datos obtenidos podemos considerar su afectación en la RP:
Con un factor adimencional calculado βo de {beta_0} obtenemos 
Un relación de productividad de {RP_grafica:.2f} de la grafica en relación a su βo y un {RP_ecuacion:.2f} por ecuación.
''')

# -------------------- Tabla de diseño de disparos ------------------------------

headers = [
    "Cima", "Base", "Espesor", "Temperatura °F", "Tipo de Formación", "Tipo de Terminación",
    "Densidad de tiro (spf)", "Fase de tiro", "Tipo de Carga", "Daño por disparos", "RP Gráfica", "RP Ecuación",
    "Tipo de disparo", "Técnica de disparo"
]

# Fila de valores (3ra fila)
row = [[
    profundidad_cima,
    profundidad_base,
    profundidad_base - profundidad_cima,
    temperatura_intervalo_centigrados_seguro,
    formacion,
    nombre_terminacion,
    tpf_estimado,
    fase,
    tipo_carga,
    f"{Sp:.4f}",
    f"{RP_grafica:.2f}",
    f"{RP_ecuacion:.2f}",
    tipo_disparo,
    "Bajo-balance"
]]

# Título de la tabla (1ra fila de la imagen)
print("\nDiseño de disparos\n")

# Mostrar tabla
print(tabulate(row, headers=headers, tablefmt="grid"))
