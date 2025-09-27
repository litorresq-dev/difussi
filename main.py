# ================================================================================================
# CALCULADORA DE DIFUSIÓN DE GASES Y LÍQUIDOS - BACKEND FLASK
# ================================================================================================
# Este programa calcula coeficientes de difusión de gases siguiendo el método de Welty
# "Fundamentals of Momentum, Heat and Mass Transfer"
# Y también calcula difusividad en líquidos usando varias correlaciones
# 
# FÓRMULAS PRINCIPALES GASES:
# - D_AB = 0.001858*T^(3/2)*(1/M_A + 1/M_B)^(0.5)/(P*Omega_D*sigma_AB^2)
# - Para mezclas: D_Amix = 1/(sum(y'i/D_Ai)) donde y'i = yi/(1-yA)
#
# FÓRMULAS PRINCIPALES LÍQUIDOS:
# - Wilke-Chang: D_AB = 7.4e-8 * (φ*M_B)^0.5 * T / (μ * V_A^0.6)
# - Hayduk-Laudie: D_AB = 13.26e-5 * μ^-1.14 * V_A^-0.589
# - Scheibel: D_AB = K * T / (μ_B * V_A^1/3)
# ================================================================================================

from flask import Flask, render_template, request, jsonify
import math
import json

app = Flask(__name__)

# ================================================================================================
# TABLA DE VALORES PHI_B PARA SOLVENTES (LÍQUIDOS)
# ================================================================================================
PHI_B_TABLA = {
    "agua": 2.26,
    "metanol": 1.9,
    "etanol": 1.5,
    "benceno": 1.0,
    "eter": 1.0,
    "heptano": 1.0,
    "otros": 1.0
}

# ================================================================================================
# CLASE SUBSTANCE - REPRESENTA UNA SUSTANCIA QUÍMICA (GASES)
# ================================================================================================
class Substance():
    """
    Clase que representa una sustancia química con todas sus propiedades necesarias
    para el cálculo de difusión.
    
    Parámetros:
    - name: Nombre de la sustancia (ej: "Metano", "Etano")
    - molecular_mass: Masa molecular en g/mol (ej: 16.04 para metano)
    - critical_volume: Volumen crítico en cm³/mol (ej: 98.6 para metano)
    - critical_temperature: Temperatura crítica en K (ej: 190.6 para metano)
    - proportion: Fracción molar en la mezcla (0-1) (ej: 0.5 = 50%)
    """
    
    def __init__(self, name, molecular_mass, critical_volume, critical_temperature, proportion):
        # Propiedades básicas de la sustancia
        self.name = name                              # Nombre de la sustancia
        self.molecular_mass = molecular_mass          # Masa molecular en g/mol
        self.critical_volume = critical_volume        # Volumen crítico en cm³/mol
        self.critical_temperature = critical_temperature  # Temperatura crítica en K
        self.proportion = proportion                  # Fracción molar (0-1)
        
        # Calcular diámetro de colisión automáticamente al crear la sustancia
        self.colision_diameter = self.calculate_colision_diameter()
    
    def calculate_colision_diameter(self):
        """
        Calcula el diámetro de colisión (sigma) de la sustancia.
        
        FÓRMULA: sigma_A = 0.841 * (V_critico)^(1/3)
        
        Donde:
        - sigma_A = diámetro de colisión en Angstroms
        - V_critico = volumen crítico en cm³/mol
        - 0.841 = constante empírica
        
        Returns:
            float: Diámetro de colisión en Angstroms
        """
        return 0.841 * (self.critical_volume) ** (1/3)

# ================================================================================================
# CLASE BINARYMIXTURE - MANEJA CÁLCULOS PARA DOS SUSTANCIAS (GASES)
# ================================================================================================
class BinaryMixture():
    """
    Clase que maneja los cálculos de difusión para una mezcla de exactamente dos sustancias.
    Implementa la ecuación de Chapman-Enskog para difusión binaria.
    """
    
    def __init__(self, substance1, substance2, temperature, pressure):
        # Las dos sustancias de la mezcla binaria
        self.substance1 = substance1    # Primera sustancia (objeto Substance)
        self.substance2 = substance2    # Segunda sustancia (objeto Substance)
        self.temperature = temperature  # Temperatura del sistema en K
        self.pressure = pressure        # Presión del sistema en atm
        
        # TABLA DE INTEGRALES DE CHOQUE (T*, Omega_D)
        # Esta tabla contiene valores experimentales para la integral de choque Omega_D
        # T* = temperatura reducida = T/sqrt(T_CA * T_CB * 0.77²)
        # Omega_D = integral de choque (adimensional)
        self.table_index = [
            # [T*, Omega_D]
            [0.30, 2.662], [0.35, 2.476], [0.40, 2.318], [0.45, 2.184], [0.50, 2.066],
            [0.55, 1.966], [0.60, 1.877], [0.65, 1.798], [0.70, 1.729], [0.75, 1.667],
            [0.80, 1.612], [0.85, 1.562], [0.90, 1.517], [0.95, 1.476], [1.00, 1.439],
            [1.05, 1.406], [1.10, 1.375], [1.15, 1.346], [1.20, 1.320], [1.25, 1.296],
            [1.30, 1.273], [1.35, 1.253], [1.40, 1.233], [1.45, 1.215], [1.50, 1.198],
            [1.55, 1.182], [1.60, 1.167], [1.65, 1.153], [1.70, 1.140], [1.80, 1.116],
            [1.85, 1.105], [1.90, 1.094], [1.95, 1.084], [2.00, 1.075], [2.10, 1.057],
            [2.20, 1.041], [2.30, 1.026], [2.40, 1.012], [2.50, 0.9996], [2.60, 0.9878],
            [2.70, 0.9770], [2.80, 0.9672], [2.90, 0.9576], [3.00, 0.9490], [3.10, 0.9406],
            [3.20, 0.9328], [3.30, 0.9256], [3.40, 0.9186], [3.50, 0.9120], [3.60, 0.9058],
            [3.70, 0.8998], [3.80, 0.8942], [3.90, 0.8888], [4.00, 0.8836], [4.10, 0.8788],
            [4.20, 0.8740], [4.30, 0.8694], [4.40, 0.8652], [4.50, 0.8610], [4.60, 0.8568],
            [4.70, 0.8530], [4.80, 0.8492], [4.90, 0.8456], [5.0, 0.8422], [6.0, 0.8124],
            [7.0, 0.7896], [8.0, 0.7712], [10.0, 0.7424], [20.0, 0.6640], [30.0, 0.6232],
            [40.0, 0.5960], [50.0, 0.5756], [60.0, 0.5596], [70.0, 0.5464], [80.0, 0.5352],
            [90.0, 0.5256]
        ]
    
    def calculate_colision_diameter(self):
        """
        Calcula el diámetro de colisión promedio para la mezcla binaria.
        
        FÓRMULA: sigma_AB = (sigma_A + sigma_B) / 2
        
        Es simplemente el promedio aritmético de los diámetros de las dos sustancias.
        
        Returns:
            float: Diámetro de colisión promedio en Angstroms
        """
        return (self.substance1.colision_diameter + self.substance2.colision_diameter) / 2
    
    def calculate_difusivity_integral(self):
        """
        Calcula la integral de difusión Omega_D usando interpolación en la tabla.
        
        PROCESO:
        1. Calcula T* = T / sqrt(T_CA * T_CB * 0.77²)
        2. Busca T* en la tabla de valores experimentales
        3. Si no encuentra valor exacto, interpola linealmente entre valores cercanos
        
        Returns:
            float: Valor de Omega_D (integral de choque)
        """
        # PASO 1: Calcular temperatura reducida T*
        # T* = T / sqrt(T_critica_A * T_critica_B * 0.77²)
        # El factor 0.77 es una constante empírica del método de Welty
        t_star = self.temperature / math.sqrt(
            self.substance1.critical_temperature * 
            self.substance2.critical_temperature * 
            (0.77 ** 2)
        )
        
        # PASO 2: Buscar T* en la tabla e interpolar si es necesario
        for i in range(len(self.table_index)):
            # Si encontramos un valor de T* mayor o igual al calculado
            if self.table_index[i][0] >= t_star:
                # Si es el primer valor de la tabla, usar ese valor directamente
                if i == 0:
                    return self.table_index[0][1]
                else:
                    # INTERPOLACIÓN LINEAL entre dos puntos de la tabla
                    # Puntos: (x1, y1) y (x2, y2)
                    x1, y1 = self.table_index[i-1]  # Punto anterior
                    x2, y2 = self.table_index[i]    # Punto actual
                    
                    # Fórmula de interpolación lineal: y = y1 + (y2-y1)*(x-x1)/(x2-x1)
                    omega_d = y1 + (y2 - y1) * (t_star - x1) / (x2 - x1)
                    return omega_d
        
        # Si T* es mayor que todos los valores en la tabla, usar el último valor
        return self.table_index[-1][1]
    
    def calculate_binary_diffusivity(self):
        """
        Calcula la difusividad binaria usando la ecuación de Chapman-Enskog.
        
        ECUACIÓN PRINCIPAL:
        D_AB = 0.001858 * T^(3/2) * (1/M_A + 1/M_B)^(0.5) / (P * Omega_D * sigma_AB²)
        
        Donde:
        - D_AB = coeficiente de difusión en cm²/s
        - T = temperatura en K
        - M_A, M_B = masas moleculares en g/mol
        - P = presión en atm
        - Omega_D = integral de choque (de la tabla)
        - sigma_AB = diámetro de colisión promedio en Angstroms
        - 0.001858 = constante de conversión de unidades
        
        Returns:
            float: Coeficiente de difusión binaria en cm²/s
        """
        # PASO 1: Obtener integral de choque de la tabla
        omega_d = self.calculate_difusivity_integral()
        
        # PASO 2: Obtener diámetro de colisión promedio
        o_ab = self.calculate_colision_diameter()
        
        # PASO 3: Calcular cada término de la ecuación
        
        # Término 1: 0.001858 * T^(3/2)
        # La potencia 1.5 equivale a (3/2)
        term1 = 0.001858 * (self.temperature ** 1.5)
        
        # Término 2: (1/M_A + 1/M_B)^(0.5)
        # Suma de inversos de masas moleculares, elevado a 0.5
        term2 = (1/self.substance1.molecular_mass + 1/self.substance2.molecular_mass)**(0.5)
        
        # Término 3: P * Omega_D * sigma_AB²
        # Denominador de la ecuación
        term3 = self.pressure * omega_d * (o_ab ** 2)
        
        # PASO 4: Calcular difusividad final
        # D_AB = (término1 * término2) / término3
        d_ab = term1 * term2 / term3
        
        return d_ab

# ================================================================================================
# CLASE MIXTURE - MANEJA CÁLCULOS PARA MEZCLAS MULTICOMPONENTE (GASES)
# ================================================================================================
class Mixture():
    """
    Clase que maneja cálculos de difusión para mezclas de múltiples sustancias.
    
    Para mezclas con más de 2 componentes, usa la aproximación:
    D_Amix = 1 / (suma de y'i/D_Ai)
    
    Donde:
    - D_Amix = difusividad de A en la mezcla
    - y'i = fracción molar corregida de componente i
    - D_Ai = difusividad binaria de A con componente i
    """
    
    def __init__(self, target_substance, substances, temperature, pressure):
        self.substances = substances        # Lista completa de todas las sustancias
        self.target_substance = target_substance  # Sustancia objetivo (A) para calcular su difusividad
        self.temperature = temperature      # Temperatura del sistema en K
        self.pressure = pressure           # Presión del sistema en atm
    
    def calculate_mixture_diffusivity(self):
        """
        Calcula la difusividad de la sustancia objetivo en la mezcla multicomponente.
        
        CASOS:
        1. Si hay solo 2 sustancias -> Usar cálculo binario directo
        2. Si hay más de 2 sustancias -> Usar aproximación multicomponente
        
        MÉTODO MULTICOMPONENTE:
        D_Amix = 1 / (suma de términos y'i/D_Ai)
        
        Donde:
        - y'i = yi / (1 - yA)  (fracción molar corregida)
        - yi = fracción molar del componente i
        - yA = fracción molar de la sustancia objetivo A
        - D_Ai = difusividad binaria entre A y componente i
        
        Returns:
            float: Coeficiente de difusión en la mezcla en cm²/s
        """
        
        # CASO 1: MEZCLA BINARIA (solo 2 sustancias)
        if len(self.substances) == 2:
            print("Detectada mezcla binaria - usando cálculo directo")
            # Para mezcla binaria, crear objeto BinaryMixture y calcular directamente
            binary_mix = BinaryMixture(self.substances[0], self.substances[1], 
                                     self.temperature, self.pressure)
            return binary_mix.calculate_binary_diffusivity()
        
        # CASO 2: MEZCLA MULTICOMPONENTE (más de 2 sustancias)
        print(f"Detectada mezcla multicomponente con {len(self.substances)} sustancias")
        
        # Inicializar suma de términos 1/D_Ai ponderados
        sum_terms = 0
        
        # Obtener fracción molar de la sustancia objetivo (yA)
        target_proportion = self.target_substance.proportion
        
        print(f"Sustancia objetivo: {self.target_substance.name} (fracción = {target_proportion})")
        
        # RECORRER TODAS LAS SUSTANCIAS EXCEPTO LA OBJETIVO
        for substance in self.substances:
            if substance.name != self.target_substance.name:
                print(f"Procesando: {substance.name}")
                
                # PASO 1: Calcular fracción molar corregida y'X
                # y'X = yX / (1 - yA)
                # Esto "reescala" las fracciones molares excluyendo la sustancia objetivo
                y_prime = substance.proportion / (1 - target_proportion)
                print(f"   - Fracción original (y): {substance.proportion}")
                print(f"   - Fracción corregida (y'): {y_prime}")
                
                # PASO 2: Calcular difusividad binaria D_AX entre objetivo y esta sustancia
                binary_mix = BinaryMixture(self.target_substance, substance, 
                                         self.temperature, self.pressure)
                d_ax = binary_mix.calculate_binary_diffusivity()
                print(f"   - Difusividad binaria D_A{substance.name}: {d_ax:.6e} cm²/s")
                
                # PASO 3: Agregar término y'X/D_AX a la suma
                term = y_prime / d_ax
                sum_terms += term
                print(f"   - Término y'/D: {term:.6e}")
        
        print(f"Suma total de términos: {sum_terms:.6e}")
        
        # PASO 4: Calcular difusividad final como 1/suma
        # D_Amix = 1 / (suma de términos)
        result = 1 / sum_terms if sum_terms > 0 else 0
        print(f"Difusividad final en mezcla: {result:.6e} cm²/s")
        
        return result

# ================================================================================================
# FUNCIONES DE DIFUSIVIDAD EN LÍQUIDOS
# ================================================================================================

def calculate_Va_from_Vc(Vc):
    """
    Calcula volumen molar Va usando correlación de Tyn y Calus.
    FÓRMULA: Va = 0.285 * (Vc^1.048)
    """
    return 0.285 * (Vc ** 1.048)

def wilke_chang(T, mu, Va, M, phi):
    """
    Ecuación de Wilke-Chang para difusividad en líquidos.
    D_AB = 7.4e-8 * (φ*M_B)^0.5 * T / (μ * V_A^0.6)
    
    Args:
        T: Temperatura en K
        mu: Viscosidad del solvente en cP
        Va: Volumen molar del soluto en cm³/mol
        M: Masa molecular del solvente en g/mol
        phi: Factor de asociación del solvente
    
    Returns:
        float: Difusividad en cm²/s
    """
    return 7.4e-8 * (phi * M)**0.5 * T / (mu * (Va**0.6))

def hayduk_laudie(mu, Va):
    """
    Ecuación de Hayduk-Laudie para difusividad en líquidos.
    D_AB = 13.26e-5 * μ^-1.14 * V_A^-0.589
    
    Args:
        mu: Viscosidad del solvente en cP
        Va: Volumen molar del soluto en cm³/mol
    
    Returns:
        float: Difusividad en cm²/s
    """
    return 13.26e-5 * (mu**-1.14) * (Va**-0.589)

def scheibel(T, muB, Va, Vb):
    """
    Ecuación de Scheibel para difusividad en líquidos.
    D_AB = K * T / (μ_B * V_A^1/3)
    
    Donde K depende de la relación Va/Vb
    
    Args:
        T: Temperatura en K
        muB: Viscosidad del solvente en cP
        Va: Volumen molar del soluto en cm³/mol
        Vb: Volumen molar del solvente en cm³/mol
    
    Returns:
        float: Difusividad en cm²/s
    """
    # Calcular la constante K según las condiciones
    ratio = Va / Vb if Vb > 0 else float('inf')
    
    if Va < 2*Vb:
        K = 18.9e-8
    elif Va < 2.5*Vb:
        K = 17.5e-8
    else:
        K = 8.2e-8 * (1 + (3 * Vb / Va)**(2/3))
    
    return K * T / (muB * Va**(1/3))

def tyne_extrapolation(DAB_T2, T1, T2, Tc, deltaHv):
    """
    Extrapolación de difusividad usando método de Tyne.
    D_AB(T1) = D_AB(T2) * ((Tc - T2)/(Tc - T1))^n
    
    Args:
        DAB_T2: Difusividad conocida a temperatura T2 en cm²/s
        T1: Temperatura objetivo en K
        T2: Temperatura conocida en K
        Tc: Temperatura crítica en K
        deltaHv: Entalpía de vaporización en kJ/kmol
    
    Returns:
        tuple: (Difusividad calculada, exponente n usado)
    """
    # Determinar el exponente n basado en deltaHv
    if 7900 <= deltaHv < 30000:
        n = 3
    elif 30000 <= deltaHv < 39000:
        n = 4
    elif 39000 <= deltaHv < 46000:
        n = 6
    elif 46000 <= deltaHv < 50000:
        n = 8
    else:  # >= 50000
        n = 10

    D_AB_T1 = DAB_T2 * ((Tc - T2)/(Tc - T1))**n
    return D_AB_T1, n

# ================================================================================================
# FUNCIONES DE CONVERSIÓN DE UNIDADES
# ================================================================================================

def convert_temperature_to_kelvin(temp, unit):
    """
    Convierte temperatura de cualquier unidad a Kelvin.
    
    CONVERSIONES:
    - Celsius a K: K = °C + 273.15
    - Fahrenheit a K: K = (°F - 32) * 5/9 + 273.15
    - Kelvin: sin conversión
    
    Args:
        temp (float): Valor de temperatura
        unit (str): Unidad ('C', 'F', 'K')
    
    Returns:
        float: Temperatura en Kelvin
    """
    if unit.lower() == 'c' or unit.lower() == 'celsius':
        # Celsius a Kelvin: sumar 273.15
        return temp + 273.15
    elif unit.lower() == 'f' or unit.lower() == 'fahrenheit':
        # Fahrenheit a Kelvin: convertir a Celsius primero, luego a Kelvin
        return (temp - 32) * 5/9 + 273.15
    else:  # Kelvin (sin conversión)
        return temp

def convert_pressure_to_atm(pressure, unit):
    """
    Convierte presión de cualquier unidad a atmósferas.
    
    FACTORES DE CONVERSIÓN:
    - 1 atm = 1 atm (sin conversión)
    - 1 bar = 0.986923 atm
    - 1 psi = 0.068046 atm
    - 1 Pa = 0.00000986923 atm
    - 1 kPa = 0.00986923 atm
    - 1 mmHg = 0.00131579 atm
    - 1 torr = 0.00131579 atm
    
    Args:
        pressure (float): Valor de presión
        unit (str): Unidad de presión
    
    Returns:
        float: Presión en atmósferas
    """
    conversions = {
        'atm': 1,           # Sin conversión
        'bar': 0.986923,    # 1 bar = 0.986923 atm
        'psi': 0.068046,    # 1 psi = 0.068046 atm
        'pa': 0.00000986923, # 1 Pa = 0.00000986923 atm
        'kpa': 0.00986923,  # 1 kPa = 0.00986923 atm
        'mmhg': 0.00131579, # 1 mmHg = 0.00131579 atm
        'torr': 0.00131579  # 1 torr = 0.00131579 atm
    }
    
    # Multiplicar por el factor de conversión correspondiente
    return pressure * conversions.get(unit.lower(), 1)

# ================================================================================================
# RUTAS DE FLASK - ENDPOINTS DE LA APLICACIÓN WEB
# ================================================================================================

@app.route('/')
def index():
    """
    Ruta principal - muestra la página web con el formulario para gases.
    
    Returns:
        str: HTML renderizado de la página principal
    """
    return render_template('index.html')

@app.route('/liquids')
def liquids():
    """
    Ruta para cálculo de difusividad en líquidos.
    
    Returns:
        str: HTML renderizado de la página de líquidos
    """
    return render_template('index2.html')

@app.route('/calculate', methods=['POST'])
def calculate_diffusivity():
    """
    Endpoint para calcular difusividad en gases.
    
    Recibe datos JSON con:
    - substances: lista de sustancias con sus propiedades
    - target_substance: nombre de la sustancia objetivo
    - temperature, temp_unit: temperatura y unidad
    - pressure, pressure_unit: presión y unidad
    
    Returns:
        JSON: Resultado del cálculo o mensaje de error
    """
    try:
        print("Iniciando cálculo de difusividad en gases...")
        
        # PASO 1: RECIBIR Y EXTRAER DATOS DEL REQUEST
        data = request.json
        print(f"Datos recibidos: {len(data.get('substances', []))} sustancias")
        
        # Extraer datos principales
        substances_data = data['substances']           # Lista de sustancias
        target_substance_name = data['target_substance'] # Sustancia objetivo
        
        # Convertir temperatura y presión a unidades estándar
        temperature = convert_temperature_to_kelvin(data['temperature'], data['temp_unit'])
        pressure = convert_pressure_to_atm(data['pressure'], data['pressure_unit'])
        
        print(f"Temperatura: {temperature:.2f} K")
        print(f"Presión: {pressure:.3f} atm")
        print(f"Sustancia objetivo: {target_substance_name}")
        
        # PASO 2: CREAR OBJETOS SUBSTANCE PARA CADA SUSTANCIA
        substances = []      # Lista de objetos Substance
        target_substance = None  # Objeto de la sustancia objetivo
        
        for i, sub_data in enumerate(substances_data):
            print(f"Procesando sustancia {i+1}: {sub_data['name']}")
            
            # Crear objeto Substance con datos convertidos
            substance = Substance(
                name=sub_data['name'],
                molecular_mass=sub_data['molecular_mass'],
                critical_volume=sub_data['critical_volume'],
                critical_temperature=convert_temperature_to_kelvin(
                    sub_data['critical_temperature'], 
                    sub_data['temp_unit']
                ),
                proportion=sub_data['proportion']
            )
            
            # Agregar a la lista
            substances.append(substance)
            
            # Si es la sustancia objetivo, guardar referencia
            if substance.name == target_substance_name:
                target_substance = substance
                print(f"Sustancia objetivo identificada: {substance.name}")
        
        # PASO 3: VALIDACIONES DE DATOS
        
        # Verificar que se encontró la sustancia objetivo
        if not target_substance:
            error_msg = f'Sustancia objetivo "{target_substance_name}" no encontrada en la lista'
            print(f"Error: {error_msg}")
            return jsonify({'error': error_msg}), 400
        
        # Verificar que las proporciones sumen 1.0
        total_proportion = sum(s.proportion for s in substances)
        print(f"Suma de proporciones: {total_proportion:.6f}")
        
        if abs(total_proportion - 1.0) > 0.001:  # Tolerancia de ±0.001
            error_msg = f'Las proporciones deben sumar 1.0 (actual: {total_proportion:.3f})'
            print(f"Error: {error_msg}")
            return jsonify({'error': error_msg}), 400
        
        # PASO 4: CREAR MEZCLA Y CALCULAR DIFUSIVIDAD
        print("Creando objeto Mixture y calculando difusividad...")
        
        mixture = Mixture(target_substance, substances, temperature, pressure)
        diffusivity = mixture.calculate_mixture_diffusivity()
        
        print(f"Cálculo completado: {diffusivity:.6e} cm²/s")
        
        # PASO 5: PREPARAR RESPUESTA JSON
        result = {
            'diffusivity': diffusivity,                    # Valor calculado
            'units': 'cm²/s',                             # Unidades del resultado
            'temperature_k': temperature,                  # Temperatura usada (K)
            'pressure_atm': pressure,                     # Presión usada (atm)
            'target_substance': target_substance_name,     # Sustancia objetivo
            'mixture_type': 'binaria' if len(substances) == 2 else 'multicomponente'  # Tipo de mezcla
        }
        
        print("Enviando resultado al cliente...")
        return jsonify(result)
        
    except Exception as e:
        # MANEJO DE ERRORES
        error_msg = str(e)
        print(f"Error en el cálculo: {error_msg}")
        return jsonify({'error': error_msg}), 500

@app.route('/calculate-liquid', methods=['POST'])
def calculate_liquid_diffusivity():
    """
    Endpoint para calcular difusividad en líquidos.
    
    Recibe datos JSON con:
    - method: método de cálculo ('wilke_chang', 'hayduk_laudie', 'scheibel', 'tyne')
    - temperature, temp_unit: temperatura y unidad
    - viscosity: viscosidad del solvente en cP
    - molecular_mass_solvent: masa molecular del solvente en g/mol (para Wilke-Chang)
    - solvent: nombre del solvente (para obtener phi_B)
    - molar_volume: volumen molar del soluto en cm³/mol
    - critical_volume: volumen crítico en cm³/mol (si no se tiene Va)
    - solvent_molar_volume: volumen molar del solvente (para Scheibel)
    - solvent_critical_volume: volumen crítico del solvente (para Scheibel)
    - Para Tyne: DAB_T2, T1, T2, Tc, deltaHv
    
    Returns:
        JSON: Resultado del cálculo o mensaje de error
    """
    try:
        print("Iniciando cálculo de difusividad en líquidos...")
        
        # PASO 1: RECIBIR Y EXTRAER DATOS DEL REQUEST
        data = request.json
        method = data['method']
        
        print(f"Método seleccionado: {method}")
        
        # Convertir temperatura a Kelvin
        temperature = convert_temperature_to_kelvin(data['temperature'], data['temp_unit'])
        print(f"Temperatura: {temperature:.2f} K")
        
        # Variables comunes
        result = {
            'temperature_k': temperature,
            'method': method,
            'units': 'cm²/s'
        }
        
        if method == 'wilke_chang':
            # MÉTODO WILKE-CHANG
            viscosity = data['viscosity']  # cP
            molecular_mass_solvent = data['molecular_mass_solvent']  # g/mol
            solvent = data['solvent'].lower()
            
            # Obtener factor de asociación phi_B
            phi_B = PHI_B_TABLA.get(solvent, 1.0)
            print(f"Factor de asociación φ_B para {solvent}: {phi_B}")
            
            # Obtener volumen molar del soluto
            if 'molar_volume' in data and data['molar_volume']:
                Va = data['molar_volume']
                print(f"Usando volumen molar directo: {Va} cm³/mol")
            elif 'critical_volume' in data and data['critical_volume']:
                Va = calculate_Va_from_Vc(data['critical_volume'])
                print(f"Calculado Va desde Vc usando Tyn-Calus: {Va:.2f} cm³/mol")
            else:
                raise ValueError("Se requiere volumen molar (Va) o volumen crítico (Vc) del soluto")
            
            # Calcular difusividad
            diffusivity = wilke_chang(temperature, viscosity, Va, molecular_mass_solvent, phi_B)
            
            result.update({
                'diffusivity': diffusivity,
                'viscosity': viscosity,
                'molecular_mass_solvent': molecular_mass_solvent,
                'solvent': data['solvent'],
                'phi_B': phi_B,
                'Va': Va
            })
            
        elif method == 'hayduk_laudie':
            # MÉTODO HAYDUK-LAUDIE
            viscosity = data['viscosity']  # cP
            
            # Obtener volumen molar del soluto
            if 'molar_volume' in data and data['molar_volume']:
                Va = data['molar_volume']
                print(f"Usando volumen molar directo: {Va} cm³/mol")
            elif 'critical_volume' in data and data['critical_volume']:
                Va = calculate_Va_from_Vc(data['critical_volume'])
                print(f"Calculado Va desde Vc usando Tyn-Calus: {Va:.2f} cm³/mol")
            else:
                raise ValueError("Se requiere volumen molar (Va) o volumen crítico (Vc) del soluto")
            
            # Calcular difusividad
            diffusivity = hayduk_laudie(viscosity, Va)
            
            result.update({
                'diffusivity': diffusivity,
                'viscosity': viscosity,
                'Va': Va
            })
            
        elif method == 'scheibel':
            # MÉTODO SCHEIBEL
            viscosity = data['viscosity']  # cP
            
            # Obtener volumen molar del soluto
            if 'molar_volume' in data and data['molar_volume']:
                Va = data['molar_volume']
                print(f"Usando volumen molar del soluto directo: {Va} cm³/mol")
            elif 'critical_volume' in data and data['critical_volume']:
                Va = calculate_Va_from_Vc(data['critical_volume'])
                print(f"Calculado Va desde Vc usando Tyn-Calus: {Va:.2f} cm³/mol")
            else:
                raise ValueError("Se requiere volumen molar (Va) o volumen crítico (Vc) del soluto")
            
            # Obtener volumen molar del solvente
            if 'solvent_molar_volume' in data and data['solvent_molar_volume']:
                Vb = data['solvent_molar_volume']
                print(f"Usando volumen molar del solvente directo: {Vb} cm³/mol")
            elif 'solvent_critical_volume' in data and data['solvent_critical_volume']:
                Vb = calculate_Va_from_Vc(data['solvent_critical_volume'])
                print(f"Calculado Vb desde Vc del solvente usando Tyn-Calus: {Vb:.2f} cm³/mol")
            else:
                raise ValueError("Se requiere volumen molar (Vb) o volumen crítico (Vc) del solvente")
            
            # Calcular difusividad
            diffusivity = scheibel(temperature, viscosity, Va, Vb)
            
            result.update({
                'diffusivity': diffusivity,
                'viscosity': viscosity,
                'Va': Va,
                'Vb': Vb
            })
            
        elif method == 'tyne':
            # MÉTODO TYNE (EXTRAPOLACIÓN)
            DAB_T2 = data['DAB_T2']  # cm²/s
            T1 = temperature  # K (temperatura objetivo)
            T2 = convert_temperature_to_kelvin(data['T2'], data['T2_unit'])  # K (temperatura conocida)
            Tc = convert_temperature_to_kelvin(data['Tc'], data['Tc_unit'])  # K (temperatura crítica)
            deltaHv = data['deltaHv']  # kJ/kmol
            
            print(f"T1 (objetivo): {T1:.2f} K")
            print(f"T2 (conocida): {T2:.2f} K")
            print(f"Tc (crítica): {Tc:.2f} K")
            print(f"ΔHv: {deltaHv} kJ/kmol")
            
            # Calcular difusividad y exponente
            diffusivity, n = tyne_extrapolation(DAB_T2, T1, T2, Tc, deltaHv)
            
            result.update({
                'diffusivity': diffusivity,
                'DAB_T2': DAB_T2,
                'T1': T1,
                'T2': T2,
                'Tc': Tc,
                'deltaHv': deltaHv,
                'exponent_n': n
            })
            
        else:
            raise ValueError(f"Método no reconocido: {method}")
        
        print(f"Cálculo completado: {diffusivity:.6e} cm²/s usando {method}")
        return jsonify(result)
        
    except Exception as e:
        # MANEJO DE ERRORES
        error_msg = str(e)
        print(f"Error en el cálculo: {error_msg}")
        return jsonify({'error': error_msg}), 500

# ================================================================================================
# PUNTO DE ENTRADA PRINCIPAL
# ================================================================================================
if __name__ == '__main__':
    """
    Ejecuta la aplicación Flask en modo de desarrollo.
    
    La aplicación estará disponible en: http://localhost:5000
    El modo debug=True permite:
    - Recarga automática cuando se modifica el código
    - Información detallada de errores en el navegador
    - Mensajes de debug en la consola
    """
    print("Iniciando servidor Flask...")
    print("La aplicación estará disponible en: http://localhost:5000")
    print("Rutas disponibles:")
    print("  - /          (Cálculo de difusividad en gases)")
    print("  - /liquids   (Cálculo de difusividad en líquidos)")
    print("Modo debug activado - el servidor se reiniciará automáticamente al detectar cambios")
    
    app.run(debug=True)