# ================================================================================================
# CALCULADORA DE DIFUSIÓN DE GASES - BACKEND FLASK
# ================================================================================================
# Este programa calcula coeficientes de difusión de gases siguiendo el método de Welty
# "Fundamentals of Momentum, Heat and Mass Transfer"
# 
# FÓRMULAS PRINCIPALES:
# - D_AB = 0.001858*T^(3/2)*(1/M_A + 1/M_B)^(0.5)/(P*Omega_D*sigma_AB^2)
# - Para mezclas: D_Amix = 1/(sum(y'i/D_Ai)) donde y'i = yi/(1-yA)
# ================================================================================================

from flask import Flask, render_template, request, jsonify
import math
import json

app = Flask(__name__)

# ================================================================================================
# CLASE SUBSTANCE - REPRESENTA UNA SUSTANCIA QUÍMICA
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
# CLASE BINARYMIXTURE - MANEJA CÁLCULOS PARA DOS SUSTANCIAS
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
# CLASE MIXTURE - MANEJA CÁLCULOS PARA MEZCLAS MULTICOMPONENTE
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
            print("🔄 Detectada mezcla binaria - usando cálculo directo")
            # Para mezcla binaria, crear objeto BinaryMixture y calcular directamente
            binary_mix = BinaryMixture(self.substances[0], self.substances[1], 
                                     self.temperature, self.pressure)
            return binary_mix.calculate_binary_diffusivity()
        
        # CASO 2: MEZCLA MULTICOMPONENTE (más de 2 sustancias)
        print(f"🔄 Detectada mezcla multicomponente con {len(self.substances)} sustancias")
        
        # Inicializar suma de términos 1/D_Ai ponderados
        sum_terms = 0
        
        # Obtener fracción molar de la sustancia objetivo (yA)
        target_proportion = self.target_substance.proportion
        
        print(f"📊 Sustancia objetivo: {self.target_substance.name} (fracción = {target_proportion})")
        
        # RECORRER TODAS LAS SUSTANCIAS EXCEPTO LA OBJETIVO
        for substance in self.substances:
            if substance.name != self.target_substance.name:
                print(f"🧮 Procesando: {substance.name}")
                
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
        
        print(f"📈 Suma total de términos: {sum_terms:.6e}")
        
        # PASO 4: Calcular difusividad final como 1/suma
        # D_Amix = 1 / (suma de términos)
        result = 1 / sum_terms if sum_terms > 0 else 0
        print(f"✅ Difusividad final en mezcla: {result:.6e} cm²/s")
        
        return result

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
    Ruta principal - muestra la página web con el formulario.
    
    Returns:
        str: HTML renderizado de la página principal
    """
    return render_template('index.html')

@app.route('/calculate', methods=['POST'])
def calculate_diffusivity():
    """
    Endpoint para calcular difusividad.
    
    Recibe datos JSON con:
    - substances: lista de sustancias con sus propiedades
    - target_substance: nombre de la sustancia objetivo
    - temperature, temp_unit: temperatura y unidad
    - pressure, pressure_unit: presión y unidad
    
    Returns:
        JSON: Resultado del cálculo o mensaje de error
    """
    try:
        print("🚀 Iniciando cálculo de difusividad...")
        
        # PASO 1: RECIBIR Y EXTRAER DATOS DEL REQUEST
        data = request.json
        print(f"📥 Datos recibidos: {len(data.get('substances', []))} sustancias")
        
        # Extraer datos principales
        substances_data = data['substances']           # Lista de sustancias
        target_substance_name = data['target_substance'] # Sustancia objetivo
        
        # Convertir temperatura y presión a unidades estándar
        temperature = convert_temperature_to_kelvin(data['temperature'], data['temp_unit'])
        pressure = convert_pressure_to_atm(data['pressure'], data['pressure_unit'])
        
        print(f"🌡️  Temperatura: {temperature:.2f} K")
        print(f"📊 Presión: {pressure:.3f} atm")
        print(f"🎯 Sustancia objetivo: {target_substance_name}")
        
        # PASO 2: CREAR OBJETOS SUBSTANCE PARA CADA SUSTANCIA
        substances = []      # Lista de objetos Substance
        target_substance = None  # Objeto de la sustancia objetivo
        
        for i, sub_data in enumerate(substances_data):
            print(f"🧪 Procesando sustancia {i+1}: {sub_data['name']}")
            
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
                print(f"✅ Sustancia objetivo identificada: {substance.name}")
        
        # PASO 3: VALIDACIONES DE DATOS
        
        # Verificar que se encontró la sustancia objetivo
        if not target_substance:
            error_msg = f'Sustancia objetivo "{target_substance_name}" no encontrada en la lista'
            print(f"❌ Error: {error_msg}")
            return jsonify({'error': error_msg}), 400
        
        # Verificar que las proporciones sumen 1.0
        total_proportion = sum(s.proportion for s in substances)
        print(f"📊 Suma de proporciones: {total_proportion:.6f}")
        
        if abs(total_proportion - 1.0) > 0.001:  # Tolerancia de ±0.001
            error_msg = f'Las proporciones deben sumar 1.0 (actual: {total_proportion:.3f})'
            print(f"❌ Error: {error_msg}")
            return jsonify({'error': error_msg}), 400
        
        # PASO 4: CREAR MEZCLA Y CALCULAR DIFUSIVIDAD
        print("🧮 Creando objeto Mixture y calculando difusividad...")
        
        mixture = Mixture(target_substance, substances, temperature, pressure)
        diffusivity = mixture.calculate_mixture_diffusivity()
        
        print(f"✅ Cálculo completado: {diffusivity:.6e} cm²/s")
        
        # PASO 5: PREPARAR RESPUESTA JSON
        result = {
            'diffusivity': diffusivity,                    # Valor calculado
            'units': 'cm²/s',                             # Unidades del resultado
            'temperature_k': temperature,                  # Temperatura usada (K)
            'pressure_atm': pressure,                     # Presión usada (atm)
            'target_substance': target_substance_name,     # Sustancia objetivo
            'mixture_type': 'binaria' if len(substances) == 2 else 'multicomponente'  # Tipo de mezcla
        }
        
        print("📤 Enviando resultado al cliente...")
        return jsonify(result)
        
    except Exception as e:
        # MANEJO DE ERRORES
        error_msg = str(e)
        print(f"💥 Error en el cálculo: {error_msg}")
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
    print("🌟 Iniciando servidor Flask...")
    print("🌐 La aplicación estará disponible en: http://localhost:5000")
    print("🔧 Modo debug activado - el servidor se reiniciará automáticamente al detectar cambios")
    
    app.run(debug=True)