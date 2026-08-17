# Proyecto Integrador III: Modelo de Saturación con Olvido
## Cátedra de Métodos Numéricos Avanzados (3.° Año)
### Ingeniería en Inteligencia Artificial — UNSTA
**Autores del Equipo 3:** barrios, Julieta L. · Sappia, Lucio A. · Tartarini, Luciano A.  
**Docente Coordinador:** Cátedra de Métodos Numéricos  
**Fecha de Presentación:** 23 de junio de 2026

---

## 1. Introducción y Justificación del Modelo

En la teoría de difusión social, propagación de rumores y mercadotecnia viral, los modelos matemáticos basados en Ecuaciones Diferenciales Ordinarias (EDO) desempeñan un rol fundamental. Este simulador implementa el **Modelo de Saturación con Olvido**, una variante no lineal inspirada en los modelos epidemiológicos clásicos (como el modelo SIS) y en los modelos de crecimiento logístico (Verhulst), adaptada para representar la dinámica de difusión de información dentro de una red de contactos cerrados.

El modelo asume una población total constante $P$ donde cada individuo puede estar en uno de dos estados respecto a un mensaje o idea:
1. **Desinformado / Receptivo:** Individuos que no conocen la información pero son susceptibles de adquirirla mediante contacto social.
2. **Informado / Activo:** Individuos que conocen la información, la difunden activamente y, tras un tiempo medio, sufren un desinterés u olvido regresando al grupo receptivo.

### 1.1 Ecuación Diferencial Gobernadora

La evolución temporal del número de personas informadas a tiempo $t$, denotado por $y(t)$, se rige por la siguiente EDO de primer orden no lineal:

$$\frac{dy}{dt} = \beta \cdot y(t) \cdot [P - y(t)] - \delta \cdot y(t)$$

Donde los parámetros representan:
*   $y(t)$: Número de personas informadas activas en el instante $t$.
*   $P$: Población total expuesta (fijada en $P = 100,000$ agentes).
*   $P - y(t)$: Fracción o número de personas desinformadas receptivas.
*   $\beta$: Coeficiente o tasa de contagio social ($1/\text{personas}\cdot\text{día}$). Regula la probabilidad de transmisión de la información por unidad de tiempo al interactuar una persona informada con una desinformada.
*   $\delta$: Coeficiente o tasa de olvido ($1/\text{día}$). Representa la velocidad con la que los agentes informados pierden el interés en la idea o rumor y regresan a ser durmientes o receptivos.

---

## 2. Análisis Cualitativo del Modelo

La ecuación diferencial gobernadora es autónoma y se puede reescribir analíticamente factorizando el término $y$:

$$\frac{dy}{dt} = y \cdot \left[ (\beta P - \delta) - \beta y \right]$$

### 2.1 Puntos de Equilibrio (Estados Estacionarios)

Para determinar los estados estacionarios de la red, resolvemos $\frac{dy}{dt} = 0$:

1.  **Equilibrio Libre de Información (Extinción):**  
    $$y_1^* = 0$$
    Siempre es un punto crítico. Representa una población donde nadie difunde el mensaje.

2.  **Equilibrio de Coexistencia Asintótica (Saturación de Red):**  
    $$y_2^* = P - \frac{\delta}{\beta}$$
    Este equilibrio físicamente realista tiene sentido únicamente si $P > \frac{\delta}{\beta}$, lo cual requiere que la tasa de contagio social sea mayor al cociente del olvido respecto a la población.

### 2.2 Número Reproductivo Básico ($R_0$)

Por analogía con los modelos matemáticos en epidemiología, definimos el **Número Reproductivo Básico ($R_0$)** como:

$$R_0 = \frac{\beta \cdot P}{\delta}$$

Este parámetro adimensional categoriza el comportamiento a largo plazo de la solución ($t \to \infty$):
*   **Caso $R_0 \le 1$ ($\beta P \le \delta$):**  
    El punto de equilibrio $y_1^* = 0$ es local y globalmente asintóticamente estable en el semiplano positivo. Físicamente, el olvido domina sobre el contagio social de contactos; la información carece de fuerza de autoreproducción y se extingue exponencialmente de la población.
*   **Caso $R_0 > 1$ ($\beta P > \delta$):**  
    El punto $y_1^* = 0$ se vuelve inestable (repulsor), y el equilibrio de persistencia estocástica/saturación $y_2^* = K^* = P - \frac{\delta}{\beta}$ es asintóticamente estable (atractor). El sistema converge asintóticamente hacia esta **Capacidad de Carga Efectiva ($K^*$)**, manteniendo un flujo constante de personas informadas y desinformadas en equilibrio dinámico.

---

## 3. Soluciones Numéricas Implementadas

Dado que el modelo aborda una ecuación no lineal donde se evalúan comportamientos dinámicos bajo mallas temporales discretas de paso de integración $h$, el simulador expone tres metodologías computacionales para comparar precisión, estabilidad y convergencia.

### 3.1 Método de Runge-Kutta de 4.° Orden (RK4)

Es un método de un solo paso de alta precisión, clasificado como de orden de convergencia global $O(h^4)$. Evalúa la pendiente de la función en cuatro puntos estratégicos del intervalo $[t_n, t_n + h]$ para minimizar el error de truncamiento local.

Para el problema de valor inicial $y' = f(t, y)$ con $y(t_n) = y_n$, las ecuaciones de recurrencia son:

$$k_1 = f(t_n, y_n)$$
$$k_2 = f\left(t_n + \frac{h}{2}, y_n + \frac{h}{2} k_1\right)$$
$$k_3 = f\left(t_n + \frac{h}{2}, y_n + \frac{h}{2} k_2\right)$$
$$k_4 = f(t_n + h, y_n + h k_3)$$
$$y_{n+1} = y_n + \frac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4)$$

*   **Propiedades:** Excelente estabilidad numérica absoluta. Actúa como el solver de referencia ("Ground Truth") en nuestra interfaz gráfica.

### 3.2 Método de Euler Explícito (Primer Orden)

El método más básico para integrar una EDO. Se fundamenta en la aproximación de la derivada mediante la diferencia finita hacia adelante (primer término de la serie de Taylor). Su orden de convergencia global es $O(h)$.

$$y_{n+1} = y_n + h \cdot f(t_n, y_n)$$

*   **Propiedades:** Es de cómputo ligero pero propenso a oscilaciones numéricas divergentes o errores cumulativos grotescos si el paso $h$ es grande o si la derivada es empinada (rigidez numérica o *stiffness*).

### 3.3 Método Multipaso de Adams-Bashforth de 2 Pasos (AB2)

Es un integrador explícito multipaso lineal con error de convergencia global de segundo orden $O(h^2)$. A diferencia de RK4, aprovecha la información de pasos ya calculados anteriormente para reconstruir un polinomio de interpolación de la derivada.

La fórmula de recurrencia utiliza los valores evaluados en $t_n$ y $t_{n-1}$:

$$y_{n+1} = y_n + h \left[ \frac{3}{2} f(t_n, y_n) - \frac{1}{2} f(t_{n-1}, y_{n-1}) \right]$$

*   **Arranque (Bootstrap):** Como el método requiere de dos puntos iniciales y solo contamos con $y(0) = y_0$, el primer paso $y_1$ se inicializa estrictamente usando un paso previo de Runge-Kutta de cuarto orden para evitar degradar la precisión global del algoritmo.

---

## 4. Método del Disparo No Lineal (Shooting Method) adaptado

Uno de los aportes prácticos de mayor relieve matemático de este simulador es el módulo de **Optimización por Método del Disparo**. 

Normalmente, el problema de valor inicial parte de condiciones conocidas a $t=0$. Sin embargo, en la planificación de campañas publicitarias es común enfrentar un **problema de contorno** o meta intermedia: *«Queremos hallar con exactitud qué tasa de contagio $\beta$ requiere la red para asegurar que en un día exacto $t_{meta}$ haya exactamente un volumen de informados igual a $y_{meta}$»*.

### 4.1 Formulación Matemática

Consiste en definir una función residual no lineal $G(\beta)$ que dependa paramétricamente de $\beta$:

$$G(\beta) = y(t_{meta}; \beta) - y_{meta} = 0$$

Donde $y(t_{meta}; \beta)$ representa la integración numérica por método RK4 partiendo de $y(0) = 100$ evaluada en el extremo tiempo $t_{meta}$ usando como coeficiente de contagio el valor de $\beta$ ensayado.

### 4.2 Algoritmo de Bisección sobre el Espacio de Parámetros

El simulador implementa el método clásico de bisección para resolver de forma segura y robusta la raíz $G(\beta) = 0$.

1.  **Definición del Intervalo Inicial:** Se establece una cota de búsqueda físicamente coherente para $\beta$: $[a, b] = [1\times 10^{-7}, 3\times 10^{-5}]$.
2.  **Verificación del Teorema de Bolzano:** Se evalúan $G(a)$ y $G(b)$. Si tienen signos opuestos ($G(a) \cdot G(b) < 0$), se garantiza la existencia de al menos una tasa óptima $\beta^*$.
3.  **Bucle de Bisección:**
    *   Se calcula el punto medio: $c = \frac{a+b}{2}$.
    *   Se efectúa una simulación completa con el estimador actual $\beta = c$ bajo el solver RK4 y paso ajustado $h$ desde $t=0$ hasta $t_{meta}$.
    *   Se extrae $y(t_{meta}; c)$ y se calcula el residuo o error absoluto respecto a la meta: $\text{Residuo} = | y(t_{meta}; c) - y_{meta} |$.
    *   Si $\text{Residuo} \le \epsilon$ (tolerancia absoluta seleccionada por el usuario), el algoritmo converge exitosamente regresando $\beta^* = c$.
    *   De lo contrario, determinamos el nuevo subintervalo basándonos en los signos de $G(a)$ y $G(c)$ para repetir la subdivisión.
4.  **Límite de Operación:** Se impidió el estancamiento limitando el proceso a un máximo de 50 iteraciones (suficientes para lograr precisión de hasta $10^{-16}$ en el flotante).

---

## 5. Arquitectura del Software y Componentes del Panel

El simulador está construido sobre una arquitectura web de alto rendimiento modular en **React + TypeScript (Vite)** impulsada estilísticamente por **Tailwind CSS** en un entorno inmersivo oscuro "Slate & Cyan" de grado profesional para uso en laboratorio de ciencias. 

Los componentes interactivos se desglosan de la siguiente forma:

1.  **Panel de Parámetros Físicos:** Permite manipular libremente el coeficiente de transferencia por contagio $\beta$, el horizonte máximo $T$ de simulación, la tasa de olvido $\delta$, y elegir con precisión el paso de discretización $h \in \{0.2, 0.05, 0.01\}$ para examinar cómo cambian los errores de aproximación.
2.  **Monitor de Red Social Interactiva (SVG Dinámico):** Mapea algebraicamente un grafo radial de 150 nodos activos estructurado con conexiones internas fijas. A medida que progresa el tiempo $t_{current}$, la porción correspondiente de nodos se enciende o se apaga proporcionalmente a la fracción de red $\frac{y(t)}{P}$. Ofrece una reproducción animada con controles de velocidad física y scrubbing temporal interactivo.
3.  **Consola en Vivo de Indicadores Matemáticos:** Computa en tiempo real propiedades críticas asociadas a la EDO autónoma, tales como el número reproductivo $R_0$, la capacidad máxima límite de carga estable $K^*$, y advierte analíticamente la convergencia del sistema basándose en la estabilidad global.
4.  **Gráfico de Trayectorias de Fase (Evolución Temporal):** Grafica de manera simultánea en formato vectorial nativo las curvas continuas arrojadas por los tres solvers numéricos (RK4, Euler y multipaso AB2) permitiendo detectar con crudeza el fenómeno de inestabilidad y divergencia del integrador explícito de Euler de primer orden ante paso de malla burdo.
5.  **Gráfica de Sensibilidad Paramétrica $y(T) \text{ vs } \beta$:** Barrido computacional que realiza un mapa en frío de cómo cambia el estado de amortiguación terminal con respecto a la magnitud de $\beta$.
6.  **Buscador Óptimo por Método del Disparo (Bisección):** Formulario avanzado que ejecuta y detalla paso a paso mediante una tabla interactiva el historial iterativo de intervalos del algoritmo de bisección, reflejando el descenso del residuo a cada iteración en milisegundos con posibilidad de inyectar instantáneamente el parámetro óptimo hallado en la simulación activa.

---

## 6. Instrucciones de Uso Académico en Laboratorio

Para corroborar de forma didáctica la teoría matemática mediante la herramienta, se sugieren los siguientes ensayos experimentales:

### Ejercicio A: Análisis de la Estabilidad Numérica (Fenómeno de Rigidez)
1.  Fije los siguientes parámetros: $\beta = 1.8 \times 10^{-5}$ y $\delta = 0.5$ (por lo tanto, $R_0 > 1$, la curva busca saturarse establemente en $K^* \approx 72,222$ personas).
2.  Ajuste el selector de **Paso de Integración Numérica (h)** en su nivel más burdo $h = 0.2$. Active la casilla de visualización del gráfico para **Euler Explícito** y **Runge-Kutta 4 (RK4)**.
3.  Observe en el gráfico de trayectorias cómo la curva azul de Euler oscila abruptamente de manera desmedida y se aparta de la referencia RK4 (acumulación sistemática del error local).
4.  Empiece a disminuir gradualmente el paso de discretización $h$ a $0.05$ primero, y finalmente a $0.01$. Verifique cómo las curvas convergen armónicamente resolviendo la inestabilidad matemática inherente de los solvers explícitos no lineales de bajo orden.

### Ejercicio B: Comprobación de Extinción de Información ($R_0 \le 1$)
1.  Fije la tasa de contagio en un valor sumamente bajo utilizando el control deslizante: $\beta = 4.0 \times 10^{-6}$.
2.  Establezca el coeficiente de olvido en un nivel alto, como por ejemplo $\delta = 0.6$.
3.  Examine cómo en el monitor interactivo de indicadores, el sistema advierte que el **Número Reproductivo Básico es $R_0 = 0.67 < 1$** y destaca la leyenda roja exponiendo que **la información se extinguirá**. 
4.  Visualice las tres curvas en el gráfico y constate cómo todas descienden suave pero asintóticamente a cero en los primeros 10 días, extinguiendo el rumor por completo de la red de nodos.

### Ejercicio C: Ejercicio Práctico de Ingeniería Química, Biológica o de Datos con Método del Disparo
1.  Usted es un analista de marketing y desea que en el **Día 15 ($t_{meta} = 15$)** la red tenga exactamente **50,000 personas informadas activas ($y_{meta} = 50,000$)**.
2.  Vaya a la consola del del **Buscador de β por Método del Disparo** e ingrese $t_{meta} = 15 y_{meta} = 50000$. Elija una precisión fina de $\pm 10$ personas.
3.  Haga clic en **Calibrar β por Bisección**.
4.  Siga con atención sobre la tabla interactiva el proceso por el cual el intervalo de búsqueda se encoge dinámicamente reduciendo secuencialmente el residuo con el teorema del valor intermedio. 
5.  Una vez exitosa la convergencia, pulse el botón verde de **Aplicar β Hallado** y observe cómo la gráfica temporal de la pantalla se expande cruzando exactamente de manera magistral el marcador virtual en el plano $(15, 50\text{k})$.
