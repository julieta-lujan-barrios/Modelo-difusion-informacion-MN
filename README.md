# Modelo de Difusión de Información con EDO

Motor de simulación numérica para modelar la propagación de información en redes sociales mediante una **Ecuación Diferencial Ordinaria no lineal** (logística con término de olvido), inspirada en modelos epidemiológicos tipo SIS.

Proyecto Final — Cátedra de **Métodos Numéricos**, Ingeniería en Inteligencia Artificial, UNSTA (2026).

🔗 **[Simulador interactivo desplegado](https://simulador-de-difusi-n-de-informaci-n-174273423336.us-west1.run.app)**

📓 **[Notebook completo (Colab)](notebook/FINAL_M_N_Equipo_3-2.ipynb)**

---

## 📖 Descripción

Se modela la evolución temporal de agentes informados $y(t)$ dentro de una red cerrada de tamaño $P$, donde la información se propaga por contacto social pero también se **olvida** con el tiempo:

$$\frac{dy}{dt} = \beta \cdot y(t)\cdot[P - y(t)] - \delta \cdot y(t), \qquad y(0) = y_0$$

| Parámetro | Valor base | Descripción |
|---|---|---|
| $P$ | 100.000 | Población total de la red |
| $\beta$ | $3.0\times10^{-6}$ | Tasa de contagio social (parámetro libre) |
| $\delta$ | $0.1$ | Tasa de olvido/desinterés (fijo) |
| $y_0$ | 100 | Agentes informados en $t=0$ |
| $T$ | 15 días | Horizonte de simulación |

El sistema tiene dos equilibrios: el trivial $\bar y_1 = 0$ (inestable) y el de saturación $\bar y_2 = K^* = P - \delta/\beta$ (estable), condicionados por el número reproductivo básico $R_0 = \beta P/\delta$.

## Fases del proyecto

**Fase 1 — Implementación de solvers.** Euler explícito, Runge-Kutta de 4° orden (RK4) y Adams-Bashforth de 2 pasos (AB2), con arranque del método multipaso vía un paso de RK4.

**Fase 2 — Estudio de convergencia.** Refinamiento de malla ($h = 0.2,\ 0.1,\ 0.05,\ 0.025$), error global máximo contra una solución de referencia y estimación del orden empírico en escala log-log. Resultado: $p_{RK4}\approx 4.03$, $p_{AB2}\approx 1.98$.

**Fase 3 — Método del Disparo.** Calibración inversa del parámetro $\beta$ mediante bisección acoplada a RK4, para que $y(10) = 0.85\,P$. Convergencia en 19 iteraciones, con $\beta^{*}\approx 1.0542\times10^{-5}$.

## Tecnologías

- **Python 3** — NumPy, Matplotlib (solvers, estudio de convergencia y bisección)
- **React + TypeScript (Vite) + Tailwind CSS** — simulador interactivo web


## Integrantes:

- Barrios, Julieta Luján
- Sappia, Lucio Agustín
- Tartarini, Luciano Aquiles

**Docente:** Cátedra de Métodos Numéricos — UNSTA
**Fecha de presentación:** 23/06/2026

## Referencias:

- Turner, P., Arildsen, T., & Kavanagh, K. (2018). *Applied Scientific Computing with Python*. Springer.
- Kermack, W. O., & McKendrick, A. G. (1927). *A contribution to the mathematical theory of epidemics*. Proc. Royal Society A.
- Murray, J. D. (2002). *Mathematical Biology: I. An Introduction* (3rd ed.). Springer.
