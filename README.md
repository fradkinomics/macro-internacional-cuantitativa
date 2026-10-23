# Macroeconomía Internacional Cuantitativa

Problem sets de la materia Macroeconomía Internacional Cuantitativa
(Maestría en Economía, Universidad de San Andrés, 2026), dictada por Francisco Roldán.
Los modelos se resuelven y simulan en [Julia](https://julialang.org); los informes están escritos en LaTeX.

| Problem Set | Tema | Código | Informe |
|---|---|---|---|
| [PS1](PS1/) | Modelos de búsqueda (McCall, 1970) | [`ps1.jl`](PS1/ps1.jl) | [PDF](PS1/informe/MacroInternacional_ps1_Fradkin.pdf) |
| [PS2](PS2/) | Problema de la torta | [`ps2.jl`](PS2/ps2.jl) | [PDF](PS2/informe/MacroInternacional_ps2_Fradkin.pdf) |

## PS1 — Modelos de búsqueda

Se resuelve por iteración de la función de valor el modelo de búsqueda laboral de McCall (1970) y se estudia:

1. **Estáticas comparadas** del salario de reserva respecto del seguro de desempleo *b* y del factor de descuento *β*, y el efecto de la resolución de la grilla.
2. **Simulaciones** del tiempo de búsqueda *T* (10.000 historias) y su comparación con la distribución geométrica teórica.
3. **Duración esperada del desempleo** en función de la paciencia.

`mccall.jl` es el código base provisto en el curso. `ps1.jl` contiene la resolución del problem set.

## PS2 — Problema de la torta

Se resuelve el problema de consumo y ahorro con utilidad CRRA por iteración de la función de valor, con y sin interpolación, y se compara con la forma cerrada:

1. **Funciones de consumo**: fracción de la torta consumida y ahorrada con grillas de 20 y 1000 puntos, y comparación de los dos métodos de resolución.
2. **Para pensar**: la regla de decisión con la función de valor exacta (*guess and verify*) y la reformulación del problema con la variable de estado en logaritmos, con y sin extrapolación.
3. **Simulador**: trayectorias de consumo y torta, y acumulación del error numérico en el tiempo.

`cakeeating.jl` e `itpcake.jl` son los códigos base provistos en el curso. `ps2.jl` contiene la resolución del problem set.

## Cómo reproducir los resultados

Requiere [Julia](https://julialang.org/downloads/) 1.10 o superior.

```bash
git clone https://github.com/fradkinomics/macro-internacional-cuantitativa.git
cd macro-internacional-cuantitativa
julia PS1/ps1.jl
julia PS2/ps2.jl
```

La primera ejecución instala automáticamente las versiones exactas de los paquetes registradas en el `Project.toml` y el `Manifest.toml` de cada carpeta. Los gráficos se guardan en `PS1/figures/` y `PS2/figures/`.

## Autor

Martín A. Fradkin