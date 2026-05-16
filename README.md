# BBN simulation with pynucastro — notebooks v2

Paquete de notebooks para simular una red BBN ligera con `pynucastro`, comparar dos arranques (`t0=10 s` y `t0=50 s`) y generar CSVs/PDFs/GIFs para la memoria.

## Orden de ejecución

Desde la carpeta `main/`:

```bash
pip install -r requirements.txt
jupyter notebook
```

Ejecuta en este orden:

1. `00_config_fuentes_y_termodinamica.ipynb`
2. `01_red_pynucastro_y_graficas_red.ipynb`
3. `02_integracion_modelos_t0_10_y_50.ipynb`
4. `03_graficas_y_comparacion_observacional.ipynb`
5. `04_tablas_heatmaps_y_gifs.ipynb`

## Salidas principales

- `data/bbn_evolution_t0_10s.csv`
- `data/bbn_evolution_t0_50s.csv`
- `data/bbn_snapshots_all_models.csv`
- `data/bbn_final_compact.csv`
- `data/comparison_models_vs_observed.csv`
- `figures/fig_all_abundances_t0_10s.pdf`
- `figures/fig_all_abundances_t0_50s.pdf`
- `figures/fig_key_ratios_evolution_vs_observed.pdf`
- `figures/fig_he4_evolution_vs_observed.pdf`
- `figures/fig_lithium7_problem_comparison.pdf`
- `figures/fig_network_t0_50_snapshot_50s.pdf`
- `figures/fig_network_t0_50_snapshot_200s.pdf`
- `figures/fig_network_t0_50_snapshot_1000s.pdf`
- `gifs/gif_abundances_t0_10s.gif`
- `gifs/gif_abundances_t0_50s.gif`

## Nota física

Esto es una red BBN reducida. No sustituye a AlterBBN ni a un cálculo completo con neutrinos, pares e±, correcciones radiativas, etc. Está pensada para el ejercicio: entender la dependencia de abundancias ligeras con una historia termodinámica impuesta y comparar los observables principales.


## Cambios v3

- La trayectoria termodinámica ya no es una interpolación log-lineal simple. Ahora usa una referencia de Universo dominado por radiación, $T_9\propto t^{-1/2}$, y una densidad de referencia $ho_b=m_u\eta n_\gamma(T)$ con $\eta_{10}=6.10$.
- Sobre esa referencia se aplica una corrección PCHIP suave en escala logarítmica que fuerza exactamente los tres snapshots del enunciado.
- Las gráficas de red usan flechas extraídas de las reacciones reales de `pynucastro`, con `FancyArrowPatch`, `shrinkA/shrinkB` y etiquetas con contorno para mejorar la legibilidad.
