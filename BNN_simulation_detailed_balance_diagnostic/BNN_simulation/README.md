# BBN simulation with pynucastro + detailed-balance diagnostic + AlterBBN benchmark

Paquete de notebooks para simular una red BBN ligera con `pynucastro`, arrancando en el primer snapshot del enunciado (`t0 = 50 s`), diagnosticar el problema de D/H mediante una red alternativa con tasas inversas (`pynucastro_db`), y comparar los resultados finales con la salida estándar de AlterBBN (`./stand_cosmo.x 5`) y con valores observacionales de referencia.

## Orden de ejecución

Desde la carpeta `main/`:

```bash
pip install -r requirements.txt
jupyter notebook
```

Ejecuta en este orden:

1. `00_config_fuentes_y_termodinamica.ipynb`
2. `01_red_pynucastro_y_graficas_red.ipynb`
3. `02_integracion_modelo_t0_50.ipynb`
4. `06_pynucastro_detailed_balance_diagnostic.ipynb`
5. `03_graficas_y_comparacion_observacional.ipynb`
6. `04_tablas_heatmaps_y_gifs.ipynb`
7. `05_alterbbn_benchmark.ipynb`

También puedes usar:

```bash
cd main
bash run_all_notebooks.sh
```

El notebook 06 es diagnóstico: intenta generar `bbn_network_db.py` con tasas inversas derivadas por balance detallado. Si la integración de `pynucastro_db` funciona, el notebook 05 la detecta automáticamente y la añade a las tablas y figuras finales. Si no funciona, deja un CSV de diagnóstico y la memoria sigue compilando con la comparación `pynucastro` vs AlterBBN vs observación.

El notebook 05 no ejecuta AlterBBN; toma los valores que salieron de `./stand_cosmo.x 5` y los guarda en CSV para que la memoria pueda leerlos de forma reproducible.

## Salidas principales

### pynucastro base

- `data/bbn_evolution_t0_50s.csv`
- `data/bbn_snapshots_all_models.csv`
- `data/bbn_final_compact.csv`
- `tables/table_snapshots_models.csv`
- `tables/table_final_models.csv`
- `figures/fig_all_abundances_t0_50s.pdf`
- `figures/fig_network_t0_50_snapshot_50s.pdf`
- `figures/fig_network_t0_50_snapshot_200s.pdf`
- `figures/fig_network_t0_50_snapshot_1000s.pdf`
- `gifs/gif_abundances_t0_50s.gif`

### Diagnóstico de balance detallado

- `main/06_make_detailed_balance_network.py`
- `main/06_pynucastro_detailed_balance_diagnostic.ipynb`
- `data/pynucastro_db_status.csv`
- `data/bbn_network_reverse_diagnostic.csv`
- `data/bbn_network_db_rates.csv` *(si se puede construir la red)*
- `data/bbn_network_db_derived_attempts.csv` *(si se puede construir la red)*
- `data/bbn_evolution_pynucastro_db.csv` *(si se integra con éxito)*
- `data/bbn_final_abundances_pynucastro_db.csv` *(si se integra con éxito)*

### AlterBBN y comparación final

- `data/alterbbn_standard_results.csv`
- `data/alterbbn_correlation_matrix.csv`
- `data/bbn_final_benchmark_combined.csv`
- `data/comparison_pynucastro_alterbbn_observed.csv`
- `data/comparison_models_vs_observed.csv`
- `tables/table_alterbbn_standard_results.csv`
- `tables/table_bbn_final_benchmark_combined.csv`
- `tables/table_comparison_models_vs_observed.csv`
- `figures/fig_final_he4_bar_models_vs_observed.pdf`
- `figures/fig_final_ratios_bar_models_vs_observed.pdf`
- `figures/fig_lithium7_problem_comparison.pdf`
- `figures/fig_combined_benchmark_pynucastro_alterbbn_observed.pdf`

## Nota física

La versión actual usa como modelo didáctico principal una red `pynucastro` con arranque en `t0 = 50 s`, `n/p = 1/7` y corte en `t = 1000 s`. El arranque a `10 s` se ha retirado de la memoria y de las gráficas principales.

La comparación con AlterBBN muestra que `pynucastro` reproduce bien la escala de `Yp`, pero quema demasiado deuterio. El notebook 06 prueba la hipótesis más natural: que la red base no representa suficientemente el equilibrio directo-inverso del cuello de botella del deuterio, especialmente el canal `d -> n + p`. Si `pynucastro_db` sube D/H hacia `1e-5`, el diagnóstico apunta a las inversas; si no, la limitación es más amplia y viene de la trayectoria cosmológica simplificada y de la física BBN que AlterBBN sí resuelve.
