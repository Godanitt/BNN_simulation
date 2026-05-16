# BBN simulation with pynucastro + AlterBBN benchmark

Paquete de notebooks para simular una red BBN ligera con `pynucastro`, arrancando en el primer snapshot del enunciado (`t0 = 50 s`), y comparar los resultados finales con la salida estándar de AlterBBN (`./stand_cosmo.x 5`) y con valores observacionales de referencia.

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
4. `03_graficas_y_comparacion_observacional.ipynb`
5. `04_tablas_heatmaps_y_gifs.ipynb`
6. `05_alterbbn_benchmark.ipynb`

También puedes usar:

```bash
cd main
bash run_all_notebooks.sh
```

El notebook 05 no ejecuta AlterBBN; toma los valores que salieron de `./stand_cosmo.x 5` y los guarda en CSV para que la memoria pueda leerlos de forma reproducible.

## Salidas principales

### pynucastro

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

El notebook 01 intenta construir una red algo más fina incluyendo `Li6` si la instalación local de `pynucastro` lo permite, y exporta un diagnóstico de una red con reacciones inversas. Para mantener una integración estable y transparente, la red principal se integra sin forzar las inversas. AlterBBN se introduce como benchmark externo de precisión.
