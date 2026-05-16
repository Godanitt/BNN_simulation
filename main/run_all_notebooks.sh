#!/usr/bin/env bash
set -euo pipefail

jupyter nbconvert --to notebook --execute 00_config_fuentes_y_termodinamica.ipynb --inplace
jupyter nbconvert --to notebook --execute 01_red_pynucastro_y_graficas_red.ipynb --inplace
jupyter nbconvert --to notebook --execute 02_integracion_modelo_t0_50.ipynb --inplace
jupyter nbconvert --to notebook --execute 03_graficas_y_comparacion_observacional.ipynb --inplace
jupyter nbconvert --to notebook --execute 04_tablas_heatmaps_y_gifs.ipynb --inplace
jupyter nbconvert --to notebook --execute 05_alterbbn_benchmark.ipynb --inplace
