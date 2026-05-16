#!/usr/bin/env bash
set -euo pipefail

# Run the BBN workflow in dependency order.
# 06 is optional/diagnostic: it builds a second pynucastro network with
# reverse rates by detailed balance if the local pynucastro installation allows it.
# 05 is run after 06 so that the final benchmark automatically includes
# pynucastro_db whenever the diagnostic integration succeeds.

NOTEBOOKS=(
  "00_config_fuentes_y_termodinamica.ipynb"
  "01_red_pynucastro_y_graficas_red.ipynb"
  "02_integracion_modelo_t0_50.ipynb"
  "06_pynucastro_detailed_balance_diagnostic.ipynb"
  "03_graficas_y_comparacion_observacional.ipynb"
  "04_tablas_heatmaps_y_gifs.ipynb"
  "05_alterbbn_benchmark.ipynb"
)

for nb in "${NOTEBOOKS[@]}"; do
  echo "[run_all] Running ${nb}"
  jupyter nbconvert --to notebook --execute --inplace "${nb}"
done

echo "[run_all] Done. Recompile memoria/main.typ if needed."
