# Memoria Typst compatible con `BNN_simulation`

Esta carpeta está pensada para vivir dentro de:

```text
BNN_simulation/
├── main/
│   ├── bbn_snapshots.csv
│   ├── bbn_network_nuclei.csv
│   ├── bbn_evolution.csv
│   ├── bbn_snapshot_table_clean.csv
│   ├── fig_thermo_history.pdf
│   ├── fig_abundances_evolution.pdf
│   └── fig_light_element_ratios.pdf
└── memoria/
    ├── main.typ
    ├── references.bib
    └── typst/csv_tools.typ
```

La memoria **no copia** CSVs ni PDFs. Los lee directamente desde `../main/`.

## Compilar

Desde la raíz del proyecto:

```bash
cd memoria
typst compile main.typ memoria_bbn.pdf
```

O en una línea:

```bash
typst compile memoria/main.typ memoria/memoria_bbn.pdf
```

## Actualizar resultados

Ejecuta de nuevo los notebooks dentro de `main/` y recompila la memoria. Como los datos se leen desde `../main/`, no hace falta copiar nada.
