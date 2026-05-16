# Memoria Typst compatible con `BNN_simulation`

La memoria vive dentro de:

```text
BNN_simulation/
├── main/       # notebooks, scripts, data/, figures/, tables/, gifs/
└── memoria/    # main.typ y references.bib
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

Ejecuta de nuevo los notebooks dentro de `main/` y recompila la memoria. El orden recomendado es:

```text
00 -> 01 -> 02 -> 06 -> 03 -> 04 -> 05
```

El notebook `06_pynucastro_detailed_balance_diagnostic.ipynb` genera una prueba opcional `pynucastro_db`. El último notebook, `05_alterbbn_benchmark.ipynb`, debe ejecutarse al final porque sobrescribe las tablas y figuras comparativas con la versión que incluye AlterBBN y, si existe, también `pynucastro_db`.
