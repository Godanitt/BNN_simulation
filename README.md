# BBN mínimo con pynucastro

Conjunto de notebooks básicos para atacar el ejercicio 1 de *Nuclear Structure and Astrophysics – Assignments 2026*:

- snapshots BBN: t≈50 s, 200 s, 1000 s
- T9 = 1.7, 1.0, 0.4
- rho = 2e-4, 2e-5, 2e-6 g/cm3
- condición inicial n/p = 1/7

Orden recomendado:

1. `00_setup_y_condiciones.ipynb`
2. `01_red_bbn_pynucastro.ipynb`
3. `02_integracion_snapshots.ipynb`
4. `03_plots_y_tablas_minimas.ipynb`

Notas:
- Esto es un modelo mínimo de red nuclear bajo una historia T(t), rho(t) prescrita.
- No sustituye a un código BBN cosmológico completo como AlterBBN.
- Las figuras se guardan en PDF.
