// =============================================================
// Memoria Typst: BBN con pynucastro y benchmark AlterBBN
// =============================================================

#set document(
  title: "Big Bang Nucleosynthesis with pynucastro, detailed-balance diagnostic and AlterBBN",
  author: "Daniel Vázquez",
)
#set page(paper: "a4", margin: (x: 2.25cm, y: 2.30cm))
#set text(size: 10.5pt, lang: "es")
#set par(justify: true, leading: 0.65em)
#set heading(numbering: "1.1")
#set figure(numbering: "1")
#set table(stroke: 0.45pt + luma(180), inset: 4pt)

// --------------------------
// Isótopos cómodos
// --------------------------
#let H1 = $""^1 "H"$
#let H2 = $""^2 "H"$
#let H3 = $""^3 "H"$
#let He3 = $""^3 "He"$
#let He4 = $""^4 "He"$
#let Li6 = $""^6 "Li"$
#let Li7 = $""^7 "Li"$
#let Be7 = $""^7 "Be"$

// --------------------------
// Herramientas internas CSV
// --------------------------
#let load-csv(path) = csv(path, row-type: dictionary)

#let first-n(rows, n) = {
  if rows.len() <= n { rows } else { rows.slice(0, n) }
}

#let cell(row, key) = {
  if row.keys().contains(key) {
    row.at(key)
  } else {
    text(fill: red)[missing: #key]
  }
}

#let csv-table(rows, columns, headers: none, small: false) = {
  let hs = if headers == none { columns } else { headers }
  let tab = table(
    columns: columns.len(),
    align: horizon,
    stroke: 0.45pt + luma(180),
    inset: 3.5pt,
    table.header(..hs.map(h => table.cell(fill: luma(235), strong(h)))),
    ..rows.map(row => {
      columns.map(col => table.cell(cell(row, col)))
    }).flatten()
  )
  if small { text(size: 7.2pt, tab) } else { tab }
}

#let note-box(body) = block(
  fill: luma(246),
  stroke: 0.45pt + luma(180),
  radius: 3pt,
  inset: 8pt,
)[#body]

// --------------------------
// CSVs generados por los notebooks
// --------------------------
#let assignment-snapshots = load-csv("../main/data/bbn_assignment_snapshots.csv")
#let thermo-check = load-csv("../main/tables/table_thermo_snapshot_check.csv")
#let nuclei = load-csv("../main/data/bbn_network_nuclei.csv")
#let rates = load-csv("../main/data/bbn_network_rates.csv")
#let obs = load-csv("../main/data/observational_abundances.csv")
#let snap-table = load-csv("../main/tables/table_snapshots_models.csv")
#let final-table = load-csv("../main/tables/table_final_models.csv")
#let alter-table = load-csv("../main/tables/table_alterbbn_standard_results.csv")
#let benchmark-table = load-csv("../main/tables/table_bbn_final_benchmark_combined.csv")
#let comparison-table = load-csv("../main/tables/table_comparison_models_vs_observed.csv")
#let reverse-diagnostic = load-csv("../main/data/bbn_network_reverse_diagnostic.csv")
#let db-status = load-csv("../main/data/pynucastro_db_status.csv")

#align(center)[
  #text(size: 18pt, weight: "bold")[Big Bang Nucleosynthesis with `pynucastro`]

  #v(0.25em)
  #text(size: 13pt)[detailed-balance diagnostic and comparison with AlterBBN]

  #v(0.35em)
  #text(size: 13pt)[Trabajo 1 --- Estructura Nuclear y Astrofísica]

  #v(0.75em)
  Daniel Vázquez

  #v(0.35em)
  #datetime.today().display("[day]/[month]/[year]")
]

#v(0.8em)
#outline(title: [Índice])

= Introducción y motivación

La nucleosíntesis primordial, o _Big Bang Nucleosynthesis_ (BBN), es una de las pruebas más limpias de la conexión entre cosmología, física nuclear y física de partículas. Durante los primeros minutos del Universo la temperatura y la densidad permitieron reacciones termonucleares entre protones, neutrones y núcleos ligeros. La expansión fue demasiado rápida para producir una tabla periódica completa; el resultado principal fue #H1, #H2, #He3, #He4 y trazas de masa 7, especialmente #Li7 y #Be7 @Peebles1966 @Wagoner1967.

El observable más robusto es la fracción másica de helio primordial. Si al comenzar la nucleosíntesis efectiva se adopta

$ N_n / N_p approx 1 / 7, $

entonces casi todos los neutrones supervivientes acaban ligados en partículas alfa. Dos neutrones requieren dos protones para formar #He4, y los protones restantes quedan mayoritariamente como hidrógeno. La estimación elemental da

$ X_alpha approx 4 / (4 + 12) = 0.25, $

muy cerca de los valores inferidos observacionalmente para $Y_p$ @Aver2021 @Steigman2010. Por eso #He4 es un control físico directo de la razón neutrón-protón inicial.

El deuterio cumple otro papel. La reacción

$ p + n -> #H2 + gamma $

es el primer paso para construir núcleos más pesados, pero el deuterio es frágil: mientras la cola de fotones energéticos pueda fotodesintegrarlo, la red queda bloqueada. Este retraso se conoce como _deuterium bottleneck_. Cuando el Universo se enfría lo suficiente, el deuterio sobrevive y se procesa rápidamente hacia #He3, #H3 y #He4. La abundancia residual D/H es por tanto muy sensible a la densidad bariónica @Cooke2018 @Steigman2010.

En este trabajo se hacen dos cosas distintas. Primero, se construye una red ligera y transparente con `pynucastro`, usando las condiciones termodinámicas del enunciado. Segundo, se añade una comparación externa con AlterBBN, ejecutado en modo cosmológico estándar mediante `stand_cosmo.x 5`. La comparación no pretende ajustar parámetros, sino separar claramente qué reproduce el modelo didáctico y qué predice un código BBN especializado @AlterBBNCode.

La versión final de la memoria usa un único modelo `pynucastro`: integración desde $t_0=50$ s hasta $t=1000$ s, con la razón inicial $N_n/N_p=1/7$ fijada por el enunciado. Se elimina el arranque a $10$ s porque introduce una condición inicial adicional que no está pedida y complica la lectura de las figuras.

= Condiciones del ejercicio y observables

El enunciado fija tres snapshots de nucleosíntesis primordial: $t approx 50$ s con $T_9=1.7$ y $rho=2 dot 10^(-4)$ g cm$""^(-3)$; $t approx 200$ s con $T_9=1.0$ y $rho=2 dot 10^(-5)$ g cm$""^(-3)$; y $t approx 1000$ s con $T_9=0.4$ y $rho=2 dot 10^(-6)$ g cm$""^(-3)$ @USCAssignments2026. Estos puntos se leen desde `../main/data/bbn_assignment_snapshots.csv` y se resumen en la tabla @tab:assignment-snapshots.

#figure(
  csv-table(
    assignment-snapshots,
    ("label", "t_s", "T9", "rho_g_cm3", "T_K"),
    headers: ([snapshot], [$t$ [s]], [$T_9$], [$rho$ [g cm$""^(-3)$]], [$T$ [K]]),
  ),
  caption: [Snapshots impuestos por el enunciado y usados como anclaje de la trayectoria termodinámica.],
) <tab:assignment-snapshots>

Para comparar con observaciones se usan D/H, #He3/H, $X(#He4)$ y el canal efectivo de masa 7. En la red `pynucastro` este canal se estima como

$ (#Li7 + #Be7) / "H", $

porque #Be7 producido durante BBN termina contribuyendo a #Li7 por captura electrónica. En AlterBBN se usa directamente la columna `Li7/H` de la salida estándar, que representa la abundancia final de litio de masa 7 para la comparación con el plateau observado. Los valores observacionales de referencia se cargan desde `../main/data/observational_abundances.csv`.

#figure(
  csv-table(
    obs,
    ("observable", "value", "sigma", "model_column", "source_key"),
    headers: ([observable], [valor], [$sigma$], [columna `pynucastro`], [fuente]),
    small: true,
  ),
  caption: [Valores observacionales empleados en las comparaciones. D/H procede de absorción en cuásares pobres en metales @Cooke2018, $Y_p$ de regiones H II pobres en metales @Aver2021, #He3/H se usa como control de escala @Steigman2010 y #Li7/H representa el plateau de Spite discutido en el problema cosmológico del litio @Fields2011.],
) <tab:observational-input>

= Modelo termodinámico

La trayectoria termodinámica se fuerza a pasar por los tres snapshots del enunciado. Para evitar una interpolación con quiebros artificiales se parte de una referencia de Universo dominado por radiación,

$ T(t) prop t^(-1/2), $

y de una densidad bariónica de referencia

$ rho_b(T) = m_u eta n_gamma(T), $

con

$ n_gamma(T) = (2 zeta(3)) / pi^2 (k_B T / (h c))^3. $

Sobre esa referencia se aplica una corrección suave en escala logarítmica mediante `PchipInterpolator`, de modo que la curva pase exactamente por los tres puntos pedidos. La tabla @tab:thermo-check comprueba que la trayectoria reproduce los snapshots.

#figure(
  csv-table(
    thermo-check,
    ("label", "t_s", "T9_assignment", "T9_model", "rho_assignment_g_cm3", "rho_model_g_cm3", "eta10_effective"),
    headers: ([snapshot], [$t$ [s]], [$T_9$ enunciado], [$T_9$ modelo], [$rho$ enunciado], [$rho$ modelo], [$eta_10$ efectiva]),
    small: true,
  ),
  caption: [Control de la trayectoria termodinámica. La $eta_10$ efectiva no se interpreta como una medida cosmológica, sino como diagnóstico de la densidad impuesta por el enunciado.],
) <tab:thermo-check>

#figure(
  image("../main/figures/fig_thermo_T9_forced_radiation.pdf", width: 85%),
  caption: [Historia de temperatura usada por la red. La curva pasa por los puntos del enunciado y mantiene como referencia la ley radiación-dominada $T_9 prop t^(-1/2)$.],
) <fig:thermo-t9>

#figure(
  image("../main/figures/fig_thermo_rho_forced_eta.pdf", width: 85%),
  caption: [Densidad bariónica impuesta. La referencia se calcula mediante $rho_b=m_u eta n_gamma(T)$ usando una $eta_10$ cosmológica de orden $6$.],
) <fig:thermo-rho>

#figure(
  image("../main/figures/fig_eta10_effective_forced.pdf", width: 85%),
  caption: [Razón barión-fotón efectiva asociada a la trayectoria forzada.],
) <fig:eta-effective>

= Red nuclear generada con `pynucastro`

La red principal incluye neutrones, protones, deuterio, tritio, #He3, #He4, #Li7 y #Be7. Esta es la red base usada para la integración `pynucastro`. Además se añade ahora un diagnóstico independiente, en el notebook `06_pynucastro_detailed_balance_diagnostic.ipynb`, que intenta construir una segunda red llamada `pynucastro_db`. Esta red añade tasas inversas mediante balance detallado cuando `pynucastro` puede generarlas y, como respaldo, revisa las tasas inversas disponibles en ReacLib. El objetivo de esa segunda red no es maquillar el resultado, sino probar si la discrepancia de D/H se debe principalmente a un tratamiento incompleto del cuello de botella del deuterio.

#figure(
  csv-table(
    nuclei,
    ("index", "name", "A", "Z", "N", "A_minus_2Z"),
    headers: ([índice], [núcleo], [$A$], [$Z$], [$N$], [$A-2Z$]),
  ),
  caption: [Núcleos incluidos en la red `pynucastro` generada.],
) <tab:nuclei>

La tabla @tab:rates muestra una selección de reacciones exportadas por la red. El fichero completo queda guardado en `../main/data/bbn_network_rates.csv`.

#figure(
  csv-table(
    first-n(rates, 18),
    ("index", "rate", "reactants", "products", "Q_MeV"),
    headers: ([índice], [rate], [reactivos], [productos], [$Q$ [MeV]]),
    small: true,
  ),
  caption: [Primeras reacciones de la red `pynucastro`.],
) <tab:rates>

#note-box[
  La red `pynucastro` no pretende sustituir un código BBN de precisión. Su valor es que permite ver de forma transparente qué canales controlan #He4, D/H, #He3/H y masa 7 bajo las condiciones impuestas. Precisamente por eso se conserva una red base simple y se añade una prueba `pynucastro_db`: si las tasas inversas elevan D/H hacia $10^(-5)$, el problema principal era el equilibrio directo-inverso; si no lo hacen, la discrepancia apunta a limitaciones más globales de la trayectoria cosmológica simplificada.
]

== Prueba de balance detallado y canales inversos

La comparación inicial con AlterBBN mostró una diferencia muy grande en D/H: la red base daba una abundancia de deuterio de orden $10^(-9)$, mientras que AlterBBN y la referencia observacional están en la escala $10^(-5)$. Esa diferencia no puede explicarse por una conversión incorrecta entre fracciones másicas y cocientes de número, porque ese error sería como mucho de un factor dos para el deuterio. Tampoco parece venir de la razón inicial $N_n/N_p$, porque $X(#He4)$ sale cerca de $0.25$, justo lo esperado si casi todos los neutrones acaban en partículas alfa.

La pista física es que #He3/H no se aleja tanto como D/H. Eso sugiere que el deuterio sí se forma, pero fluye demasiado pronto hacia #He3, #H3 y #He4. La hipótesis natural es entonces la ausencia, o insuficiente peso, de reacciones inversas de fotodesintegración. En BBN real, el canal

$ p + n <-> #H2 + gamma $

no funciona solo hacia la derecha: la enorme población de fotones retrasa la supervivencia neta de #H2 hasta que baja la temperatura. Si la red no representa bien ese equilibrio, el deuterio se quema antes de tiempo y queda una D/H artificialmente baja.

El notebook `06` implementa una prueba explícita. Construye `bbn_network_db.py`, exporta todas sus tasas a `../main/data/bbn_network_db_rates.csv`, guarda los intentos de derivación en `../main/data/bbn_network_db_derived_attempts.csv` y resume los canales inversos clave en la tabla @tab:reverse-diagnostic. Si la integración de `pynucastro_db` termina correctamente, sus resultados se añaden automáticamente a las tablas finales y al benchmark generado por el notebook `05`.

#figure(
  csv-table(
    db-status,
    ("status", "message"),
    headers: ([estado], [mensaje]),
    small: true,
  ),
  caption: [Estado de la prueba `pynucastro_db`. En un entorno sin `pynucastro` puede aparecer como `not_run`; al ejecutarlo localmente se actualiza con el resultado de la integración y los valores finales.],
) <tab:db-status>

#figure(
  csv-table(
    reverse-diagnostic,
    ("channel", "reactants_expected", "products_expected", "found", "rate_class", "derived_from_inverse", "diagnostic"),
    headers: ([canal], [reactivos], [productos], [encontrado], [clase], [derivada], [diagnóstico]),
    small: true,
  ),
  caption: [Diagnóstico de canales inversos relevantes para el cuello de botella del deuterio y el flujo de masa 7. El canal crítico es `d_to_n_p`, equivalente a #H2$(gamma,n)$p.],
) <tab:reverse-diagnostic>

== Grafos de red en los snapshots

Para visualizar la evolución de la composición se genera un grafo de red en los tres instantes del enunciado. Los nodos se colocan en el plano $(Z, A-2Z)$, el color codifica $log_10 X_i$ y las flechas muestran conexiones nucleares extraídas de las tasas reales de `pynucastro`.

#figure(
  image("../main/figures/fig_network_t0_50_snapshot_50s.pdf", width: 82%),
  caption: [Grafo de red evaluado en $t=50$ s.],
) <fig:network-50>

#figure(
  image("../main/figures/fig_network_t0_50_snapshot_200s.pdf", width: 82%),
  caption: [Grafo de red evaluado en $t=200$ s.],
) <fig:network-200>

#figure(
  image("../main/figures/fig_network_t0_50_snapshot_1000s.pdf", width: 82%),
  caption: [Grafo de red evaluado en $t=1000$ s.],
) <fig:network-1000>

= Integración numérica con `pynucastro`

La integración base resuelve un modelo con la condición inicial

$ Y_n(t_0)=1/8, quad Y_p(t_0)=7/8, $

que corresponde a $N_n/N_p=1/7$. Todas las demás abundancias se inicializan a cero o a un valor despreciable. La integración comienza en $t_0=50$ s y termina en $t=1000$ s, el último snapshot pedido. Las ecuaciones tienen la forma

$ d Y_i / d t = F_i(Y, rho(t), T(t)), $

con $T(t)$ y $rho(t)$ dados por la trayectoria forzada. La evolución completa se exporta a

```text
../main/data/bbn_evolution_t0_50s.csv
../main/data/bbn_evolution_all_models.csv
```

y las tablas compactas se guardan en `../main/tables/` para ser leídas directamente por esta memoria.

Si el notebook `06` se ejecuta con éxito, se añade una segunda integración `pynucastro_db` a `../main/data/bbn_evolution_all_models.csv` y a `../main/data/bbn_final_abundances_all_models.csv`. El benchmark del notebook `05` detecta automáticamente esa fila y la incorpora a las figuras finales.

= Resultados `pynucastro`

== Evolución completa

La figura @fig:all-abundances-50 muestra la evolución de las fracciones másicas $X_i=A_iY_i$. Se observa el paso desde una composición dominada por protones y neutrones hacia una mezcla final dominada por #H1 y #He4, con residuos minoritarios de D, #He3 y masa 7.

#figure(
  image("../main/figures/fig_all_abundances_t0_50s.pdf", width: 88%),
  caption: [Evolución de las fracciones másicas para el modelo `pynucastro` iniciado en $t_0=50$ s.],
) <fig:all-abundances-50>

== Abundancias en los snapshots

La tabla @tab:snapshots-models resume las magnitudes principales en los tres tiempos pedidos por el enunciado.

#figure(
  csv-table(
    snap-table,
    ("model", "t0_s", "t_s", "T9", "rho_g_cm3", "eta10", "D/H", "He3/H", "X_He4", "Li7/H", "Be7/H", "Li7_plus_Be7_over_H", "baryon_sum"),
    headers: ([modelo], [$t_0$], [$t$], [$T_9$], [$rho$], [$eta_10$], [D/H], [#He3/H], [$X(#He4)$], [#Li7/H], [#Be7/H], [(#Li7+#Be7)/H], [$sum X_i$]),
    small: true,
  ),
  caption: [Resultados del modelo en los snapshots del enunciado. La columna $sum X_i$ sirve como control de conservación bariónica numérica.],
) <tab:snapshots-models>

== Valores finales

La tabla @tab:final-models recoge el valor final a $t=1000$ s. No se extiende la integración a tiempos mucho mayores porque las condiciones termodinámicas del ejercicio solo están especificadas hasta ese instante.

#figure(
  csv-table(
    final-table,
    ("model", "t0_s", "t_s", "T9", "rho_g_cm3", "eta10", "D/H", "He3/H", "X_He4", "Li7/H", "Be7/H", "Li7_plus_Be7_over_H", "baryon_sum"),
    headers: ([modelo], [$t_0$], [$t$ final], [$T_9$], [$rho$], [$eta_10$], [D/H], [#He3/H], [$X(#He4)$], [#Li7/H], [#Be7/H], [(#Li7+#Be7)/H], [$sum X_i$]),
    small: true,
  ),
  caption: [Valores finales del modelo principal `pynucastro`.],
) <tab:final-models>

= Benchmark externo con AlterBBN

Se ejecutó AlterBBN en modo cosmológico estándar con

```text
./stand_cosmo.x 5
```

La salida central y sus incertidumbres se guardaron en `../main/data/alterbbn_standard_results.csv`. La tabla @tab:alterbbn-results reproduce los valores usados después en las comparaciones. La matriz de correlación completa queda en `../main/data/alterbbn_correlation_matrix.csv`.

#figure(
  csv-table(
    alter-table,
    ("observable", "value", "sigma", "low", "cent", "high"),
    headers: ([observable], [valor], [$sigma$], [low], [cent], [high]),
    small: true,
  ),
  caption: [Resultados de AlterBBN introducidos a partir de la salida `stand_cosmo.x 5`.],
) <tab:alterbbn-results>

La tabla @tab:benchmark-combined junta en un único fichero los resultados finales de `pynucastro`, la prueba opcional `pynucastro_db` si se ha ejecutado correctamente, y AlterBBN. Para las redes `pynucastro`, la columna #Li7/H se interpreta como la suma posterior #Li7+#Be7. Para AlterBBN, se usa la columna final #Li7/H que entrega el propio código.

#figure(
  csv-table(
    benchmark-table,
    ("model", "Yp_He4", "D/H", "He3/H", "Li7/H", "Li6/H", "Be7/H", "source"),
    headers: ([modelo], [$Y_p$], [D/H], [#He3/H], [#Li7/H], [#Li6/H], [#Be7/H], [fuente]),
    small: true,
  ),
  caption: [Resumen final combinado `pynucastro`--AlterBBN.],
) <tab:benchmark-combined>

= Comparación con observaciones

La figura @fig:key-ratios compara la evolución temporal de `pynucastro` con los valores de AlterBBN y con las bandas observacionales. La diferencia más visible aparece en D/H: el modelo didáctico procesa demasiado deuterio bajo la trayectoria impuesta, mientras que AlterBBN queda cerca del valor observado.

#figure(
  image("../main/figures/fig_key_ratios_evolution_vs_observed.pdf", width: 88%),
  caption: [Evolución de D/H, #He3/H y masa 7 para el modelo `pynucastro`, comparada con AlterBBN y con las bandas observacionales usadas en el notebook.],
) <fig:key-ratios>

El caso de #He4 se representa por separado en @fig:he4-evolution porque es una fracción másica de orden $0.25$, no un cociente pequeño.

#figure(
  image("../main/figures/fig_he4_evolution_vs_observed.pdf", width: 84%),
  caption: [Evolución de la fracción másica $X(#He4)$ y comparación con AlterBBN y con una banda observacional de $Y_p$.],
) <fig:he4-evolution>

#figure(
  image("../main/figures/fig_final_he4_bar_models_vs_observed.pdf", width: 70%),
  caption: [Comparación final de $X(#He4)$ entre `pynucastro`, AlterBBN y el valor observacional adoptado.],
) <fig:he4-bars>

La comparación final de cocientes se condensa en @fig:final-ratios-bars y en la tabla @tab:comparison-observed.

#figure(
  image("../main/figures/fig_final_ratios_bar_models_vs_observed.pdf", width: 82%),
  caption: [Comparación final de cocientes de abundancia entre `pynucastro`, AlterBBN y los valores observacionales.],
) <fig:final-ratios-bars>

#figure(
  image("../main/figures/fig_combined_benchmark_pynucastro_alterbbn_observed.pdf", width: 80%),
  caption: [Resumen visual de los cocientes pequeños D/H, #He3/H y #Li7/H.],
) <fig:combined-benchmark>

#figure(
  csv-table(
    comparison-table,
    ("observable", "kind", "value", "sigma", "source"),
    headers: ([observable], [caso], [valor], [$sigma$], [fuente]),
    small: true,
  ),
  caption: [Tabla completa de comparación entre `pynucastro`, AlterBBN y los valores observacionales.],
) <tab:comparison-observed>

= Diagnóstico del problema del litio

El canal de masa 7 se trata separadamente porque #Be7 producido durante BBN se transforma posteriormente en #Li7. En BBN estándar, el problema cosmológico del litio aparece porque las predicciones de #Li7 quedan por encima del plateau observado en estrellas pobres en metales @Fields2011. Aquí se observa la misma tensión en AlterBBN: su valor de #Li7/H es mayor que la referencia observacional adoptada, aunque mucho menor que el canal de masa 7 producido por la red `pynucastro` simplificada.

#figure(
  image("../main/figures/fig_lithium7_problem_comparison.pdf", width: 78%),
  caption: [Diagnóstico del litio: comparación de la abundancia final de masa 7 entre `pynucastro`, AlterBBN y la banda observacional usada para #Li7/H.],
) <fig:lithium-problem>

= Mapa de abundancias y animación suplementaria

La figura @fig:heatmap-50 muestra un mapa de $log_10 X_i$ en función del tiempo y de la especie. Es una forma compacta de ver qué núcleos dominan y cuándo aparecen los canales minoritarios.

#figure(
  image("../main/figures/fig_heatmap_abundances_t0_50s.pdf", width: 88%),
  caption: [Mapa de abundancias para el modelo `pynucastro`.],
) <fig:heatmap-50>

También se exporta un GIF suplementario:

```text
../main/gifs/gif_abundances_t0_50s.gif
```

El GIF no se incrusta en el PDF porque la memoria final es estática; para impresión o entrega, el mapa @fig:heatmap-50 contiene la misma información de forma compacta.

= Discusión

La simulación `pynucastro` reproduce bien la escala física básica de #He4. La razón inicial $N_n/N_p=1/7$ hace que casi todos los neutrones disponibles se incorporen a #He4, dando una fracción másica cercana a $0.25$. Este resultado coincide en escala tanto con AlterBBN como con la estimación elemental. Por tanto, el primer control físico de la red funciona: no estamos fallando en el balance bruto de neutrones.

El problema aparece al mirar D/H. La discrepancia con AlterBBN es de varios órdenes de magnitud, así que no es un detalle de red ni una cuestión de redondeo. El proceso seguido para diagnosticarlo fue:

+ comprobar que #He4 es correcto en escala, por lo que la condición inicial $N_n/N_p=1/7$ no parece ser el origen principal;
+ comprobar que el error de conversión entre fracciones másicas y cocientes de número no puede explicar una diferencia de $10^3$--$10^4$;
+ observar que #He3/H no está tan hundido como D/H, lo que indica que el deuterio se produce pero se procesa demasiado pronto;
+ identificar el cuello de botella del deuterio como sospechoso natural, porque requiere competencia entre capturas nucleares y fotodesintegraciones inversas;
+ implementar la prueba `pynucastro_db` con canales inversos derivados por balance detallado, dejando un CSV diagnóstico para ver si aparece realmente #H2$(gamma,n)$p.

Si `pynucastro_db` eleva D/H hacia la escala $10^(-5)$, la interpretación es clara: la red base estaba quemando deuterio demasiado pronto al no representar adecuadamente el equilibrio directo-inverso. Si incluso `pynucastro_db` sigue dando D/H muy baja, entonces la limitación no es solo la lista de reacciones: también entra la trayectoria termodinámica impuesta, la falta de una cosmología BBN autoconsistente, el tratamiento completo de tasas débiles y el baño de fotones, pares $e^+e^-$ y neutrinos. En ambos casos el diagnóstico es útil, porque separa una hipótesis nuclear concreta de una limitación cosmológica más global.

El canal de masa 7 es todavía más sensible. La red simplificada sobreproduce #Be7 en la trayectoria usada, mientras que AlterBBN predice la escala estándar de #Li7/H. Aun así, incluso AlterBBN queda por encima del valor observado adoptado, ilustrando el problema cosmológico del litio @Fields2011 @Vangioni2000.

= Limitaciones del cálculo

Las principales limitaciones son:

+ La trayectoria termodinámica de `pynucastro` se fuerza a los snapshots del enunciado. Aunque se usa una referencia radiación-dominada para suavizar la curva, no se resuelven las ecuaciones de Friedmann acopladas a la red nuclear.
+ La red `pynucastro` base es reducida y está diseñada para ser transparente y estable. La prueba `pynucastro_db` intenta añadir inversas por balance detallado, pero sigue sin convertir el cálculo en un código BBN completo.
+ La razón inicial $N_n/N_p=1/7$ se impone como condición inicial. En un cálculo cosmológico completo debería salir del freeze-out débil y del decaimiento beta del neutrón.
+ La comparación observacional usa valores de referencia representativos, pero no realiza una inferencia estadística de parámetros.
+ AlterBBN se usa aquí como benchmark externo a partir de una única ejecución estándar. No se hace todavía un barrido en $eta$, vida media del neutrón o tasas nucleares.

Estas limitaciones no eliminan la utilidad del ejercicio: el objetivo es construir una simulación mínima, documentada y reproducible, y después contrastarla con un código especializado para identificar qué parte de la física falta en el modelo didáctico.

= Conclusiones

Se ha construido una simulación de BBN con `pynucastro` usando una red de núcleos ligeros y una trayectoria termodinámica forzada a los tres snapshots del enunciado. La integración principal comienza en $t_0=50$ s, con $N_n/N_p=1/7$, y termina en $t=1000$ s.

La memoria incorpora ahora un benchmark con AlterBBN usando los resultados obtenidos con `stand_cosmo.x 5`. Todos los valores se guardan en CSV y se combinan en una única tabla `pynucastro`--`pynucastro_db`--AlterBBN--observación cuando la prueba de balance detallado puede ejecutarse. Esto permite presentar el trabajo de forma más honesta: `pynucastro` se usa como red transparente para entender canales nucleares, mientras que AlterBBN proporciona la escala estándar de precisión.

El resultado más robusto del modelo `pynucastro` es #He4, cercano a $X(#He4) approx 0.25$. Las discrepancias en D/H y masa 7 muestran precisamente las limitaciones esperables de una red reducida con trayectoria impuesta. La nueva prueba `pynucastro_db` permite comprobar si D/H mejora al introducir inversas por balance detallado; si no mejora, la memoria ya deja claro que el problema es más amplio que una reacción aislada. AlterBBN recupera la escala estándar de D/H, #He3/H y $Y_p$, y además reproduce la tensión conocida del litio primordial.

= Apéndice: ficheros generados

Los notebooks generan las siguientes salidas principales:

```text
../main/data/bbn_thermo_curve_forced.csv
../main/data/bbn_thermo_snapshot_check.csv
../main/data/bbn_network_nuclei.csv
../main/data/bbn_network_rates.csv
../main/data/bbn_network_reverse_diagnostic.csv
../main/data/bbn_network_db_rates.csv
../main/data/bbn_network_db_derived_attempts.csv
../main/data/pynucastro_db_status.csv
../main/data/bbn_evolution_pynucastro_db.csv
../main/data/bbn_final_abundances_pynucastro_db.csv
../main/data/bbn_evolution_t0_50s.csv
../main/data/bbn_snapshots_all_models.csv
../main/data/bbn_final_compact.csv
../main/data/alterbbn_standard_results.csv
../main/data/alterbbn_correlation_matrix.csv
../main/data/bbn_final_benchmark_combined.csv
../main/data/comparison_pynucastro_alterbbn_observed.csv
../main/data/comparison_models_vs_observed.csv
../main/tables/table_thermo_snapshot_check.csv
../main/tables/table_snapshots_models.csv
../main/tables/table_final_models.csv
../main/tables/table_alterbbn_standard_results.csv
../main/tables/table_bbn_final_benchmark_combined.csv
../main/tables/table_comparison_models_vs_observed.csv
../main/figures/*.pdf
../main/gifs/gif_abundances_t0_50s.gif
```

La memoria no copia esos ficheros: los lee directamente desde `../main/`. Basta con volver a ejecutar los notebooks y recompilar `main.typ` para actualizar tablas y figuras.

= Referencias

#bibliography("references.bib", title: none)
