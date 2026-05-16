// =============================================================
// Memoria Typst: BBN con pynucastro
// Ubicación esperada:
//   BNN_SIMULATION/
//   ├── main/       -> notebooks, data/, figures/, tables/, gifs/
//   └── memoria/    -> este main.typ y references.bib
// =============================================================

#set document(
  title: "Big Bang Nucleosynthesis with pynucastro",
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
#let Li7 = $""^7 "Li"$
#let Be7 = $""^7 "Be"$

// --------------------------
// Herramientas internas CSV
// --------------------------
#let load-csv(path) = csv(path, row-type: dictionary)

#let first-n(rows, n) = {
  if rows.len() <= n { rows } else { rows.slice(0, n) }
}

#let last-row(rows) = {
  if rows.len() == 0 { panic("CSV vacío.") }
  rows.at(rows.len() - 1)
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

#let key-value-table(row, entries) = {
  table(
    columns: (auto, 1fr),
    align: horizon,
    stroke: 0.45pt + luma(180),
    inset: 4pt,
    table.header(
      table.cell(fill: luma(235), strong[Magnitud]),
      table.cell(fill: luma(235), strong[Valor]),
    ),
    ..entries.map(entry => {
      (
        table.cell(entry.at(0)),
        table.cell(cell(row, entry.at(1))),
      )
    }).flatten()
  )
}

#let note-box(body) = block(
  fill: luma(246),
  stroke: 0.45pt + luma(180),
  radius: 3pt,
  inset: 8pt,
)[#body]

// --------------------------
// Lectura de CSVs generados por los notebooks
// --------------------------
#let assignment-snapshots = load-csv("../main/data/bbn_assignment_snapshots.csv")
#let thermo-check = load-csv("../main/data/bbn_thermo_snapshot_check.csv")
#let thermo-curve = load-csv("../main/data/bbn_thermo_curve_forced.csv")
#let nuclei = load-csv("../main/data/bbn_network_nuclei.csv")
#let rates = load-csv("../main/data/bbn_network_rates.csv")
#let obs = load-csv("../main/data/observational_abundances.csv")
#let obs-refs = load-csv("../main/data/references_observational_data.csv")
#let snap-table = load-csv("../main/tables/table_snapshots_models.csv")
#let final-table = load-csv("../main/tables/table_final_models.csv")
#let comparison-table = load-csv("../main/tables/table_comparison_models_vs_observed.csv")
#let final-row = last-row(final-table)

#align(center)[
  #text(size: 18pt, weight: "bold")[Big Bang Nucleosynthesis with `pynucastro`]

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

La nucleosíntesis primordial, o _Big Bang Nucleosynthesis_ (BBN), constituye una de las primeras situaciones físicas en las que la cosmología, la física nuclear y la física de partículas quedan acopladas de forma directa. Durante los primeros minutos de evolución del Universo, la temperatura y la densidad todavía eran suficientemente altas como para permitir reacciones termonucleares entre protones, neutrones y núcleos ligeros. Sin embargo, la rápida expansión cosmológica hizo que este proceso durase poco tiempo: no se trató de una nucleosíntesis prolongada como en las estrellas, sino de una ventana temporal breve, aproximadamente entre decenas y algunos cientos de segundos después del Big Bang. Los primeros cálculos modernos de abundancias primordiales ya mostraron que, en un Universo caliente y en expansión, la producción significativa queda esencialmente limitada a D, #He3, #He4 y #Li7 @Peebles1966 @Wagoner1967.

En ese intervalo se sintetizaron principalmente #H2, #He3, #He4 y pequeñas cantidades de #Li7 y #Be7. La importancia de la BBN reside en que sus productos pueden compararse con abundancias observadas en entornos astrofísicos poco procesados químicamente. En particular, la abundancia primordial de #He4 puede estimarse mediante observaciones de regiones H II pobres en metales, típicamente en galaxias enanas, extrapolando la abundancia de helio a metalicidad nula. Valores observacionales modernos sitúan la fracción másica primordial de helio alrededor de

$ X_alpha equiv X(#He4) approx 0.246 - 0.247, $

con pequeñas diferencias según el análisis espectroscópico y el tratamiento de sistemáticos @Aver2021 @Steigman2010. Esta cercanía a $0.25$ no es accidental, sino que refleja directamente la razón neutrón-protón disponible cuando comienza la nucleosíntesis efectiva.

La razón por la que aparece naturalmente una fracción de helio de ese orden es sencilla. Si en el momento en que comienzan eficientemente las reacciones nucleares tenemos aproximadamente

$ N_n / N_p approx 1 / 7, $

entonces, por cada neutrón disponible, hay siete protones. Como casi todos los neutrones que sobreviven acaban ligados en #He4, dos neutrones y dos protones forman una partícula alfa, quedando el exceso de protones como hidrógeno. En una cuenta mínima, dos neutrones van acompañados de catorce protones: dos de esos protones entran en #He4 y doce quedan libres. Así, la fracción másica esperada de helio es

$ X_alpha approx 4 / (4 + 12) = 0.25. $

Esta estimación no depende de los detalles finos de la red nuclear, sino de una idea física muy robusta: en BBN casi todos los neutrones supervivientes terminan en #He4. Por ello, #He4 es uno de los observables más limpios del modelo estándar de nucleosíntesis primordial.

El deuterio juega un papel diferente. Aunque la reacción

$ p + n -> #H2 + gamma $

es el primer paso necesario para construir núcleos más pesados, el deuterio es muy frágil frente a fotodesintegración. Mientras la cola de alta energía de la distribución de Planck contenga suficientes fotones con energía comparable o superior a la energía de ligadura del deuterio, $B_d approx 2.225 "MeV"$, la abundancia de #H2 permanece muy baja. Este retraso se conoce como _deuterium bottleneck_. Solo cuando el Universo se enfría lo suficiente, el deuterio puede sobrevivir; entonces la red nuclear avanza rápidamente hacia #He3, #H3 y finalmente #He4. Por eso la abundancia final de deuterio es especialmente sensible a la densidad bariónica: si hay más bariones por fotón, el deuterio se procesa más eficientemente hacia helio y su abundancia residual disminuye @Cooke2018 @Steigman2010.

Desde el punto de vista de la física nuclear, el problema se formula mediante una red de ecuaciones acopladas para las abundancias molares $Y_i$,

$ (d Y_i) / (d t) = F_i(Y_j, rho(t), T(t), lambda_k), $

donde las funciones $F_i$ contienen las tasas termonucleares de captura, intercambio de partículas y decaimiento beta. Las tasas dependen de promedios térmicos del tipo

$ chevron.l sigma v chevron.r = integral_0^infinity sigma(E) v(E) P(E, T) dif E, $

por lo que la evolución de abundancias queda determinada por la historia termodinámica $T(t)$, $rho(t)$, las secciones eficaces nucleares y las condiciones iniciales. En este sentido, la BBN es un banco de pruebas especialmente limpio: una red nuclear relativamente pequeña produce observables cosmológicos medibles.

Además, el cálculo está fuertemente conectado con la física de partículas. Antes de la nucleosíntesis efectiva, las interacciones débiles mantienen aproximadamente el equilibrio entre protones y neutrones mediante procesos como

$ n + nu_e <-> p + e^-, quad n + e^+ <-> p + overline(nu)_e. $

En equilibrio térmico, la razón neutrón-protón está controlada por la diferencia de masas entre el neutrón y el protón,

$ N_n / N_p approx exp(-Q / (k_B T)), quad Q = (m_n - m_p)c^2 approx 1.293 "MeV". $

Cuando el ritmo de expansión del Universo supera al ritmo de estas interacciones débiles, se produce el _freeze-out_ débil. A partir de ese momento, la razón $N_n / N_p$ deja de seguir el equilibrio térmico y solo evoluciona de forma apreciable por el decaimiento beta del neutrón. Por esta razón, observables como $X(#He4)$ son sensibles a parámetros de física de partículas como la vida media del neutrón, el número efectivo de especies relativistas y la tasa de expansión del Universo @Peebles1966 @Steigman2010.

Otro parámetro central es la razón barión-fotón,

$ eta = n_b / n_gamma. $

En BBN estándar, una vez fijadas las tasas nucleares y el contenido relativista del Universo, $eta$ controla la densidad de bariones disponible para formar núcleos. Este parámetro puede inferirse de forma independiente a partir del fondo cósmico de microondas, lo que convierte a la comparación entre BBN y CMB en una prueba cruzada del modelo cosmológico estándar @Planck2018. En particular, el deuterio actúa como un barómetro muy sensible: un valor mayor de $eta$ implica una mayor densidad bariónica a temperatura dada, favorece que las reacciones nucleares compitan antes con la expansión y reduce la abundancia final de D/H.

En este trabajo se construye una simulación numérica sencilla de la BBN con una red de núcleos ligeros generada con `pynucastro`. El objetivo no es competir con códigos especializados como AlterBBN, sino entender de forma controlada qué física entra en el problema y hasta dónde puede llegar un modelo mínimo. Para ello se estudian dos condiciones iniciales: una evolución iniciada alrededor de $t = 10$ s y otra iniciada alrededor de $t = 50$ s, ambas con una razón inicial

$ N_n / N_p = 1 / 7. $

La comparación entre ambos casos permite separar parcialmente dos efectos: por un lado, la evolución nuclear una vez fijada la composición inicial; por otro, la sensibilidad del resultado a la elección del instante inicial y, por tanto, al tiempo disponible para que actúen las reacciones.

Los resultados principales que se buscarán son las abundancias finales de #H1, #H2, #He3, #He4, #Li7 y #Be7, así como los cocientes observacionalmente relevantes

$ "D/H", quad #He3 "/" "H", quad X(#He4), quad #Li7 "/" "H", quad #Be7 "/" "H". $

Estos valores se comparan con abundancias primordiales observadas y con expectativas de BBN estándar. En particular, se comprueba si el modelo reproduce correctamente la escala de $X(#He4) approx 0.25$ y si obtiene una abundancia de deuterio razonable. También se analiza el canal de masa 7, teniendo en cuenta que #Be7 decae posteriormente a #Li7 por captura electrónica, de modo que el observable relevante no es solo #Li7 producido directamente, sino la suma efectiva

$ (#Li7 + #Be7) / "H". $

Finalmente, se presta atención al conocido problema cosmológico del litio. Mientras que la BBN estándar reproduce de forma notable las abundancias de #He4 y deuterio, la abundancia predicha de #Li7 suele quedar por encima de la observada en estrellas pobres en metales. Esta discrepancia puede deberse a una combinación de efectos astrofísicos, como depleción de litio en atmósferas estelares, incertidumbres sistemáticas observacionales, detalles de evolución química galáctica, o incluso a física no estándar @Fields2011. Además, la propia producción de LiBeB está limitada por la estructura nuclear: la ausencia de núcleos estables con $A = 5$ y $A = 8$ impide que la BBN avance eficientemente hacia elementos más pesados, dejando al litio y berilio como productos finales minoritarios y delicados @Vangioni2000.

Así, el propósito del trabajo es doble. Primero, reproducir de forma transparente la física nuclear básica de la BBN a partir de una red de reacciones ligeras y una historia termodinámica dada. Segundo, usar las abundancias finales como diagnóstico: verificar la robustez de la producción de #He4, estudiar la sensibilidad del deuterio a la densidad bariónica y observar hasta qué punto el canal de #Li7/#Be7 se desvía de las abundancias inferidas observacionalmente.

= Condiciones del ejercicio y observables

El enunciado del trabajo fija tres snapshots de nucleosíntesis primordial: $t approx 50$ s con $T_9=1.7$ y $rho=2 dot 10^(-4)$ g cm$""^(-3)$; $t approx 200$ s con $T_9=1$ y $rho=2 dot 10^(-5)$ g cm$""^(-3)$; y $t approx 1000$ s con $T_9=0.4$ y $rho=2 dot 10^(-6)$ g cm$""^(-3)$. Además, pide asumir una razón inicial neutrón-protón de $1/7$ @USCAssignments2026. Los puntos de entrada utilizados por los notebooks se leen directamente desde `../main/data/bbn_assignment_snapshots.csv` y se resumen en la tabla @tab:assignment-snapshots.

#figure(
  csv-table(
    assignment-snapshots,
    ("label", "t_s", "T9", "rho_g_cm3", "T_K"),
    headers: ([snapshot], [$t$ [s]], [$T_9$], [$rho$ [g cm$""^(-3)$]], [$T$ [K]]),
  ),
  caption: [Snapshots impuestos por el enunciado y usados como puntos de anclaje de la trayectoria termodinámica.],
) <tab:assignment-snapshots>

Para poder comparar las salidas de la red con datos observacionales, se define un conjunto de observables mínimos: D/H, #He3/H, $X(#He4)$ y el canal de masa 7, $((#Li7 + #Be7)/H)$. Los valores usados en las gráficas no se escriben a mano en la memoria, sino que se cargan desde `../main/data/observational_abundances.csv`, mostrado en la tabla @tab:observational-input.

#figure(
  csv-table(
    obs,
    ("observable", "value", "sigma", "model_column", "source_key"),
    headers: ([observable], [valor], [$sigma$], [columna del modelo], [fuente]),
    small: true,
  ),
  caption: [Valores observacionales empleados en las comparaciones. D/H procede de absorción en cuásares pobres en metales @Cooke2018, $Y_p$ de regiones H II pobres en metales @Aver2021, #He3/H se usa solo como control débil @Steigman2010, y #Li7/H corresponde a la escala del plateau de Spite discutida en el problema cosmológico del litio @Fields2011.],
) <tab:observational-input>

= Modelo termodinámico

== Trayectoria forzada con referencia cosmológica

Una interpolación lineal en escala logarítmica entre los tres snapshots es suficiente para pasar por los puntos del enunciado, pero introduce quiebros artificiales en la pendiente de $T(t)$ y $rho(t)$. Por ello se utiliza una construcción ligeramente más física. Se parte de una referencia de Universo dominado por radiación, para la cual

$ H^2 = (8 pi G) / 3 rho_r, quad rho_r = pi^2 / 30 g_* T^4, $

y, usando $H approx 1/(2t)$ durante dominación radiativa, se obtiene la escala aproximada

$ T(t) prop t^(-1/2). $

Para la densidad bariónica se toma como referencia

$ rho_b(T) = m_u n_b(T) = m_u eta n_gamma(T), $

con

$ n_gamma(T) = (2 zeta(3)) / pi^2 (k_B T / (h c))^3. $

En la práctica, los notebooks construyen una referencia suave de este tipo y aplican una corrección multiplicativa mediante `PchipInterpolator` en escala logarítmica. Así se consigue que la trayectoria pase exactamente por los tres snapshots del enunciado, sin perder una forma global suave. El resultado se guarda en `../main/data/bbn_thermo_curve_forced.csv` y se comprueba en `../main/data/bbn_thermo_snapshot_check.csv`.

#figure(
  csv-table(
    thermo-check,
    ("label", "t_s", "T9_assignment", "T9_model", "rho_assignment_g_cm3", "rho_model_g_cm3", "eta10_effective"),
    headers: ([snapshot], [$t$ [s]], [$T_9$ enunciado], [$T_9$ modelo], [$rho$ enunciado], [$rho$ modelo], [$eta_10$ efectiva]),
    small: true,
  ),
  caption: [Control de la trayectoria termodinámica: la curva forzada reproduce los snapshots del enunciado y proporciona una $eta_10$ efectiva asociada a la densidad bariónica usada.],
) <tab:thermo-check>

La figura @fig:thermo-t9 compara la curva $T_9(t)$ utilizada con la referencia pura $T_9 prop t^(-1/2)$. La figura @fig:thermo-rho muestra la densidad bariónica forzada frente a la densidad de referencia calculada a partir de $rho_b=m_u eta n_gamma(T)$. Finalmente, @fig:eta-effective muestra la razón barión-fotón efectiva que correspondería a los snapshots impuestos.

#figure(
  image("../main/figures/fig_thermo_T9_forced_radiation.pdf", width: 85%),
  caption: [Historia de temperatura usada por la red. La curva pasa por los puntos del enunciado, pero mantiene como referencia la ley radiación-dominada $T_9 prop t^(-1/2)$.],
) <fig:thermo-t9>

#figure(
  image("../main/figures/fig_thermo_rho_forced_eta.pdf", width: 85%),
  caption: [Densidad bariónica impuesta. La referencia se calcula mediante $rho_b=m_u eta n_gamma(T)$ usando una $eta_10$ cosmológica de orden $6$, mientras que la curva final se fuerza a los snapshots del enunciado.],
) <fig:thermo-rho>

#figure(
  image("../main/figures/fig_eta10_effective_forced.pdf", width: 85%),
  caption: [Razón barión-fotón efectiva asociada a la trayectoria forzada. No se interpreta como una medida de $eta$, sino como diagnóstico de consistencia entre los snapshots del enunciado y la referencia cosmológica usada.],
) <fig:eta-effective>

== Comentario sobre $eta$

La razón barión-fotón medida por el fondo cósmico de microondas proporciona una restricción independiente de la densidad bariónica cosmológica. En este trabajo no se ajusta $eta$ para reproducir las abundancias, sino que se usa como referencia para construir una curva de densidad físicamente interpretable. La comparación entre `rho_forced` y `rho_eta_reference` permite ver hasta qué punto los snapshots del ejercicio son compatibles con una ley bariónica estándar basada en $rho_b prop eta T^3$.

= Red nuclear generada con `pynucastro`

== Especies incluidas

La red se genera con `pynucastro` a partir de una lista de núcleos ligeros. Para el objetivo de BBN mínima se incluyen neutrones, protones, deuterio, tritio, #He3, #He4, #Li7 y #Be7. La tabla @tab:nuclei se lee directamente desde `../main/data/bbn_network_nuclei.csv`.

#figure(
  csv-table(
    nuclei,
    ("index", "name", "A", "Z", "N", "A_minus_2Z"),
    headers: ([índice], [núcleo], [$A$], [$Z$], [$N$], [$A-2Z$]),
  ),
  caption: [Núcleos incluidos en la red BBN mínima generada con `pynucastro`.],
) <tab:nuclei>

La tabla @tab:rates muestra una selección de las reacciones exportadas por la red. El fichero completo se guarda en `../main/data/bbn_network_rates.csv`. Esta tabla es útil porque permite documentar exactamente qué canales nucleares entran en la simulación, evitando que la memoria dependa de una descripción verbal incompleta.

#figure(
  csv-table(
    first-n(rates, 18),
    ("index", "rate", "reactants", "products", "Q_MeV"),
    headers: ([índice], [rate], [reactivos], [productos], [$Q$ [MeV]]),
    small: true,
  ),
  caption: [Primeras reacciones de la red `pynucastro`. La tabla completa queda exportada en CSV.],
) <tab:rates>

#note-box[
  En esta versión se usa una red ligera y estable numéricamente. El propósito es estudiar la evolución básica y comparar escalas de abundancia. Un cálculo de precisión requeriría activar y controlar cuidadosamente reacciones inversas, fotodesintegraciones, tasas débiles y una evolución cosmológica completa, como hacen códigos especializados de BBN.
]

== Grafos de red en los snapshots del enunciado

Para visualizar la evolución de la composición se genera un grafo de red en los tres instantes del enunciado, usando únicamente la simulación que arranca en $t_0=50$ s. Los nodos representan núcleos en el plano $(Z, A-2Z)$, el color codifica $log_10 X_i$ y las flechas muestran conexiones nucleares extraídas de las tasas reales de `pynucastro`. Las tres figuras se exportan como PDF y se llaman directamente desde la memoria.

#figure(
  image("../main/figures/fig_network_t0_50_snapshot_50s.pdf", width: 82%),
  caption: [Grafo de red para la simulación $t_0=50$ s evaluado en el primer snapshot, $t=50$ s.],
) <fig:network-50>

#figure(
  image("../main/figures/fig_network_t0_50_snapshot_200s.pdf", width: 82%),
  caption: [Grafo de red para la simulación $t_0=50$ s evaluado en $t=200$ s.],
) <fig:network-200>

#figure(
  image("../main/figures/fig_network_t0_50_snapshot_1000s.pdf", width: 82%),
  caption: [Grafo de red para la simulación $t_0=50$ s evaluado en $t=1000$ s.],
) <fig:network-1000>

= Integración numérica

Se resuelven dos modelos con la misma red nuclear y la misma trayectoria termodinámica forzada, pero con dos tiempos iniciales distintos:

+ `t0_10s`: integración desde $t_0=10$ s.
+ `t0_50s`: integración desde $t_0=50$ s.

En ambos casos se impone la condición inicial

$ Y_n(t_0)=1/8, quad Y_p(t_0)=7/8, $

que corresponde a $N_n/N_p=1/7$. Las demás abundancias se inicializan a un valor nulo o despreciable. Las ecuaciones integradas tienen la forma

$ d Y_i / d t = F_i(Y, rho(t), T(t)), $

con $T(t)$ y $rho(t)$ dados por el modelo termodinámico descrito en la sección anterior. La integración se realiza con un método para sistemas rígidos (`BDF`/`LSODA`, según la versión del notebook), y se exporta la evolución completa de abundancias a:

```text
../main/data/bbn_evolution_t0_10s.csv
../main/data/bbn_evolution_t0_50s.csv
../main/data/bbn_evolution_all_models.csv
```

Además, se exportan tablas compactas para los snapshots y los valores finales, que son las que se leen a continuación desde la carpeta `../main/tables/`.

= Resultados de abundancias

== Evolución completa de fracciones másicas

Las figuras @fig:all-abundances-10 y @fig:all-abundances-50 muestran la evolución de todas las fracciones másicas $X_i=A_i Y_i$ para las dos condiciones iniciales. Estas gráficas permiten ver de forma global cuándo se forman las especies ligeras y qué núcleos dominan la composición final.

#figure(
  image("../main/figures/fig_all_abundances_t0_10s.pdf", width: 88%),
  caption: [Evolución de todas las fracciones másicas para la simulación iniciada en $t_0=10$ s.],
) <fig:all-abundances-10>

#figure(
  image("../main/figures/fig_all_abundances_t0_50s.pdf", width: 88%),
  caption: [Evolución de todas las fracciones másicas para la simulación iniciada en $t_0=50$ s.],
) <fig:all-abundances-50>

El resultado esperado de forma cualitativa es que una fracción importante de los neutrones termine en #He4, mientras que el exceso de protones permanezca como #H1. En cambio, D/H y #He3/H son residuos más sensibles a la historia térmica y a los canales de destrucción posteriores. El canal de masa 7, formado por #Li7 y #Be7, es todavía más delicado y se analiza por separado.

== Abundancias en los snapshots

La tabla @tab:snapshots-models resume las magnitudes principales en los tres tiempos pedidos por el enunciado para ambos modelos. Es la tabla central de la memoria, ya que responde directamente a la petición del ejercicio: evaluar la composición en $50$, $200$ y $1000$ s.

#figure(
  csv-table(
    snap-table,
    ("model", "t0_s", "t_s", "T9", "rho_g_cm3", "eta10", "D/H", "He3/H", "X_He4", "Li7/H", "Be7/H", "Li7_plus_Be7_over_H", "baryon_sum"),
    headers: ([modelo], [$t_0$], [$t$], [$T_9$], [$rho$], [$eta_10$], [D/H], [#He3/H], [$X(#He4)$], [#Li7/H], [#Be7/H], [(#Li7+#Be7)/H], [$sum X_i$]),
    small: true,
  ),
  caption: [Resultados de los dos modelos en los snapshots del enunciado. La columna $sum X_i$ sirve como control de conservación bariónica numérica.],
) <tab:snapshots-models>

== Valores finales

La tabla @tab:final-models recoge los valores finales usados para comparar con los observables. Esta tabla es más compacta que la anterior y se usa como entrada para las gráficas de comparación modelo-observación.

#figure(
  csv-table(
    final-table,
    ("model", "t0_s", "t_s", "T9", "rho_g_cm3", "eta10", "D/H", "He3/H", "X_He4", "Li7/H", "Be7/H", "Li7_plus_Be7_over_H", "baryon_sum"),
    headers: ([modelo], [$t_0$], [$t$ final], [$T_9$], [$rho$], [$eta_10$], [D/H], [#He3/H], [$X(#He4)$], [#Li7/H], [#Be7/H], [(#Li7+#Be7)/H], [$sum X_i$]),
    small: true,
  ),
  caption: [Valores finales de los modelos.],
) <tab:final-models>

= Comparación con abundancias observacionales

== Cocientes ligeros D/H, #He3/H y masa 7

La figura @fig:key-ratios compara la evolución temporal de los cocientes más relevantes con bandas observacionales horizontales. Se usa el mismo color para cada observable y distintos estilos de línea para distinguir los dos modelos: línea continua para $t_0=10$ s y línea discontinua para $t_0=50$ s.

#figure(
  image("../main/figures/fig_key_ratios_evolution_vs_observed.pdf", width: 88%),
  caption: [Evolución de D/H, #He3/H y $(#Li7+#Be7)/H$ para los dos modelos, comparada con las bandas observacionales usadas en el notebook.],
) <fig:key-ratios>

La comparación con D/H es especialmente relevante porque el deuterio es un barómetro muy sensible de la densidad bariónica. En BBN estándar, un mayor valor de $eta$ produce una destrucción más eficiente del deuterio y, por tanto, una menor D/H residual @Cooke2018 @Steigman2010. En cambio, #He3/H es menos limpio observacionalmente, por lo que aquí se usa sobre todo como control de orden de magnitud.

== Helio-4

El caso de #He4 se representa por separado en @fig:he4-evolution, ya que no es un cociente numérico pequeño sino una fracción másica de orden $0.25$. Su interés es que depende fuertemente de la razón neutrón-protón inicial y de la vida media del neutrón, pero es menos sensible a detalles finos de la red nuclear que el deuterio o el litio.

#figure(
  image("../main/figures/fig_he4_evolution_vs_observed.pdf", width: 84%),
  caption: [Evolución de la fracción másica $X(#He4)$ para los dos modelos y comparación con una banda observacional de $Y_p$.],
) <fig:he4-evolution>

La comparación final en formato de barras se muestra en @fig:he4-bars. Si el modelo conserva la idea física básica, debería producir una fracción de #He4 cercana a la estimación elemental $X_alpha approx 0.25$. Las diferencias respecto al valor observado deben interpretarse con cautela: este cálculo no incluye toda la física de BBN de precisión, sino que está construido como una red reducida para estudiar tendencias.

#figure(
  image("../main/figures/fig_final_he4_bar_models_vs_observed.pdf", width: 76%),
  caption: [Comparación final de $X(#He4)$ entre los dos modelos y el valor observacional adoptado.],
) <fig:he4-bars>

== Comparación final de cocientes

La figura @fig:final-ratios-bars condensa los cocientes D/H, #He3/H y #Li7/H en una gráfica de barras con escala logarítmica. En la tabla @tab:comparison-observed se exportan además todos los valores usados en la figura, incluyendo modelos y referencias observacionales.

#figure(
  image("../main/figures/fig_final_ratios_bar_models_vs_observed.pdf", width: 84%),
  caption: [Comparación final de cocientes de abundancia entre los dos modelos y los valores observacionales.],
) <fig:final-ratios-bars>

#figure(
  csv-table(
    comparison-table,
    ("observable", "kind", "value", "sigma", "source"),
    headers: ([observable], [caso], [valor], [$sigma$], [fuente]),
    small: true,
  ),
  caption: [Tabla completa de comparación entre modelos y valores observacionales.],
) <tab:comparison-observed>

= Diagnóstico del problema del litio

El canal de masa 7 se trata separadamente porque #Be7 producido durante BBN se transforma posteriormente en #Li7 por captura electrónica. Por tanto, la cantidad relevante que debe compararse con el litio observado no es únicamente #Li7 directo, sino aproximadamente

$ (#Li7 + #Be7) / H. $

La figura @fig:lithium-problem compara la evolución de esta cantidad con el valor observacional adoptado para #Li7/H. En BBN estándar, el problema del litio aparece porque las predicciones de #Li7 suelen exceder el plateau observado en estrellas pobres en metales @Fields2011. En este trabajo no se pretende resolver esa discrepancia, pero sí comprobar si el canal de masa 7 se comporta como un observable delicado y difícil de ajustar simultáneamente con D/H y $X(#He4)$.

#figure(
  image("../main/figures/fig_lithium7_problem_comparison.pdf", width: 86%),
  caption: [Diagnóstico del canal de masa 7. Se compara la producción de #Li7 directo y la suma #Li7+#Be7 frente a la banda observacional usada para #Li7/H.],
) <fig:lithium-problem>

= Mapas de abundancia y animaciones

Además de las curvas individuales, se generan mapas estáticos de $log_10 X_i$ en función del tiempo y de la especie. Estos mapas son útiles para identificar de un vistazo qué núcleos dominan en cada época y cuándo aparecen especies minoritarias.

#figure(
  image("../main/figures/fig_heatmap_abundances_t0_10s.pdf", width: 88%),
  caption: [Mapa de abundancias para el modelo $t_0=10$ s.],
) <fig:heatmap-10>

#figure(
  image("../main/figures/fig_heatmap_abundances_t0_50s.pdf", width: 88%),
  caption: [Mapa de abundancias para el modelo $t_0=50$ s.],
) <fig:heatmap-50>

También se exportan GIFs con la evolución temporal de las abundancias:

```text
../main/gifs/gif_abundances_t0_10s.gif
../main/gifs/gif_abundances_t0_50s.gif
```

// #figure(
//   image("../main/gifs/gif_abundances_t0_10s.gif", width: 80%),
//   caption: [Animación de la simulación.]
// )
// #figure(
//   image("../main/gifs/gif_abundances_t0_50s.gif", width: 80%),
//   caption: [Animación de la simulación.]
// )

Estos GIFs se consideran material suplementario, ya que el PDF final no conserva animación. Para la memoria impresa, las figuras @fig:heatmap-10 y @fig:heatmap-50 contienen la misma información de forma estática.

= Discusión

La simulación reproduce de forma explícita la lógica física básica de la BBN: la composición inicial está dominada por protones y neutrones libres, el deuterio actúa como cuello de botella, y la mayor parte de los neutrones disponibles termina ligada en #He4. Por ello, el primer control del cálculo es si la escala de $X(#He4)$ queda cerca de $0.25$, como sugiere la estimación elemental basada en $N_n/N_p=1/7$.

La diferencia entre los modelos $t_0=10$ s y $t_0=50$ s permite comprobar la sensibilidad a la elección del tiempo inicial. Si ambos modelos se aproximan a valores similares para #He4 pero difieren en especies minoritarias, la interpretación natural es que la producción de helio está dominada por la razón neutrón-protón impuesta, mientras que D/H, #He3/H y el canal de masa 7 dependen más de la historia térmica y de los detalles de la red.

El deuterio es el observable más útil para conectar el cálculo con la densidad bariónica. En la memoria se ha incluido la razón $eta_10$ efectiva asociada a la trayectoria forzada para dejar claro que los snapshots del enunciado no equivalen necesariamente a una trayectoria cosmológica estándar exacta. Esto es importante: al forzar simultáneamente $T(t)$ y $rho(t)$, se puede obtener una $eta_10$ efectiva que varía con el tiempo. Por tanto, la simulación debe entenderse como una trayectoria termodinámica impuesta para responder al ejercicio, no como una solución autoconsistente de Friedmann con $eta$ constante.

En cuanto al litio, el resultado debe interpretarse todavía con más cautela. La producción de masa 7 en BBN depende de canales minoritarios y de la competencia entre producción y destrucción de #Li7 y #Be7. Además, el problema cosmológico del litio es una discrepancia conocida entre predicciones estándar y observaciones estelares, no algo que una red reducida pueda resolver por sí sola. En esta memoria se usa como diagnóstico: si el canal de masa 7 no se ajusta tan bien como #He4 o D/H, eso no invalida automáticamente la simulación básica, sino que refleja precisamente la dificultad física de este observable @Fields2011 @Vangioni2000.

= Limitaciones del cálculo

Las principales limitaciones son:

+ La trayectoria termodinámica se fuerza a los snapshots del enunciado. Aunque se usa una referencia radiación-dominada para suavizar la curva, no se resuelven de forma completa las ecuaciones de Friedmann junto con la red nuclear.
+ La red es reducida y está diseñada para ser transparente y estable. Un cálculo de precisión debe incluir más canales, reacciones inversas, fotodesintegraciones y tasas débiles tratadas con mayor detalle.
+ La razón inicial $N_n/N_p=1/7$ se impone como condición inicial en ambos modelos. En una BBN completa, esta razón debería salir del freeze-out débil y del decaimiento beta del neutrón.
+ La comparación observacional usa valores de referencia representativos, pero no realiza un análisis estadístico de incertidumbres ni una inferencia bayesiana de parámetros.
+ El valor de #Li7/H observado procede de atmósferas estelares y está sujeto a posibles efectos de depleción, difusión y sistemáticos astrofísicos.

Estas limitaciones no eliminan la utilidad del ejercicio: el objetivo es construir una simulación mínima, documentada y reproducible, que permita relacionar directamente tasas nucleares, evolución térmica y abundancias finales.

= Conclusiones

Se ha construido una simulación de BBN con `pynucastro` usando una red de núcleos ligeros y una trayectoria termodinámica forzada a los tres snapshots del enunciado. La trayectoria se suaviza mediante una referencia de Universo dominado por radiación y una corrección en escala logarítmica, de modo que se evitan quiebros artificiales sin dejar de reproducir exactamente las condiciones pedidas.

Se han comparado dos modelos, iniciados en $t_0=10$ s y $t_0=50$ s, ambos con $N_n/N_p=1/7$. Para cada modelo se han exportado las abundancias completas, las tablas de snapshots y los valores finales. La memoria lee esos CSVs directamente con funciones de Typst, por lo que las tablas se actualizan automáticamente al volver a ejecutar los notebooks.

La producción de #He4 constituye el control físico principal: una razón inicial $N_n/N_p=1/7$ conduce naturalmente a una escala $X(#He4) approx 0.25$. El deuterio sirve como diagnóstico de la densidad bariónica efectiva, mientras que el canal #Li7+#Be7 permite ilustrar el problema cosmológico del litio y la sensibilidad de los elementos minoritarios a los detalles de la red nuclear.

En conjunto, el cálculo responde al ejercicio planteado: reproduce las condiciones termodinámicas de BBN, integra una red nuclear ligera, extrae abundancias en $50$, $200$ y $1000$ s, y compara los observables principales con valores experimentales o astrofísicos de referencia.

= Apéndice: ficheros generados

Los notebooks generan las siguientes salidas principales:

```text
../main/data/bbn_thermo_curve_forced.csv
../main/data/bbn_thermo_snapshot_check.csv
../main/data/bbn_network_nuclei.csv
../main/data/bbn_network_rates.csv
../main/data/bbn_evolution_t0_10s.csv
../main/data/bbn_evolution_t0_50s.csv
../main/data/bbn_snapshots_all_models.csv
../main/data/bbn_final_compact.csv
../main/data/comparison_models_vs_observed.csv
../main/tables/table_snapshots_models.csv
../main/tables/table_final_models.csv
../main/tables/table_comparison_models_vs_observed.csv
../main/figures/*.pdf
../main/gifs/*.gif
```

La memoria no copia esos ficheros: los lee directamente desde `../main/`, de modo que basta con volver a ejecutar los notebooks y recompilar `main.typ` para actualizar tablas y figuras.

= Referencias

#bibliography("references.bib", title: none)
