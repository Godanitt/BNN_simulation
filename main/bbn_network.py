import numba
import numpy as np
from pynucastro.constants import constants
from numba.experimental import jitclass

from pynucastro.rates import (TableIndex, TableInterpolator, TabularRate,
                              TempTableInterpolator, TemperatureTabularRate,
                              Tfactors)
from pynucastro.screening import PlasmaState, ScreenFactors

jn = 0
jp = 1
jd = 2
jt = 3
jhe3 = 4
jhe4 = 5
jli6 = 6
jli7 = 7
jbe7 = 8
nnuc = 9

A = np.zeros((nnuc), dtype=np.int32)

A[jn] = 1
A[jp] = 1
A[jd] = 2
A[jt] = 3
A[jhe3] = 3
A[jhe4] = 4
A[jli6] = 6
A[jli7] = 7
A[jbe7] = 7

Z = np.zeros((nnuc), dtype=np.int32)

Z[jn] = 0
Z[jp] = 1
Z[jd] = 1
Z[jt] = 1
Z[jhe3] = 2
Z[jhe4] = 2
Z[jli6] = 3
Z[jli7] = 3
Z[jbe7] = 4

# masses in ergs
mass = np.zeros((nnuc), dtype=np.float64)

mass[jn] = 0.001505349762871528
mass[jp] = 0.0015040963047307696
mass[jd] = 0.0030058819195053215
mass[jt] = 0.004501206494525079
mass[jhe3] = 0.004501176706825056
mass[jhe4] = 0.0059735574859708365
mass[jli6] = 0.008977078184259593
mass[jli7] = 0.010470810414554471
mass[jbe7] = 0.010472191322584432

names = []
names.append("n")
names.append("H1")
names.append("H2")
names.append("H3")
names.append("He3")
names.append("He4")
names.append("Li6")
names.append("Li7")
names.append("Be7")

def to_composition(Y):
    """Convert an array of molar fractions to a Composition object."""
    from pynucastro import Composition, Nucleus
    nuclei = [Nucleus.from_cache(name) for name in names]
    comp = Composition(nuclei)
    for i, nuc in enumerate(nuclei):
        comp.X[nuc] = Y[i] * A[i]
    return comp


def energy_release(dY):
    """return the energy release in erg/g (/s if dY is actually dY/dt)"""
    enuc = 0.0
    for i, y in enumerate(dY):
        enuc += y * mass[i]
    enuc *= -1*constants.N_A
    return enuc

@jitclass([
    ("n_to_p_reaclib", numba.float64),
    ("t_to_He3_reaclib", numba.float64),
    ("He3_to_t_reaclib", numba.float64),
    ("Be7_to_Li7_reaclib", numba.float64),
    ("n_p_to_d_reaclib", numba.float64),
    ("p_p_to_d_reaclib_beta_pos", numba.float64),
    ("p_p_to_d_reaclib_electron_capture", numba.float64),
    ("n_d_to_t_reaclib", numba.float64),
    ("p_d_to_He3_reaclib", numba.float64),
    ("d_d_to_He4_reaclib", numba.float64),
    ("He4_d_to_Li6_reaclib", numba.float64),
    ("p_t_to_He4_reaclib", numba.float64),
    ("He4_t_to_Li7_reaclib", numba.float64),
    ("n_He3_to_He4_reaclib", numba.float64),
    ("p_He3_to_He4_reaclib", numba.float64),
    ("He4_He3_to_Be7_reaclib", numba.float64),
    ("n_Li6_to_Li7_reaclib", numba.float64),
    ("p_Li6_to_Be7_reaclib", numba.float64),
    ("d_d_to_n_He3_reaclib", numba.float64),
    ("d_d_to_p_t_reaclib", numba.float64),
    ("d_t_to_n_He4_reaclib", numba.float64),
    ("n_He3_to_p_t_reaclib", numba.float64),
    ("d_He3_to_p_He4_reaclib", numba.float64),
    ("t_He3_to_d_He4_reaclib", numba.float64),
    ("n_Li6_to_He4_t_reaclib", numba.float64),
    ("p_Li6_to_He4_He3_reaclib", numba.float64),
    ("d_Li6_to_n_Be7_reaclib", numba.float64),
    ("d_Li6_to_p_Li7_reaclib", numba.float64),
    ("p_Li7_to_He4_He4_reaclib", numba.float64),
    ("n_Be7_to_p_Li7_reaclib", numba.float64),
    ("n_Be7_to_He4_He4_reaclib", numba.float64),
    ("t_t_to_n_n_He4_reaclib", numba.float64),
    ("t_He3_to_n_p_He4_reaclib", numba.float64),
    ("He3_He3_to_p_p_He4_reaclib", numba.float64),
    ("d_Li7_to_n_He4_He4_reaclib", numba.float64),
    ("d_Be7_to_p_He4_He4_reaclib", numba.float64),
    ("t_Li7_to_n_n_He4_He4_reaclib", numba.float64),
    ("He3_Li7_to_n_p_He4_He4_reaclib", numba.float64),
    ("t_Be7_to_n_p_He4_He4_reaclib", numba.float64),
    ("He3_Be7_to_p_p_He4_He4_reaclib", numba.float64),
    ("n_p_He4_to_Li6_reaclib", numba.float64),
    ("n_p_p_to_p_d_reaclib", numba.float64),
])
class RateEval:
    def __init__(self):
        self.n_to_p_reaclib = np.nan
        self.t_to_He3_reaclib = np.nan
        self.He3_to_t_reaclib = np.nan
        self.Be7_to_Li7_reaclib = np.nan
        self.n_p_to_d_reaclib = np.nan
        self.p_p_to_d_reaclib_beta_pos = np.nan
        self.p_p_to_d_reaclib_electron_capture = np.nan
        self.n_d_to_t_reaclib = np.nan
        self.p_d_to_He3_reaclib = np.nan
        self.d_d_to_He4_reaclib = np.nan
        self.He4_d_to_Li6_reaclib = np.nan
        self.p_t_to_He4_reaclib = np.nan
        self.He4_t_to_Li7_reaclib = np.nan
        self.n_He3_to_He4_reaclib = np.nan
        self.p_He3_to_He4_reaclib = np.nan
        self.He4_He3_to_Be7_reaclib = np.nan
        self.n_Li6_to_Li7_reaclib = np.nan
        self.p_Li6_to_Be7_reaclib = np.nan
        self.d_d_to_n_He3_reaclib = np.nan
        self.d_d_to_p_t_reaclib = np.nan
        self.d_t_to_n_He4_reaclib = np.nan
        self.n_He3_to_p_t_reaclib = np.nan
        self.d_He3_to_p_He4_reaclib = np.nan
        self.t_He3_to_d_He4_reaclib = np.nan
        self.n_Li6_to_He4_t_reaclib = np.nan
        self.p_Li6_to_He4_He3_reaclib = np.nan
        self.d_Li6_to_n_Be7_reaclib = np.nan
        self.d_Li6_to_p_Li7_reaclib = np.nan
        self.p_Li7_to_He4_He4_reaclib = np.nan
        self.n_Be7_to_p_Li7_reaclib = np.nan
        self.n_Be7_to_He4_He4_reaclib = np.nan
        self.t_t_to_n_n_He4_reaclib = np.nan
        self.t_He3_to_n_p_He4_reaclib = np.nan
        self.He3_He3_to_p_p_He4_reaclib = np.nan
        self.d_Li7_to_n_He4_He4_reaclib = np.nan
        self.d_Be7_to_p_He4_He4_reaclib = np.nan
        self.t_Li7_to_n_n_He4_He4_reaclib = np.nan
        self.He3_Li7_to_n_p_He4_He4_reaclib = np.nan
        self.t_Be7_to_n_p_He4_He4_reaclib = np.nan
        self.He3_Be7_to_p_p_He4_He4_reaclib = np.nan
        self.n_p_He4_to_Li6_reaclib = np.nan
        self.n_p_p_to_p_d_reaclib = np.nan

@numba.njit()
def ye(Y):
    return np.sum(Z * Y)/np.sum(A * Y)

@numba.njit()
def n_to_p_reaclib(rate_eval, tf, log_scor=0.0):
    # n --> p
    rate = 0.0

    # wc12w
    ln_set_rate =  -6.78161

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.n_to_p_reaclib = rate

@numba.njit()
def t_to_He3_reaclib(rate_eval, tf, log_scor=0.0):
    # t --> He3
    rate = 0.0

    # wc12w
    ln_set_rate =  -20.1456

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.t_to_He3_reaclib = rate

@numba.njit()
def He3_to_t_reaclib(rate_eval, tf, log_scor=0.0):
    # He3 --> t
    rate = 0.0

    #   ecw
    ln_set_rate =  -32.462 + -0.21338*tf.T9i + -0.821581*tf.T913i + 11.1241*tf.T913 \
                         + -0.577338*tf.T9 + 0.0290471*tf.T953 + -0.262705*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.He3_to_t_reaclib = rate

@numba.njit()
def Be7_to_Li7_reaclib(rate_eval, tf, log_scor=0.0):
    # Be7 --> Li7
    rate = 0.0

    #   ecw
    ln_set_rate =  -23.8328 + 3.02033*tf.T913 \
                         + -0.0742132*tf.T9 + -0.00792386*tf.T953 + -0.650113*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.Be7_to_Li7_reaclib = rate

@numba.njit()
def n_p_to_d_reaclib(rate_eval, tf, log_scor=0.0):
    # n + p --> d
    rate = 0.0

    # an06n
    ln_set_rate =  12.3687 + -2.70618*tf.T913 \
                         + 0.11718*tf.T9 + -0.00312788*tf.T953 + 0.469127*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    # an06n
    ln_set_rate =  10.7548 + -2.30472*tf.T913 \
                         + -0.887862*tf.T9 + 0.137663*tf.T953

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    # an06n
    ln_set_rate =  8.84688 + -0.0102082*tf.T913 \
                         + -0.0893959*tf.T9 + 0.00696704*tf.T953 + 1.0*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.n_p_to_d_reaclib = rate

@numba.njit()
def p_p_to_d_reaclib_beta_pos(rate_eval, tf, log_scor=0.0):
    # p + p --> d
    rate = 0.0

    # bet+w
    ln_set_rate =  -34.7863 + -3.51193*tf.T913i + 3.10086*tf.T913 \
                         + -0.198314*tf.T9 + 0.0126251*tf.T953 + -1.02517*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.p_p_to_d_reaclib_beta_pos = rate

@numba.njit()
def p_p_to_d_reaclib_electron_capture(rate_eval, tf, log_scor=0.0):
    # p + p --> d
    rate = 0.0

    #   ecw
    ln_set_rate =  -43.6499 + -0.00246064*tf.T9i + -2.7507*tf.T913i + -0.424877*tf.T913 \
                         + 0.015987*tf.T9 + -0.000690875*tf.T953 + -0.207625*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.p_p_to_d_reaclib_electron_capture = rate

@numba.njit()
def n_d_to_t_reaclib(rate_eval, tf, log_scor=0.0):
    # d + n --> t
    rate = 0.0

    # nk06n
    ln_set_rate =  6.60935 \
                         + 1.0*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    # nk06n
    ln_set_rate =  5.36598 \
                         + 0.075*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.n_d_to_t_reaclib = rate

@numba.njit()
def p_d_to_He3_reaclib(rate_eval, tf, log_scor=0.0):
    # d + p --> He3
    rate = 0.0

    # de04 
    ln_set_rate =  8.93525 + -3.7208*tf.T913i + 0.198654*tf.T913 \
                         + 0.333333*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    # de04n
    ln_set_rate =  7.52898 + -3.7208*tf.T913i + 0.871782*tf.T913 \
                         + -0.666667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.p_d_to_He3_reaclib = rate

@numba.njit()
def d_d_to_He4_reaclib(rate_eval, tf, log_scor=0.0):
    # d + d --> He4
    rate = 0.0

    # nacrn
    ln_set_rate =  3.78177 + -4.26166*tf.T913i + -0.119233*tf.T913 \
                         + 0.778829*tf.T9 + -0.0925203*tf.T953 + -0.666667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.d_d_to_He4_reaclib = rate

@numba.njit()
def He4_d_to_Li6_reaclib(rate_eval, tf, log_scor=0.0):
    # d + He4 --> Li6
    rate = 0.0

    # tu19r
    ln_set_rate =  4.12313 + -7.889*tf.T9i \
                         + -1.5*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    # tu19n
    ln_set_rate =  -0.676485 + 6.3911e-05*tf.T9i + -7.55198*tf.T913i + 5.77546*tf.T913 \
                         + -0.487854*tf.T9 + 0.032833*tf.T953 + -1.12305*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.He4_d_to_Li6_reaclib = rate

@numba.njit()
def p_t_to_He4_reaclib(rate_eval, tf, log_scor=0.0):
    # t + p --> He4
    rate = 0.0

    # cf88n
    ln_set_rate =  9.76526 + -3.869*tf.T913i + 1.45482*tf.T913 \
                         + 0.577246*tf.T9 + -0.112199*tf.T953 + -0.666667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.p_t_to_He4_reaclib = rate

@numba.njit()
def He4_t_to_Li7_reaclib(rate_eval, tf, log_scor=0.0):
    # t + He4 --> Li7
    rate = 0.0

    # de04 
    ln_set_rate =  13.6162 + -8.0805*tf.T913i + -0.217514*tf.T913 \
                         + -0.114859*tf.T9 + 0.0470043*tf.T953 + -0.666667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.He4_t_to_Li7_reaclib = rate

@numba.njit()
def n_He3_to_He4_reaclib(rate_eval, tf, log_scor=0.0):
    # He3 + n --> He4
    rate = 0.0

    # ka02n
    ln_set_rate =  9.04572 + -1.50147*tf.T913 \
                         + 1.0*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    # ka02n
    ln_set_rate =  5.51711

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.n_He3_to_He4_reaclib = rate

@numba.njit()
def p_He3_to_He4_reaclib(rate_eval, tf, log_scor=0.0):
    # He3 + p --> He4
    rate = 0.0

    # bet+w
    ln_set_rate =  -27.7611 + -4.30107e-12*tf.T9i + -6.141*tf.T913i + -1.93473e-09*tf.T913 \
                         + 2.04145e-10*tf.T9 + -1.80372e-11*tf.T953 + -0.666667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.p_He3_to_He4_reaclib = rate

@numba.njit()
def He4_He3_to_Be7_reaclib(rate_eval, tf, log_scor=0.0):
    # He3 + He4 --> Be7
    rate = 0.0

    # cd08n
    ln_set_rate =  17.7075 + -12.8271*tf.T913i + -3.8126*tf.T913 \
                         + 0.0942285*tf.T9 + -0.00301018*tf.T953 + 1.33333*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    # cd08n
    ln_set_rate =  15.6099 + -12.8271*tf.T913i + -0.0308225*tf.T913 \
                         + -0.654685*tf.T9 + 0.0896331*tf.T953 + -0.666667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.He4_He3_to_Be7_reaclib = rate

@numba.njit()
def n_Li6_to_Li7_reaclib(rate_eval, tf, log_scor=0.0):
    # Li6 + n --> Li7
    rate = 0.0

    # jz10n
    ln_set_rate =  9.04782

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.n_Li6_to_Li7_reaclib = rate

@numba.njit()
def p_Li6_to_Be7_reaclib(rate_eval, tf, log_scor=0.0):
    # Li6 + p --> Be7
    rate = 0.0

    # nacrn
    ln_set_rate =  14.2792 + -8.4372*tf.T913i + -0.515473*tf.T913 \
                         + 0.0285578*tf.T9 + 0.00879731*tf.T953 + -0.666667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.p_Li6_to_Be7_reaclib = rate

@numba.njit()
def d_d_to_n_He3_reaclib(rate_eval, tf, log_scor=0.0):
    # d + d --> n + He3
    rate = 0.0

    # gi17n
    ln_set_rate =  19.0876 + -0.00019002*tf.T9i + -4.2292*tf.T913i + 1.6932*tf.T913 \
                         + -0.0855529*tf.T9 + -1.35709e-25*tf.T953 + -0.734513*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.d_d_to_n_He3_reaclib = rate

@numba.njit()
def d_d_to_p_t_reaclib(rate_eval, tf, log_scor=0.0):
    # d + d --> p + t
    rate = 0.0

    # go17n
    ln_set_rate =  18.8052 + 4.36209e-05*tf.T9i + -4.32296*tf.T913i + 1.91572*tf.T913 \
                         + -0.081562*tf.T9 + -3.28804e-22*tf.T953 + -0.879518*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.d_d_to_p_t_reaclib = rate

@numba.njit()
def d_t_to_n_He4_reaclib(rate_eval, tf, log_scor=0.0):
    # t + d --> n + He4
    rate = 0.0

    # de04 
    ln_set_rate =  39.3457 + -4.5244*tf.T913i + -16.4028*tf.T913 \
                         + 1.73103*tf.T9 + -0.122966*tf.T953 + 2.31304*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    # de04 
    ln_set_rate =  25.1794 + -4.5244*tf.T913i + 0.350337*tf.T913 \
                         + 0.58747*tf.T9 + -8.84909*tf.T953 + -0.666667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.d_t_to_n_He4_reaclib = rate

@numba.njit()
def n_He3_to_p_t_reaclib(rate_eval, tf, log_scor=0.0):
    # He3 + n --> p + t
    rate = 0.0

    # de04 
    ln_set_rate =  19.2762 + 0.0438557*tf.T913 \
                         + -0.201527*tf.T9 + 0.0153433*tf.T953 + 1.0*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    # de04 
    ln_set_rate =  20.3787 + -0.332788*tf.T913 \
                         + -0.700485*tf.T9 + 0.0976521*tf.T953

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.n_He3_to_p_t_reaclib = rate

@numba.njit()
def d_He3_to_p_He4_reaclib(rate_eval, tf, log_scor=0.0):
    # He3 + d --> p + He4
    rate = 0.0

    # de04 
    ln_set_rate =  41.2969 + -7.182*tf.T913i + -17.1349*tf.T913 \
                         + 1.36908*tf.T9 + -0.0814423*tf.T953 + 3.35395*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    # de04 
    ln_set_rate =  24.6839 + -7.182*tf.T913i + 0.473288*tf.T913 \
                         + 1.46847*tf.T9 + -27.9603*tf.T953 + -0.666667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.d_He3_to_p_He4_reaclib = rate

@numba.njit()
def t_He3_to_d_He4_reaclib(rate_eval, tf, log_scor=0.0):
    # He3 + t --> d + He4
    rate = 0.0

    # cf88n
    ln_set_rate =  22.4207 + -7.733*tf.T913i + -0.133473*tf.T913 \
                         + -0.294412*tf.T9 + 0.0310968*tf.T953 + -0.666667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.t_He3_to_d_He4_reaclib = rate

@numba.njit()
def n_Li6_to_He4_t_reaclib(rate_eval, tf, log_scor=0.0):
    # Li6 + n --> He4 + t
    rate = 0.0

    # cf88r
    ln_set_rate =  21.665 + -2.39128*tf.T9i \
                         + -1.5*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    # cf88n
    ln_set_rate =  18.9496 + -0.001281*tf.T9i

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.n_Li6_to_He4_t_reaclib = rate

@numba.njit()
def p_Li6_to_He4_He3_reaclib(rate_eval, tf, log_scor=0.0):
    # Li6 + p --> He4 + He3
    rate = 0.0

    # pt05n
    ln_set_rate =  24.3475 + -8.39481*tf.T913i + -0.165254*tf.T913 \
                         + -0.16936*tf.T9 + 0.0533676*tf.T953 + -0.666667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.p_Li6_to_He4_He3_reaclib = rate

@numba.njit()
def d_Li6_to_n_Be7_reaclib(rate_eval, tf, log_scor=0.0):
    # Li6 + d --> n + Be7
    rate = 0.0

    # mafon
    ln_set_rate =  28.0095 + -4.77456e-12*tf.T9i + -10.259*tf.T913i + -2.01559e-09*tf.T913 \
                         + 1.99542e-10*tf.T9 + -1.65595e-11*tf.T953 + -0.666667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.d_Li6_to_n_Be7_reaclib = rate

@numba.njit()
def d_Li6_to_p_Li7_reaclib(rate_eval, tf, log_scor=0.0):
    # Li6 + d --> p + Li7
    rate = 0.0

    # mafon
    ln_set_rate =  28.0231 + -10.135*tf.T913i \
                         + -0.666667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.d_Li6_to_p_Li7_reaclib = rate

@numba.njit()
def p_Li7_to_He4_He4_reaclib(rate_eval, tf, log_scor=0.0):
    # Li7 + p --> He4 + He4
    rate = 0.0

    # de04r
    ln_set_rate =  21.8999 + -26.1527*tf.T9i \
                         + -1.5*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    # de04 
    ln_set_rate =  20.4438 + -8.4727*tf.T913i + 0.297934*tf.T913 \
                         + 0.0582335*tf.T9 + -0.00413383*tf.T953 + -0.666667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    # de04r
    ln_set_rate =  14.2538 + -4.478*tf.T9i \
                         + -1.5*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    # de04 
    ln_set_rate =  11.9576 + -8.4727*tf.T913i + 0.417943*tf.T913 \
                         + 5.34565*tf.T9 + -4.8684*tf.T953 + -0.666667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.p_Li7_to_He4_He4_reaclib = rate

@numba.njit()
def n_Be7_to_p_Li7_reaclib(rate_eval, tf, log_scor=0.0):
    # Be7 + n --> p + Li7
    rate = 0.0

    # db18 
    ln_set_rate =  21.7899 + 0.000728098*tf.T9i + -0.30254*tf.T913i + -0.3602*tf.T913 \
                         + 0.17472*tf.T9 + -0.0223*tf.T953 + -0.4581*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.n_Be7_to_p_Li7_reaclib = rate

@numba.njit()
def n_Be7_to_He4_He4_reaclib(rate_eval, tf, log_scor=0.0):
    # Be7 + n --> He4 + He4
    rate = 0.0

    #  wagn
    ln_set_rate =  18.1614 + -0.00210045*tf.T913 \
                         + 0.000176541*tf.T9 + -1.36797e-05*tf.T953 + 1.00083*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.n_Be7_to_He4_He4_reaclib = rate

@numba.njit()
def t_t_to_n_n_He4_reaclib(rate_eval, tf, log_scor=0.0):
    # t + t --> n + n + He4
    rate = 0.0

    # cf88r
    ln_set_rate =  21.2361 + -4.872*tf.T913i + -1.72398*tf.T913 \
                         + 0.684775*tf.T9 + -0.0702582*tf.T953 + 0.333333*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    # cf88n
    ln_set_rate =  21.2361 + -4.872*tf.T913i + -0.0328579*tf.T913 \
                         + -1.13588*tf.T9 + 0.250064*tf.T953 + -0.666667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.t_t_to_n_n_He4_reaclib = rate

@numba.njit()
def t_He3_to_n_p_He4_reaclib(rate_eval, tf, log_scor=0.0):
    # He3 + t --> n + p + He4
    rate = 0.0

    # cf88n
    ln_set_rate =  22.7658 + -7.733*tf.T913i + -0.118902*tf.T913 \
                         + -0.267393*tf.T9 + 0.0275387*tf.T953 + -0.666667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.t_He3_to_n_p_He4_reaclib = rate

@numba.njit()
def He3_He3_to_p_p_He4_reaclib(rate_eval, tf, log_scor=0.0):
    # He3 + He3 --> p + p + He4
    rate = 0.0

    # nacrn
    ln_set_rate =  24.7788 + -12.277*tf.T913i + -0.103699*tf.T913 \
                         + -0.0649967*tf.T9 + 0.0168191*tf.T953 + -0.666667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.He3_He3_to_p_p_He4_reaclib = rate

@numba.njit()
def d_Li7_to_n_He4_He4_reaclib(rate_eval, tf, log_scor=0.0):
    # Li7 + d --> n + He4 + He4
    rate = 0.0

    # cf88n
    ln_set_rate =  26.4 + -10.259*tf.T913i \
                         + -0.666667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.d_Li7_to_n_He4_He4_reaclib = rate

@numba.njit()
def d_Be7_to_p_He4_He4_reaclib(rate_eval, tf, log_scor=0.0):
    # Be7 + d --> p + He4 + He4
    rate = 0.0

    # cf88n
    ln_set_rate =  27.6987 + -12.428*tf.T913i \
                         + -0.666667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.d_Be7_to_p_He4_He4_reaclib = rate

@numba.njit()
def t_Li7_to_n_n_He4_He4_reaclib(rate_eval, tf, log_scor=0.0):
    # Li7 + t --> n + n + He4 + He4
    rate = 0.0

    # mafon
    ln_set_rate =  27.5043 + -5.31692e-12*tf.T9i + -11.333*tf.T913i + -2.24192e-09*tf.T913 \
                         + 2.21773e-10*tf.T9 + -1.83941e-11*tf.T953 + -0.666667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.t_Li7_to_n_n_He4_He4_reaclib = rate

@numba.njit()
def He3_Li7_to_n_p_He4_He4_reaclib(rate_eval, tf, log_scor=0.0):
    # Li7 + He3 --> n + p + He4 + He4
    rate = 0.0

    # mafon
    ln_set_rate =  30.038 + -4.24733e-12*tf.T9i + -17.989*tf.T913i + -1.57523e-09*tf.T913 \
                         + 1.45934e-10*tf.T9 + -1.15341e-11*tf.T953 + -0.666667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.He3_Li7_to_n_p_He4_He4_reaclib = rate

@numba.njit()
def t_Be7_to_n_p_He4_He4_reaclib(rate_eval, tf, log_scor=0.0):
    # Be7 + t --> n + p + He4 + He4
    rate = 0.0

    # mafon
    ln_set_rate =  28.6992 + -6.9004e-12*tf.T9i + -13.792*tf.T913i + -2.92021e-09*tf.T913 \
                         + 2.89378e-10*tf.T9 + -2.40287e-11*tf.T953 + -0.666667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.t_Be7_to_n_p_He4_He4_reaclib = rate

@numba.njit()
def He3_Be7_to_p_p_He4_He4_reaclib(rate_eval, tf, log_scor=0.0):
    # Be7 + He3 --> p + p + He4 + He4
    rate = 0.0

    # mafon
    ln_set_rate =  31.7435 + -5.45213e-12*tf.T9i + -21.793*tf.T913i + -1.98126e-09*tf.T913 \
                         + 1.84204e-10*tf.T9 + -1.46403e-11*tf.T953 + -0.666667*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.He3_Be7_to_p_p_He4_He4_reaclib = rate

@numba.njit()
def n_p_He4_to_Li6_reaclib(rate_eval, tf, log_scor=0.0):
    # n + p + He4 --> Li6
    rate = 0.0

    # cf88r
    ln_set_rate =  -12.2851 + -19.353*tf.T9i + 1.44987*tf.T913i + -1.42759*tf.T913 \
                         + 0.0454035*tf.T9 + 0.00471161*tf.T953 + -1.0*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.n_p_He4_to_Li6_reaclib = rate

@numba.njit()
def n_p_p_to_p_d_reaclib(rate_eval, tf, log_scor=0.0):
    # n + p + p --> p + d
    rate = 0.0

    # cf88n
    ln_set_rate =  -4.24034 + -3.72*tf.T913i + 0.946313*tf.T913 \
                         + 0.105406*tf.T9 + -0.0149431*tf.T953 + -1.5*tf.lnT9

    ln_set_rate += log_scor
    set_rate = np.exp(ln_set_rate)
    rate += set_rate

    rate_eval.n_p_p_to_p_d_reaclib = rate

def rhs(t, Y, rho, T, screen_func=None):
    return rhs_eq(t, Y, rho, T, screen_func)

@numba.njit()
def do_rate_eval(t, Y, rho, T, screen_func):

    tf = Tfactors(T)
    rate_eval = RateEval()

    log_scor_d_t = 0.0
    log_scor_p_Li7 = 0.0
    log_scor_He3_He3 = 0.0
    log_scor_d_Li7 = 0.0
    log_scor_He3_He4 = 0.0
    log_scor_He3_Be7 = 0.0
    log_scor_p_p = 0.0
    log_scor_t_He3 = 0.0
    log_scor_p_d = 0.0
    log_scor_t_He4 = 0.0
    log_scor_p_He3 = 0.0
    log_scor_t_Be7 = 0.0
    log_scor_d_d = 0.0
    log_scor_He3_Li7 = 0.0
    log_scor_d_He3 = 0.0
    log_scor_p_He4 = 0.0
    log_scor_p_Li6 = 0.0
    log_scor_d_He4 = 0.0
    log_scor_d_Li6 = 0.0
    log_scor_d_Be7 = 0.0
    log_scor_t_t = 0.0
    log_scor_t_Li7 = 0.0
    log_scor_p_t = 0.0

    if screen_func is not None:
        plasma_state = PlasmaState(T, rho, Y, Z)

        scn_fac = ScreenFactors(1, 2, 1, 3)
        log_scor_d_t = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 1, 3, 7)
        log_scor_p_Li7 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(2, 3, 2, 3)
        log_scor_He3_He3 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 2, 3, 7)
        log_scor_d_Li7 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(2, 3, 2, 4)
        log_scor_He3_He4 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(2, 3, 4, 7)
        log_scor_He3_Be7 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 1, 1, 1)
        log_scor_p_p = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 3, 2, 3)
        log_scor_t_He3 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 1, 1, 2)
        log_scor_p_d = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 3, 2, 4)
        log_scor_t_He4 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 1, 2, 3)
        log_scor_p_He3 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 3, 4, 7)
        log_scor_t_Be7 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 2, 1, 2)
        log_scor_d_d = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(2, 3, 3, 7)
        log_scor_He3_Li7 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 2, 2, 3)
        log_scor_d_He3 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 1, 2, 4)
        log_scor_p_He4 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 1, 3, 6)
        log_scor_p_Li6 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 2, 2, 4)
        log_scor_d_He4 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 2, 3, 6)
        log_scor_d_Li6 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 2, 4, 7)
        log_scor_d_Be7 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 3, 1, 3)
        log_scor_t_t = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 3, 3, 7)
        log_scor_t_Li7 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 1, 1, 3)
        log_scor_p_t = screen_func(plasma_state, scn_fac)

    # reaclib rates
    n_to_p_reaclib(rate_eval, tf)
    t_to_He3_reaclib(rate_eval, tf)
    He3_to_t_reaclib(rate_eval, tf)
    Be7_to_Li7_reaclib(rate_eval, tf)
    n_p_to_d_reaclib(rate_eval, tf)
    p_p_to_d_reaclib_beta_pos(rate_eval, tf, log_scor=log_scor_p_p)
    p_p_to_d_reaclib_electron_capture(rate_eval, tf, log_scor=log_scor_p_p)
    n_d_to_t_reaclib(rate_eval, tf)
    p_d_to_He3_reaclib(rate_eval, tf, log_scor=log_scor_p_d)
    d_d_to_He4_reaclib(rate_eval, tf, log_scor=log_scor_d_d)
    He4_d_to_Li6_reaclib(rate_eval, tf, log_scor=log_scor_d_He4)
    p_t_to_He4_reaclib(rate_eval, tf, log_scor=log_scor_p_t)
    He4_t_to_Li7_reaclib(rate_eval, tf, log_scor=log_scor_t_He4)
    n_He3_to_He4_reaclib(rate_eval, tf)
    p_He3_to_He4_reaclib(rate_eval, tf, log_scor=log_scor_p_He3)
    He4_He3_to_Be7_reaclib(rate_eval, tf, log_scor=log_scor_He3_He4)
    n_Li6_to_Li7_reaclib(rate_eval, tf)
    p_Li6_to_Be7_reaclib(rate_eval, tf, log_scor=log_scor_p_Li6)
    d_d_to_n_He3_reaclib(rate_eval, tf, log_scor=log_scor_d_d)
    d_d_to_p_t_reaclib(rate_eval, tf, log_scor=log_scor_d_d)
    d_t_to_n_He4_reaclib(rate_eval, tf, log_scor=log_scor_d_t)
    n_He3_to_p_t_reaclib(rate_eval, tf)
    d_He3_to_p_He4_reaclib(rate_eval, tf, log_scor=log_scor_d_He3)
    t_He3_to_d_He4_reaclib(rate_eval, tf, log_scor=log_scor_t_He3)
    n_Li6_to_He4_t_reaclib(rate_eval, tf)
    p_Li6_to_He4_He3_reaclib(rate_eval, tf, log_scor=log_scor_p_Li6)
    d_Li6_to_n_Be7_reaclib(rate_eval, tf, log_scor=log_scor_d_Li6)
    d_Li6_to_p_Li7_reaclib(rate_eval, tf, log_scor=log_scor_d_Li6)
    p_Li7_to_He4_He4_reaclib(rate_eval, tf, log_scor=log_scor_p_Li7)
    n_Be7_to_p_Li7_reaclib(rate_eval, tf)
    n_Be7_to_He4_He4_reaclib(rate_eval, tf)
    t_t_to_n_n_He4_reaclib(rate_eval, tf, log_scor=log_scor_t_t)
    t_He3_to_n_p_He4_reaclib(rate_eval, tf, log_scor=log_scor_t_He3)
    He3_He3_to_p_p_He4_reaclib(rate_eval, tf, log_scor=log_scor_He3_He3)
    d_Li7_to_n_He4_He4_reaclib(rate_eval, tf, log_scor=log_scor_d_Li7)
    d_Be7_to_p_He4_He4_reaclib(rate_eval, tf, log_scor=log_scor_d_Be7)
    t_Li7_to_n_n_He4_He4_reaclib(rate_eval, tf, log_scor=log_scor_t_Li7)
    He3_Li7_to_n_p_He4_He4_reaclib(rate_eval, tf, log_scor=log_scor_He3_Li7)
    t_Be7_to_n_p_He4_He4_reaclib(rate_eval, tf, log_scor=log_scor_t_Be7)
    He3_Be7_to_p_p_He4_He4_reaclib(rate_eval, tf, log_scor=log_scor_He3_Be7)
    n_p_He4_to_Li6_reaclib(rate_eval, tf, log_scor=log_scor_p_He4)
    n_p_p_to_p_d_reaclib(rate_eval, tf, log_scor=log_scor_p_p)

    return rate_eval

@numba.njit()
def rhs_eq(t, Y, rho, T, screen_func):

    rate_eval = do_rate_eval(t, Y, rho, T, screen_func)
    dYdt = np.zeros((nnuc), dtype=np.float64)

    dYdt[jn] = (
          -Y[jn]*rate_eval.n_to_p_reaclib  +
          -rho*Y[jn]*Y[jp]*rate_eval.n_p_to_d_reaclib  +
          -rho*Y[jn]*Y[jd]*rate_eval.n_d_to_t_reaclib  +
          -rho*Y[jn]*Y[jhe3]*rate_eval.n_He3_to_He4_reaclib  +
          -rho*Y[jn]*Y[jli6]*rate_eval.n_Li6_to_Li7_reaclib  +
          +5.00000000000000e-01*rho*Y[jd]**2*rate_eval.d_d_to_n_He3_reaclib  +
          +rho*Y[jd]*Y[jt]*rate_eval.d_t_to_n_He4_reaclib  +
          -rho*Y[jn]*Y[jhe3]*rate_eval.n_He3_to_p_t_reaclib  +
          -rho*Y[jn]*Y[jli6]*rate_eval.n_Li6_to_He4_t_reaclib  +
          +rho*Y[jd]*Y[jli6]*rate_eval.d_Li6_to_n_Be7_reaclib  +
          -rho*Y[jn]*Y[jbe7]*rate_eval.n_Be7_to_p_Li7_reaclib  +
          -rho*Y[jn]*Y[jbe7]*rate_eval.n_Be7_to_He4_He4_reaclib  +
          + 2*5.00000000000000e-01*rho*Y[jt]**2*rate_eval.t_t_to_n_n_He4_reaclib  +
          +rho*Y[jt]*Y[jhe3]*rate_eval.t_He3_to_n_p_He4_reaclib  +
          +rho*Y[jd]*Y[jli7]*rate_eval.d_Li7_to_n_He4_He4_reaclib  +
          + 2*rho*Y[jt]*Y[jli7]*rate_eval.t_Li7_to_n_n_He4_He4_reaclib  +
          +rho*Y[jhe3]*Y[jli7]*rate_eval.He3_Li7_to_n_p_He4_He4_reaclib  +
          +rho*Y[jt]*Y[jbe7]*rate_eval.t_Be7_to_n_p_He4_He4_reaclib  +
          -rho**2*Y[jn]*Y[jp]*Y[jhe4]*rate_eval.n_p_He4_to_Li6_reaclib  +
          -5.00000000000000e-01*rho**2*Y[jn]*Y[jp]**2*rate_eval.n_p_p_to_p_d_reaclib
       )

    dYdt[jp] = (
          +Y[jn]*rate_eval.n_to_p_reaclib  +
          -rho*Y[jn]*Y[jp]*rate_eval.n_p_to_d_reaclib  +
          + -2*5.00000000000000e-01*rho*Y[jp]**2*rate_eval.p_p_to_d_reaclib_beta_pos  +
          + -2*5.00000000000000e-01*rho**2*ye(Y)*Y[jp]**2*rate_eval.p_p_to_d_reaclib_electron_capture  +
          -rho*Y[jp]*Y[jd]*rate_eval.p_d_to_He3_reaclib  +
          -rho*Y[jp]*Y[jt]*rate_eval.p_t_to_He4_reaclib  +
          -rho*Y[jp]*Y[jhe3]*rate_eval.p_He3_to_He4_reaclib  +
          -rho*Y[jp]*Y[jli6]*rate_eval.p_Li6_to_Be7_reaclib  +
          +5.00000000000000e-01*rho*Y[jd]**2*rate_eval.d_d_to_p_t_reaclib  +
          +rho*Y[jn]*Y[jhe3]*rate_eval.n_He3_to_p_t_reaclib  +
          +rho*Y[jd]*Y[jhe3]*rate_eval.d_He3_to_p_He4_reaclib  +
          -rho*Y[jp]*Y[jli6]*rate_eval.p_Li6_to_He4_He3_reaclib  +
          +rho*Y[jd]*Y[jli6]*rate_eval.d_Li6_to_p_Li7_reaclib  +
          -rho*Y[jp]*Y[jli7]*rate_eval.p_Li7_to_He4_He4_reaclib  +
          +rho*Y[jn]*Y[jbe7]*rate_eval.n_Be7_to_p_Li7_reaclib  +
          +rho*Y[jt]*Y[jhe3]*rate_eval.t_He3_to_n_p_He4_reaclib  +
          + 2*5.00000000000000e-01*rho*Y[jhe3]**2*rate_eval.He3_He3_to_p_p_He4_reaclib  +
          +rho*Y[jd]*Y[jbe7]*rate_eval.d_Be7_to_p_He4_He4_reaclib  +
          +rho*Y[jhe3]*Y[jli7]*rate_eval.He3_Li7_to_n_p_He4_He4_reaclib  +
          +rho*Y[jt]*Y[jbe7]*rate_eval.t_Be7_to_n_p_He4_He4_reaclib  +
          + 2*rho*Y[jhe3]*Y[jbe7]*rate_eval.He3_Be7_to_p_p_He4_He4_reaclib  +
          -rho**2*Y[jn]*Y[jp]*Y[jhe4]*rate_eval.n_p_He4_to_Li6_reaclib  +
          -5.00000000000000e-01*rho**2*Y[jn]*Y[jp]**2*rate_eval.n_p_p_to_p_d_reaclib
       )

    dYdt[jd] = (
          +rho*Y[jn]*Y[jp]*rate_eval.n_p_to_d_reaclib  +
          +5.00000000000000e-01*rho*Y[jp]**2*rate_eval.p_p_to_d_reaclib_beta_pos  +
          +5.00000000000000e-01*rho**2*ye(Y)*Y[jp]**2*rate_eval.p_p_to_d_reaclib_electron_capture  +
          -rho*Y[jn]*Y[jd]*rate_eval.n_d_to_t_reaclib  +
          -rho*Y[jp]*Y[jd]*rate_eval.p_d_to_He3_reaclib  +
          + -2*5.00000000000000e-01*rho*Y[jd]**2*rate_eval.d_d_to_He4_reaclib  +
          -rho*Y[jd]*Y[jhe4]*rate_eval.He4_d_to_Li6_reaclib  +
          + -2*5.00000000000000e-01*rho*Y[jd]**2*rate_eval.d_d_to_n_He3_reaclib  +
          + -2*5.00000000000000e-01*rho*Y[jd]**2*rate_eval.d_d_to_p_t_reaclib  +
          -rho*Y[jd]*Y[jt]*rate_eval.d_t_to_n_He4_reaclib  +
          -rho*Y[jd]*Y[jhe3]*rate_eval.d_He3_to_p_He4_reaclib  +
          +rho*Y[jt]*Y[jhe3]*rate_eval.t_He3_to_d_He4_reaclib  +
          -rho*Y[jd]*Y[jli6]*rate_eval.d_Li6_to_n_Be7_reaclib  +
          -rho*Y[jd]*Y[jli6]*rate_eval.d_Li6_to_p_Li7_reaclib  +
          -rho*Y[jd]*Y[jli7]*rate_eval.d_Li7_to_n_He4_He4_reaclib  +
          -rho*Y[jd]*Y[jbe7]*rate_eval.d_Be7_to_p_He4_He4_reaclib  +
          +5.00000000000000e-01*rho**2*Y[jn]*Y[jp]**2*rate_eval.n_p_p_to_p_d_reaclib
       )

    dYdt[jt] = (
          ( -Y[jt]*rate_eval.t_to_He3_reaclib +rho*ye(Y)*Y[jhe3]*rate_eval.He3_to_t_reaclib ) +
          +rho*Y[jn]*Y[jd]*rate_eval.n_d_to_t_reaclib  +
          -rho*Y[jp]*Y[jt]*rate_eval.p_t_to_He4_reaclib  +
          -rho*Y[jt]*Y[jhe4]*rate_eval.He4_t_to_Li7_reaclib  +
          +5.00000000000000e-01*rho*Y[jd]**2*rate_eval.d_d_to_p_t_reaclib  +
          -rho*Y[jd]*Y[jt]*rate_eval.d_t_to_n_He4_reaclib  +
          +rho*Y[jn]*Y[jhe3]*rate_eval.n_He3_to_p_t_reaclib  +
          -rho*Y[jt]*Y[jhe3]*rate_eval.t_He3_to_d_He4_reaclib  +
          +rho*Y[jn]*Y[jli6]*rate_eval.n_Li6_to_He4_t_reaclib  +
          + -2*5.00000000000000e-01*rho*Y[jt]**2*rate_eval.t_t_to_n_n_He4_reaclib  +
          -rho*Y[jt]*Y[jhe3]*rate_eval.t_He3_to_n_p_He4_reaclib  +
          -rho*Y[jt]*Y[jli7]*rate_eval.t_Li7_to_n_n_He4_He4_reaclib  +
          -rho*Y[jt]*Y[jbe7]*rate_eval.t_Be7_to_n_p_He4_He4_reaclib
       )

    dYdt[jhe3] = (
          ( +Y[jt]*rate_eval.t_to_He3_reaclib -rho*ye(Y)*Y[jhe3]*rate_eval.He3_to_t_reaclib ) +
          +rho*Y[jp]*Y[jd]*rate_eval.p_d_to_He3_reaclib  +
          -rho*Y[jn]*Y[jhe3]*rate_eval.n_He3_to_He4_reaclib  +
          -rho*Y[jp]*Y[jhe3]*rate_eval.p_He3_to_He4_reaclib  +
          -rho*Y[jhe3]*Y[jhe4]*rate_eval.He4_He3_to_Be7_reaclib  +
          +5.00000000000000e-01*rho*Y[jd]**2*rate_eval.d_d_to_n_He3_reaclib  +
          -rho*Y[jn]*Y[jhe3]*rate_eval.n_He3_to_p_t_reaclib  +
          -rho*Y[jd]*Y[jhe3]*rate_eval.d_He3_to_p_He4_reaclib  +
          -rho*Y[jt]*Y[jhe3]*rate_eval.t_He3_to_d_He4_reaclib  +
          +rho*Y[jp]*Y[jli6]*rate_eval.p_Li6_to_He4_He3_reaclib  +
          -rho*Y[jt]*Y[jhe3]*rate_eval.t_He3_to_n_p_He4_reaclib  +
          + -2*5.00000000000000e-01*rho*Y[jhe3]**2*rate_eval.He3_He3_to_p_p_He4_reaclib  +
          -rho*Y[jhe3]*Y[jli7]*rate_eval.He3_Li7_to_n_p_He4_He4_reaclib  +
          -rho*Y[jhe3]*Y[jbe7]*rate_eval.He3_Be7_to_p_p_He4_He4_reaclib
       )

    dYdt[jhe4] = (
          +5.00000000000000e-01*rho*Y[jd]**2*rate_eval.d_d_to_He4_reaclib  +
          -rho*Y[jd]*Y[jhe4]*rate_eval.He4_d_to_Li6_reaclib  +
          +rho*Y[jp]*Y[jt]*rate_eval.p_t_to_He4_reaclib  +
          -rho*Y[jt]*Y[jhe4]*rate_eval.He4_t_to_Li7_reaclib  +
          +rho*Y[jn]*Y[jhe3]*rate_eval.n_He3_to_He4_reaclib  +
          +rho*Y[jp]*Y[jhe3]*rate_eval.p_He3_to_He4_reaclib  +
          -rho*Y[jhe3]*Y[jhe4]*rate_eval.He4_He3_to_Be7_reaclib  +
          +rho*Y[jd]*Y[jt]*rate_eval.d_t_to_n_He4_reaclib  +
          +rho*Y[jd]*Y[jhe3]*rate_eval.d_He3_to_p_He4_reaclib  +
          +rho*Y[jt]*Y[jhe3]*rate_eval.t_He3_to_d_He4_reaclib  +
          +rho*Y[jn]*Y[jli6]*rate_eval.n_Li6_to_He4_t_reaclib  +
          +rho*Y[jp]*Y[jli6]*rate_eval.p_Li6_to_He4_He3_reaclib  +
          + 2*rho*Y[jp]*Y[jli7]*rate_eval.p_Li7_to_He4_He4_reaclib  +
          + 2*rho*Y[jn]*Y[jbe7]*rate_eval.n_Be7_to_He4_He4_reaclib  +
          +5.00000000000000e-01*rho*Y[jt]**2*rate_eval.t_t_to_n_n_He4_reaclib  +
          +rho*Y[jt]*Y[jhe3]*rate_eval.t_He3_to_n_p_He4_reaclib  +
          +5.00000000000000e-01*rho*Y[jhe3]**2*rate_eval.He3_He3_to_p_p_He4_reaclib  +
          + 2*rho*Y[jd]*Y[jli7]*rate_eval.d_Li7_to_n_He4_He4_reaclib  +
          + 2*rho*Y[jd]*Y[jbe7]*rate_eval.d_Be7_to_p_He4_He4_reaclib  +
          + 2*rho*Y[jt]*Y[jli7]*rate_eval.t_Li7_to_n_n_He4_He4_reaclib  +
          + 2*rho*Y[jhe3]*Y[jli7]*rate_eval.He3_Li7_to_n_p_He4_He4_reaclib  +
          + 2*rho*Y[jt]*Y[jbe7]*rate_eval.t_Be7_to_n_p_He4_He4_reaclib  +
          + 2*rho*Y[jhe3]*Y[jbe7]*rate_eval.He3_Be7_to_p_p_He4_He4_reaclib  +
          -rho**2*Y[jn]*Y[jp]*Y[jhe4]*rate_eval.n_p_He4_to_Li6_reaclib
       )

    dYdt[jli6] = (
          +rho*Y[jd]*Y[jhe4]*rate_eval.He4_d_to_Li6_reaclib  +
          -rho*Y[jn]*Y[jli6]*rate_eval.n_Li6_to_Li7_reaclib  +
          -rho*Y[jp]*Y[jli6]*rate_eval.p_Li6_to_Be7_reaclib  +
          -rho*Y[jn]*Y[jli6]*rate_eval.n_Li6_to_He4_t_reaclib  +
          -rho*Y[jp]*Y[jli6]*rate_eval.p_Li6_to_He4_He3_reaclib  +
          -rho*Y[jd]*Y[jli6]*rate_eval.d_Li6_to_n_Be7_reaclib  +
          -rho*Y[jd]*Y[jli6]*rate_eval.d_Li6_to_p_Li7_reaclib  +
          +rho**2*Y[jn]*Y[jp]*Y[jhe4]*rate_eval.n_p_He4_to_Li6_reaclib
       )

    dYdt[jli7] = (
          +rho*ye(Y)*Y[jbe7]*rate_eval.Be7_to_Li7_reaclib  +
          +rho*Y[jt]*Y[jhe4]*rate_eval.He4_t_to_Li7_reaclib  +
          +rho*Y[jn]*Y[jli6]*rate_eval.n_Li6_to_Li7_reaclib  +
          +rho*Y[jd]*Y[jli6]*rate_eval.d_Li6_to_p_Li7_reaclib  +
          -rho*Y[jp]*Y[jli7]*rate_eval.p_Li7_to_He4_He4_reaclib  +
          +rho*Y[jn]*Y[jbe7]*rate_eval.n_Be7_to_p_Li7_reaclib  +
          -rho*Y[jd]*Y[jli7]*rate_eval.d_Li7_to_n_He4_He4_reaclib  +
          -rho*Y[jt]*Y[jli7]*rate_eval.t_Li7_to_n_n_He4_He4_reaclib  +
          -rho*Y[jhe3]*Y[jli7]*rate_eval.He3_Li7_to_n_p_He4_He4_reaclib
       )

    dYdt[jbe7] = (
          -rho*ye(Y)*Y[jbe7]*rate_eval.Be7_to_Li7_reaclib  +
          +rho*Y[jhe3]*Y[jhe4]*rate_eval.He4_He3_to_Be7_reaclib  +
          +rho*Y[jp]*Y[jli6]*rate_eval.p_Li6_to_Be7_reaclib  +
          +rho*Y[jd]*Y[jli6]*rate_eval.d_Li6_to_n_Be7_reaclib  +
          -rho*Y[jn]*Y[jbe7]*rate_eval.n_Be7_to_p_Li7_reaclib  +
          -rho*Y[jn]*Y[jbe7]*rate_eval.n_Be7_to_He4_He4_reaclib  +
          -rho*Y[jd]*Y[jbe7]*rate_eval.d_Be7_to_p_He4_He4_reaclib  +
          -rho*Y[jt]*Y[jbe7]*rate_eval.t_Be7_to_n_p_He4_He4_reaclib  +
          -rho*Y[jhe3]*Y[jbe7]*rate_eval.He3_Be7_to_p_p_He4_He4_reaclib
       )

    return dYdt

def jacobian(t, Y, rho, T, screen_func=None):
    return jacobian_eq(t, Y, rho, T, screen_func)

@numba.njit()
def jacobian_eq(t, Y, rho, T, screen_func):

    tf = Tfactors(T)
    rate_eval = RateEval()

    log_scor_d_t = 0.0
    log_scor_p_Li7 = 0.0
    log_scor_He3_He3 = 0.0
    log_scor_d_Li7 = 0.0
    log_scor_He3_He4 = 0.0
    log_scor_He3_Be7 = 0.0
    log_scor_p_p = 0.0
    log_scor_t_He3 = 0.0
    log_scor_p_d = 0.0
    log_scor_t_He4 = 0.0
    log_scor_p_He3 = 0.0
    log_scor_t_Be7 = 0.0
    log_scor_d_d = 0.0
    log_scor_He3_Li7 = 0.0
    log_scor_d_He3 = 0.0
    log_scor_p_He4 = 0.0
    log_scor_p_Li6 = 0.0
    log_scor_d_He4 = 0.0
    log_scor_d_Li6 = 0.0
    log_scor_d_Be7 = 0.0
    log_scor_t_t = 0.0
    log_scor_t_Li7 = 0.0
    log_scor_p_t = 0.0

    if screen_func is not None:
        plasma_state = PlasmaState(T, rho, Y, Z)

        scn_fac = ScreenFactors(1, 2, 1, 3)
        log_scor_d_t = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 1, 3, 7)
        log_scor_p_Li7 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(2, 3, 2, 3)
        log_scor_He3_He3 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 2, 3, 7)
        log_scor_d_Li7 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(2, 3, 2, 4)
        log_scor_He3_He4 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(2, 3, 4, 7)
        log_scor_He3_Be7 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 1, 1, 1)
        log_scor_p_p = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 3, 2, 3)
        log_scor_t_He3 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 1, 1, 2)
        log_scor_p_d = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 3, 2, 4)
        log_scor_t_He4 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 1, 2, 3)
        log_scor_p_He3 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 3, 4, 7)
        log_scor_t_Be7 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 2, 1, 2)
        log_scor_d_d = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(2, 3, 3, 7)
        log_scor_He3_Li7 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 2, 2, 3)
        log_scor_d_He3 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 1, 2, 4)
        log_scor_p_He4 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 1, 3, 6)
        log_scor_p_Li6 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 2, 2, 4)
        log_scor_d_He4 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 2, 3, 6)
        log_scor_d_Li6 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 2, 4, 7)
        log_scor_d_Be7 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 3, 1, 3)
        log_scor_t_t = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 3, 3, 7)
        log_scor_t_Li7 = screen_func(plasma_state, scn_fac)
        scn_fac = ScreenFactors(1, 1, 1, 3)
        log_scor_p_t = screen_func(plasma_state, scn_fac)

    # reaclib rates
    n_to_p_reaclib(rate_eval, tf)
    t_to_He3_reaclib(rate_eval, tf)
    He3_to_t_reaclib(rate_eval, tf)
    Be7_to_Li7_reaclib(rate_eval, tf)
    n_p_to_d_reaclib(rate_eval, tf)
    p_p_to_d_reaclib_beta_pos(rate_eval, tf, log_scor=log_scor_p_p)
    p_p_to_d_reaclib_electron_capture(rate_eval, tf, log_scor=log_scor_p_p)
    n_d_to_t_reaclib(rate_eval, tf)
    p_d_to_He3_reaclib(rate_eval, tf, log_scor=log_scor_p_d)
    d_d_to_He4_reaclib(rate_eval, tf, log_scor=log_scor_d_d)
    He4_d_to_Li6_reaclib(rate_eval, tf, log_scor=log_scor_d_He4)
    p_t_to_He4_reaclib(rate_eval, tf, log_scor=log_scor_p_t)
    He4_t_to_Li7_reaclib(rate_eval, tf, log_scor=log_scor_t_He4)
    n_He3_to_He4_reaclib(rate_eval, tf)
    p_He3_to_He4_reaclib(rate_eval, tf, log_scor=log_scor_p_He3)
    He4_He3_to_Be7_reaclib(rate_eval, tf, log_scor=log_scor_He3_He4)
    n_Li6_to_Li7_reaclib(rate_eval, tf)
    p_Li6_to_Be7_reaclib(rate_eval, tf, log_scor=log_scor_p_Li6)
    d_d_to_n_He3_reaclib(rate_eval, tf, log_scor=log_scor_d_d)
    d_d_to_p_t_reaclib(rate_eval, tf, log_scor=log_scor_d_d)
    d_t_to_n_He4_reaclib(rate_eval, tf, log_scor=log_scor_d_t)
    n_He3_to_p_t_reaclib(rate_eval, tf)
    d_He3_to_p_He4_reaclib(rate_eval, tf, log_scor=log_scor_d_He3)
    t_He3_to_d_He4_reaclib(rate_eval, tf, log_scor=log_scor_t_He3)
    n_Li6_to_He4_t_reaclib(rate_eval, tf)
    p_Li6_to_He4_He3_reaclib(rate_eval, tf, log_scor=log_scor_p_Li6)
    d_Li6_to_n_Be7_reaclib(rate_eval, tf, log_scor=log_scor_d_Li6)
    d_Li6_to_p_Li7_reaclib(rate_eval, tf, log_scor=log_scor_d_Li6)
    p_Li7_to_He4_He4_reaclib(rate_eval, tf, log_scor=log_scor_p_Li7)
    n_Be7_to_p_Li7_reaclib(rate_eval, tf)
    n_Be7_to_He4_He4_reaclib(rate_eval, tf)
    t_t_to_n_n_He4_reaclib(rate_eval, tf, log_scor=log_scor_t_t)
    t_He3_to_n_p_He4_reaclib(rate_eval, tf, log_scor=log_scor_t_He3)
    He3_He3_to_p_p_He4_reaclib(rate_eval, tf, log_scor=log_scor_He3_He3)
    d_Li7_to_n_He4_He4_reaclib(rate_eval, tf, log_scor=log_scor_d_Li7)
    d_Be7_to_p_He4_He4_reaclib(rate_eval, tf, log_scor=log_scor_d_Be7)
    t_Li7_to_n_n_He4_He4_reaclib(rate_eval, tf, log_scor=log_scor_t_Li7)
    He3_Li7_to_n_p_He4_He4_reaclib(rate_eval, tf, log_scor=log_scor_He3_Li7)
    t_Be7_to_n_p_He4_He4_reaclib(rate_eval, tf, log_scor=log_scor_t_Be7)
    He3_Be7_to_p_p_He4_He4_reaclib(rate_eval, tf, log_scor=log_scor_He3_Be7)
    n_p_He4_to_Li6_reaclib(rate_eval, tf, log_scor=log_scor_p_He4)
    n_p_p_to_p_d_reaclib(rate_eval, tf, log_scor=log_scor_p_p)

    jac = np.zeros((nnuc, nnuc), dtype=np.float64)

    jac[jn, jn] = (
       -rate_eval.n_to_p_reaclib
       -rho*Y[jp]*rate_eval.n_p_to_d_reaclib
       -rho*Y[jd]*rate_eval.n_d_to_t_reaclib
       -rho*Y[jhe3]*rate_eval.n_He3_to_He4_reaclib
       -rho*Y[jli6]*rate_eval.n_Li6_to_Li7_reaclib
       -rho*Y[jhe3]*rate_eval.n_He3_to_p_t_reaclib
       -rho*Y[jli6]*rate_eval.n_Li6_to_He4_t_reaclib
       -rho*Y[jbe7]*rate_eval.n_Be7_to_p_Li7_reaclib
       -rho*Y[jbe7]*rate_eval.n_Be7_to_He4_He4_reaclib
       -rho**2*Y[jp]*Y[jhe4]*rate_eval.n_p_He4_to_Li6_reaclib
       -5.00000000000000e-01*rho**2*Y[jp]**2*rate_eval.n_p_p_to_p_d_reaclib
       )

    jac[jn, jp] = (
       -rho*Y[jn]*rate_eval.n_p_to_d_reaclib
       -rho**2*Y[jn]*Y[jhe4]*rate_eval.n_p_He4_to_Li6_reaclib
       -5.00000000000000e-01*rho**2*Y[jn]*2*Y[jp]*rate_eval.n_p_p_to_p_d_reaclib
       )

    jac[jn, jd] = (
       -rho*Y[jn]*rate_eval.n_d_to_t_reaclib
       +5.00000000000000e-01*rho*2*Y[jd]*rate_eval.d_d_to_n_He3_reaclib
       +rho*Y[jt]*rate_eval.d_t_to_n_He4_reaclib
       +rho*Y[jli6]*rate_eval.d_Li6_to_n_Be7_reaclib
       +rho*Y[jli7]*rate_eval.d_Li7_to_n_He4_He4_reaclib
       )

    jac[jn, jt] = (
       +rho*Y[jd]*rate_eval.d_t_to_n_He4_reaclib
       +2*5.00000000000000e-01*rho*2*Y[jt]*rate_eval.t_t_to_n_n_He4_reaclib
       +rho*Y[jhe3]*rate_eval.t_He3_to_n_p_He4_reaclib
       +2*rho*Y[jli7]*rate_eval.t_Li7_to_n_n_He4_He4_reaclib
       +rho*Y[jbe7]*rate_eval.t_Be7_to_n_p_He4_He4_reaclib
       )

    jac[jn, jhe3] = (
       -rho*Y[jn]*rate_eval.n_He3_to_He4_reaclib
       -rho*Y[jn]*rate_eval.n_He3_to_p_t_reaclib
       +rho*Y[jt]*rate_eval.t_He3_to_n_p_He4_reaclib
       +rho*Y[jli7]*rate_eval.He3_Li7_to_n_p_He4_He4_reaclib
       )

    jac[jn, jhe4] = (
       -rho**2*Y[jn]*Y[jp]*rate_eval.n_p_He4_to_Li6_reaclib
       )

    jac[jn, jli6] = (
       -rho*Y[jn]*rate_eval.n_Li6_to_Li7_reaclib
       -rho*Y[jn]*rate_eval.n_Li6_to_He4_t_reaclib
       +rho*Y[jd]*rate_eval.d_Li6_to_n_Be7_reaclib
       )

    jac[jn, jli7] = (
       +rho*Y[jd]*rate_eval.d_Li7_to_n_He4_He4_reaclib
       +2*rho*Y[jt]*rate_eval.t_Li7_to_n_n_He4_He4_reaclib
       +rho*Y[jhe3]*rate_eval.He3_Li7_to_n_p_He4_He4_reaclib
       )

    jac[jn, jbe7] = (
       -rho*Y[jn]*rate_eval.n_Be7_to_p_Li7_reaclib
       -rho*Y[jn]*rate_eval.n_Be7_to_He4_He4_reaclib
       +rho*Y[jt]*rate_eval.t_Be7_to_n_p_He4_He4_reaclib
       )

    jac[jp, jn] = (
       -rho*Y[jp]*rate_eval.n_p_to_d_reaclib
       -rho**2*Y[jp]*Y[jhe4]*rate_eval.n_p_He4_to_Li6_reaclib
       -2*5.00000000000000e-01*rho**2*Y[jp]**2*rate_eval.n_p_p_to_p_d_reaclib
       +rate_eval.n_to_p_reaclib
       +rho*Y[jhe3]*rate_eval.n_He3_to_p_t_reaclib
       +rho*Y[jbe7]*rate_eval.n_Be7_to_p_Li7_reaclib
       +5.00000000000000e-01*rho**2*Y[jp]**2*rate_eval.n_p_p_to_p_d_reaclib
       )

    jac[jp, jp] = (
       -rho*Y[jn]*rate_eval.n_p_to_d_reaclib
       -2*5.00000000000000e-01*rho*2*Y[jp]*rate_eval.p_p_to_d_reaclib_beta_pos
       -2*5.00000000000000e-01*rho**2*ye(Y)*2*Y[jp]*rate_eval.p_p_to_d_reaclib_electron_capture
       -rho*Y[jd]*rate_eval.p_d_to_He3_reaclib
       -rho*Y[jt]*rate_eval.p_t_to_He4_reaclib
       -rho*Y[jhe3]*rate_eval.p_He3_to_He4_reaclib
       -rho*Y[jli6]*rate_eval.p_Li6_to_Be7_reaclib
       -rho*Y[jli6]*rate_eval.p_Li6_to_He4_He3_reaclib
       -rho*Y[jli7]*rate_eval.p_Li7_to_He4_He4_reaclib
       -rho**2*Y[jn]*Y[jhe4]*rate_eval.n_p_He4_to_Li6_reaclib
       -2*5.00000000000000e-01*rho**2*Y[jn]*2*Y[jp]*rate_eval.n_p_p_to_p_d_reaclib
       +5.00000000000000e-01*rho**2*Y[jn]*2*Y[jp]*rate_eval.n_p_p_to_p_d_reaclib
       )

    jac[jp, jd] = (
       -rho*Y[jp]*rate_eval.p_d_to_He3_reaclib
       +5.00000000000000e-01*rho*2*Y[jd]*rate_eval.d_d_to_p_t_reaclib
       +rho*Y[jhe3]*rate_eval.d_He3_to_p_He4_reaclib
       +rho*Y[jli6]*rate_eval.d_Li6_to_p_Li7_reaclib
       +rho*Y[jbe7]*rate_eval.d_Be7_to_p_He4_He4_reaclib
       )

    jac[jp, jt] = (
       -rho*Y[jp]*rate_eval.p_t_to_He4_reaclib
       +rho*Y[jhe3]*rate_eval.t_He3_to_n_p_He4_reaclib
       +rho*Y[jbe7]*rate_eval.t_Be7_to_n_p_He4_He4_reaclib
       )

    jac[jp, jhe3] = (
       -rho*Y[jp]*rate_eval.p_He3_to_He4_reaclib
       +rho*Y[jn]*rate_eval.n_He3_to_p_t_reaclib
       +rho*Y[jd]*rate_eval.d_He3_to_p_He4_reaclib
       +rho*Y[jt]*rate_eval.t_He3_to_n_p_He4_reaclib
       +2*5.00000000000000e-01*rho*2*Y[jhe3]*rate_eval.He3_He3_to_p_p_He4_reaclib
       +rho*Y[jli7]*rate_eval.He3_Li7_to_n_p_He4_He4_reaclib
       +2*rho*Y[jbe7]*rate_eval.He3_Be7_to_p_p_He4_He4_reaclib
       )

    jac[jp, jhe4] = (
       -rho**2*Y[jn]*Y[jp]*rate_eval.n_p_He4_to_Li6_reaclib
       )

    jac[jp, jli6] = (
       -rho*Y[jp]*rate_eval.p_Li6_to_Be7_reaclib
       -rho*Y[jp]*rate_eval.p_Li6_to_He4_He3_reaclib
       +rho*Y[jd]*rate_eval.d_Li6_to_p_Li7_reaclib
       )

    jac[jp, jli7] = (
       -rho*Y[jp]*rate_eval.p_Li7_to_He4_He4_reaclib
       +rho*Y[jhe3]*rate_eval.He3_Li7_to_n_p_He4_He4_reaclib
       )

    jac[jp, jbe7] = (
       +rho*Y[jn]*rate_eval.n_Be7_to_p_Li7_reaclib
       +rho*Y[jd]*rate_eval.d_Be7_to_p_He4_He4_reaclib
       +rho*Y[jt]*rate_eval.t_Be7_to_n_p_He4_He4_reaclib
       +2*rho*Y[jhe3]*rate_eval.He3_Be7_to_p_p_He4_He4_reaclib
       )

    jac[jd, jn] = (
       -rho*Y[jd]*rate_eval.n_d_to_t_reaclib
       +rho*Y[jp]*rate_eval.n_p_to_d_reaclib
       +5.00000000000000e-01*rho**2*Y[jp]**2*rate_eval.n_p_p_to_p_d_reaclib
       )

    jac[jd, jp] = (
       -rho*Y[jd]*rate_eval.p_d_to_He3_reaclib
       +rho*Y[jn]*rate_eval.n_p_to_d_reaclib
       +5.00000000000000e-01*rho*2*Y[jp]*rate_eval.p_p_to_d_reaclib_beta_pos
       +5.00000000000000e-01*rho**2*ye(Y)*2*Y[jp]*rate_eval.p_p_to_d_reaclib_electron_capture
       +5.00000000000000e-01*rho**2*Y[jn]*2*Y[jp]*rate_eval.n_p_p_to_p_d_reaclib
       )

    jac[jd, jd] = (
       -rho*Y[jn]*rate_eval.n_d_to_t_reaclib
       -rho*Y[jp]*rate_eval.p_d_to_He3_reaclib
       -2*5.00000000000000e-01*rho*2*Y[jd]*rate_eval.d_d_to_He4_reaclib
       -rho*Y[jhe4]*rate_eval.He4_d_to_Li6_reaclib
       -2*5.00000000000000e-01*rho*2*Y[jd]*rate_eval.d_d_to_n_He3_reaclib
       -2*5.00000000000000e-01*rho*2*Y[jd]*rate_eval.d_d_to_p_t_reaclib
       -rho*Y[jt]*rate_eval.d_t_to_n_He4_reaclib
       -rho*Y[jhe3]*rate_eval.d_He3_to_p_He4_reaclib
       -rho*Y[jli6]*rate_eval.d_Li6_to_n_Be7_reaclib
       -rho*Y[jli6]*rate_eval.d_Li6_to_p_Li7_reaclib
       -rho*Y[jli7]*rate_eval.d_Li7_to_n_He4_He4_reaclib
       -rho*Y[jbe7]*rate_eval.d_Be7_to_p_He4_He4_reaclib
       )

    jac[jd, jt] = (
       -rho*Y[jd]*rate_eval.d_t_to_n_He4_reaclib
       +rho*Y[jhe3]*rate_eval.t_He3_to_d_He4_reaclib
       )

    jac[jd, jhe3] = (
       -rho*Y[jd]*rate_eval.d_He3_to_p_He4_reaclib
       +rho*Y[jt]*rate_eval.t_He3_to_d_He4_reaclib
       )

    jac[jd, jhe4] = (
       -rho*Y[jd]*rate_eval.He4_d_to_Li6_reaclib
       )

    jac[jd, jli6] = (
       -rho*Y[jd]*rate_eval.d_Li6_to_n_Be7_reaclib
       -rho*Y[jd]*rate_eval.d_Li6_to_p_Li7_reaclib
       )

    jac[jd, jli7] = (
       -rho*Y[jd]*rate_eval.d_Li7_to_n_He4_He4_reaclib
       )

    jac[jd, jbe7] = (
       -rho*Y[jd]*rate_eval.d_Be7_to_p_He4_He4_reaclib
       )

    jac[jt, jn] = (
       +rho*Y[jd]*rate_eval.n_d_to_t_reaclib
       +rho*Y[jhe3]*rate_eval.n_He3_to_p_t_reaclib
       +rho*Y[jli6]*rate_eval.n_Li6_to_He4_t_reaclib
       )

    jac[jt, jp] = (
       -rho*Y[jt]*rate_eval.p_t_to_He4_reaclib
       )

    jac[jt, jd] = (
       -rho*Y[jt]*rate_eval.d_t_to_n_He4_reaclib
       +rho*Y[jn]*rate_eval.n_d_to_t_reaclib
       +5.00000000000000e-01*rho*2*Y[jd]*rate_eval.d_d_to_p_t_reaclib
       )

    jac[jt, jt] = (
       -rate_eval.t_to_He3_reaclib
       -rho*Y[jp]*rate_eval.p_t_to_He4_reaclib
       -rho*Y[jhe4]*rate_eval.He4_t_to_Li7_reaclib
       -rho*Y[jd]*rate_eval.d_t_to_n_He4_reaclib
       -rho*Y[jhe3]*rate_eval.t_He3_to_d_He4_reaclib
       -2*5.00000000000000e-01*rho*2*Y[jt]*rate_eval.t_t_to_n_n_He4_reaclib
       -rho*Y[jhe3]*rate_eval.t_He3_to_n_p_He4_reaclib
       -rho*Y[jli7]*rate_eval.t_Li7_to_n_n_He4_He4_reaclib
       -rho*Y[jbe7]*rate_eval.t_Be7_to_n_p_He4_He4_reaclib
       )

    jac[jt, jhe3] = (
       -rho*Y[jt]*rate_eval.t_He3_to_d_He4_reaclib
       -rho*Y[jt]*rate_eval.t_He3_to_n_p_He4_reaclib
       +rho*ye(Y)*rate_eval.He3_to_t_reaclib
       +rho*Y[jn]*rate_eval.n_He3_to_p_t_reaclib
       )

    jac[jt, jhe4] = (
       -rho*Y[jt]*rate_eval.He4_t_to_Li7_reaclib
       )

    jac[jt, jli6] = (
       +rho*Y[jn]*rate_eval.n_Li6_to_He4_t_reaclib
       )

    jac[jt, jli7] = (
       -rho*Y[jt]*rate_eval.t_Li7_to_n_n_He4_He4_reaclib
       )

    jac[jt, jbe7] = (
       -rho*Y[jt]*rate_eval.t_Be7_to_n_p_He4_He4_reaclib
       )

    jac[jhe3, jn] = (
       -rho*Y[jhe3]*rate_eval.n_He3_to_He4_reaclib
       -rho*Y[jhe3]*rate_eval.n_He3_to_p_t_reaclib
       )

    jac[jhe3, jp] = (
       -rho*Y[jhe3]*rate_eval.p_He3_to_He4_reaclib
       +rho*Y[jd]*rate_eval.p_d_to_He3_reaclib
       +rho*Y[jli6]*rate_eval.p_Li6_to_He4_He3_reaclib
       )

    jac[jhe3, jd] = (
       -rho*Y[jhe3]*rate_eval.d_He3_to_p_He4_reaclib
       +rho*Y[jp]*rate_eval.p_d_to_He3_reaclib
       +5.00000000000000e-01*rho*2*Y[jd]*rate_eval.d_d_to_n_He3_reaclib
       )

    jac[jhe3, jt] = (
       -rho*Y[jhe3]*rate_eval.t_He3_to_d_He4_reaclib
       -rho*Y[jhe3]*rate_eval.t_He3_to_n_p_He4_reaclib
       +rate_eval.t_to_He3_reaclib
       )

    jac[jhe3, jhe3] = (
       -rho*ye(Y)*rate_eval.He3_to_t_reaclib
       -rho*Y[jn]*rate_eval.n_He3_to_He4_reaclib
       -rho*Y[jp]*rate_eval.p_He3_to_He4_reaclib
       -rho*Y[jhe4]*rate_eval.He4_He3_to_Be7_reaclib
       -rho*Y[jn]*rate_eval.n_He3_to_p_t_reaclib
       -rho*Y[jd]*rate_eval.d_He3_to_p_He4_reaclib
       -rho*Y[jt]*rate_eval.t_He3_to_d_He4_reaclib
       -rho*Y[jt]*rate_eval.t_He3_to_n_p_He4_reaclib
       -2*5.00000000000000e-01*rho*2*Y[jhe3]*rate_eval.He3_He3_to_p_p_He4_reaclib
       -rho*Y[jli7]*rate_eval.He3_Li7_to_n_p_He4_He4_reaclib
       -rho*Y[jbe7]*rate_eval.He3_Be7_to_p_p_He4_He4_reaclib
       )

    jac[jhe3, jhe4] = (
       -rho*Y[jhe3]*rate_eval.He4_He3_to_Be7_reaclib
       )

    jac[jhe3, jli6] = (
       +rho*Y[jp]*rate_eval.p_Li6_to_He4_He3_reaclib
       )

    jac[jhe3, jli7] = (
       -rho*Y[jhe3]*rate_eval.He3_Li7_to_n_p_He4_He4_reaclib
       )

    jac[jhe3, jbe7] = (
       -rho*Y[jhe3]*rate_eval.He3_Be7_to_p_p_He4_He4_reaclib
       )

    jac[jhe4, jn] = (
       -rho**2*Y[jp]*Y[jhe4]*rate_eval.n_p_He4_to_Li6_reaclib
       +rho*Y[jhe3]*rate_eval.n_He3_to_He4_reaclib
       +rho*Y[jli6]*rate_eval.n_Li6_to_He4_t_reaclib
       +2*rho*Y[jbe7]*rate_eval.n_Be7_to_He4_He4_reaclib
       )

    jac[jhe4, jp] = (
       -rho**2*Y[jn]*Y[jhe4]*rate_eval.n_p_He4_to_Li6_reaclib
       +rho*Y[jt]*rate_eval.p_t_to_He4_reaclib
       +rho*Y[jhe3]*rate_eval.p_He3_to_He4_reaclib
       +rho*Y[jli6]*rate_eval.p_Li6_to_He4_He3_reaclib
       +2*rho*Y[jli7]*rate_eval.p_Li7_to_He4_He4_reaclib
       )

    jac[jhe4, jd] = (
       -rho*Y[jhe4]*rate_eval.He4_d_to_Li6_reaclib
       +5.00000000000000e-01*rho*2*Y[jd]*rate_eval.d_d_to_He4_reaclib
       +rho*Y[jt]*rate_eval.d_t_to_n_He4_reaclib
       +rho*Y[jhe3]*rate_eval.d_He3_to_p_He4_reaclib
       +2*rho*Y[jli7]*rate_eval.d_Li7_to_n_He4_He4_reaclib
       +2*rho*Y[jbe7]*rate_eval.d_Be7_to_p_He4_He4_reaclib
       )

    jac[jhe4, jt] = (
       -rho*Y[jhe4]*rate_eval.He4_t_to_Li7_reaclib
       +rho*Y[jp]*rate_eval.p_t_to_He4_reaclib
       +rho*Y[jd]*rate_eval.d_t_to_n_He4_reaclib
       +rho*Y[jhe3]*rate_eval.t_He3_to_d_He4_reaclib
       +5.00000000000000e-01*rho*2*Y[jt]*rate_eval.t_t_to_n_n_He4_reaclib
       +rho*Y[jhe3]*rate_eval.t_He3_to_n_p_He4_reaclib
       +2*rho*Y[jli7]*rate_eval.t_Li7_to_n_n_He4_He4_reaclib
       +2*rho*Y[jbe7]*rate_eval.t_Be7_to_n_p_He4_He4_reaclib
       )

    jac[jhe4, jhe3] = (
       -rho*Y[jhe4]*rate_eval.He4_He3_to_Be7_reaclib
       +rho*Y[jn]*rate_eval.n_He3_to_He4_reaclib
       +rho*Y[jp]*rate_eval.p_He3_to_He4_reaclib
       +rho*Y[jd]*rate_eval.d_He3_to_p_He4_reaclib
       +rho*Y[jt]*rate_eval.t_He3_to_d_He4_reaclib
       +rho*Y[jt]*rate_eval.t_He3_to_n_p_He4_reaclib
       +5.00000000000000e-01*rho*2*Y[jhe3]*rate_eval.He3_He3_to_p_p_He4_reaclib
       +2*rho*Y[jli7]*rate_eval.He3_Li7_to_n_p_He4_He4_reaclib
       +2*rho*Y[jbe7]*rate_eval.He3_Be7_to_p_p_He4_He4_reaclib
       )

    jac[jhe4, jhe4] = (
       -rho*Y[jd]*rate_eval.He4_d_to_Li6_reaclib
       -rho*Y[jt]*rate_eval.He4_t_to_Li7_reaclib
       -rho*Y[jhe3]*rate_eval.He4_He3_to_Be7_reaclib
       -rho**2*Y[jn]*Y[jp]*rate_eval.n_p_He4_to_Li6_reaclib
       )

    jac[jhe4, jli6] = (
       +rho*Y[jn]*rate_eval.n_Li6_to_He4_t_reaclib
       +rho*Y[jp]*rate_eval.p_Li6_to_He4_He3_reaclib
       )

    jac[jhe4, jli7] = (
       +2*rho*Y[jp]*rate_eval.p_Li7_to_He4_He4_reaclib
       +2*rho*Y[jd]*rate_eval.d_Li7_to_n_He4_He4_reaclib
       +2*rho*Y[jt]*rate_eval.t_Li7_to_n_n_He4_He4_reaclib
       +2*rho*Y[jhe3]*rate_eval.He3_Li7_to_n_p_He4_He4_reaclib
       )

    jac[jhe4, jbe7] = (
       +2*rho*Y[jn]*rate_eval.n_Be7_to_He4_He4_reaclib
       +2*rho*Y[jd]*rate_eval.d_Be7_to_p_He4_He4_reaclib
       +2*rho*Y[jt]*rate_eval.t_Be7_to_n_p_He4_He4_reaclib
       +2*rho*Y[jhe3]*rate_eval.He3_Be7_to_p_p_He4_He4_reaclib
       )

    jac[jli6, jn] = (
       -rho*Y[jli6]*rate_eval.n_Li6_to_Li7_reaclib
       -rho*Y[jli6]*rate_eval.n_Li6_to_He4_t_reaclib
       +rho**2*Y[jp]*Y[jhe4]*rate_eval.n_p_He4_to_Li6_reaclib
       )

    jac[jli6, jp] = (
       -rho*Y[jli6]*rate_eval.p_Li6_to_Be7_reaclib
       -rho*Y[jli6]*rate_eval.p_Li6_to_He4_He3_reaclib
       +rho**2*Y[jn]*Y[jhe4]*rate_eval.n_p_He4_to_Li6_reaclib
       )

    jac[jli6, jd] = (
       -rho*Y[jli6]*rate_eval.d_Li6_to_n_Be7_reaclib
       -rho*Y[jli6]*rate_eval.d_Li6_to_p_Li7_reaclib
       +rho*Y[jhe4]*rate_eval.He4_d_to_Li6_reaclib
       )

    jac[jli6, jhe4] = (
       +rho*Y[jd]*rate_eval.He4_d_to_Li6_reaclib
       +rho**2*Y[jn]*Y[jp]*rate_eval.n_p_He4_to_Li6_reaclib
       )

    jac[jli6, jli6] = (
       -rho*Y[jn]*rate_eval.n_Li6_to_Li7_reaclib
       -rho*Y[jp]*rate_eval.p_Li6_to_Be7_reaclib
       -rho*Y[jn]*rate_eval.n_Li6_to_He4_t_reaclib
       -rho*Y[jp]*rate_eval.p_Li6_to_He4_He3_reaclib
       -rho*Y[jd]*rate_eval.d_Li6_to_n_Be7_reaclib
       -rho*Y[jd]*rate_eval.d_Li6_to_p_Li7_reaclib
       )

    jac[jli7, jn] = (
       +rho*Y[jli6]*rate_eval.n_Li6_to_Li7_reaclib
       +rho*Y[jbe7]*rate_eval.n_Be7_to_p_Li7_reaclib
       )

    jac[jli7, jp] = (
       -rho*Y[jli7]*rate_eval.p_Li7_to_He4_He4_reaclib
       )

    jac[jli7, jd] = (
       -rho*Y[jli7]*rate_eval.d_Li7_to_n_He4_He4_reaclib
       +rho*Y[jli6]*rate_eval.d_Li6_to_p_Li7_reaclib
       )

    jac[jli7, jt] = (
       -rho*Y[jli7]*rate_eval.t_Li7_to_n_n_He4_He4_reaclib
       +rho*Y[jhe4]*rate_eval.He4_t_to_Li7_reaclib
       )

    jac[jli7, jhe3] = (
       -rho*Y[jli7]*rate_eval.He3_Li7_to_n_p_He4_He4_reaclib
       )

    jac[jli7, jhe4] = (
       +rho*Y[jt]*rate_eval.He4_t_to_Li7_reaclib
       )

    jac[jli7, jli6] = (
       +rho*Y[jn]*rate_eval.n_Li6_to_Li7_reaclib
       +rho*Y[jd]*rate_eval.d_Li6_to_p_Li7_reaclib
       )

    jac[jli7, jli7] = (
       -rho*Y[jp]*rate_eval.p_Li7_to_He4_He4_reaclib
       -rho*Y[jd]*rate_eval.d_Li7_to_n_He4_He4_reaclib
       -rho*Y[jt]*rate_eval.t_Li7_to_n_n_He4_He4_reaclib
       -rho*Y[jhe3]*rate_eval.He3_Li7_to_n_p_He4_He4_reaclib
       )

    jac[jli7, jbe7] = (
       +rho*ye(Y)*rate_eval.Be7_to_Li7_reaclib
       +rho*Y[jn]*rate_eval.n_Be7_to_p_Li7_reaclib
       )

    jac[jbe7, jn] = (
       -rho*Y[jbe7]*rate_eval.n_Be7_to_p_Li7_reaclib
       -rho*Y[jbe7]*rate_eval.n_Be7_to_He4_He4_reaclib
       )

    jac[jbe7, jp] = (
       +rho*Y[jli6]*rate_eval.p_Li6_to_Be7_reaclib
       )

    jac[jbe7, jd] = (
       -rho*Y[jbe7]*rate_eval.d_Be7_to_p_He4_He4_reaclib
       +rho*Y[jli6]*rate_eval.d_Li6_to_n_Be7_reaclib
       )

    jac[jbe7, jt] = (
       -rho*Y[jbe7]*rate_eval.t_Be7_to_n_p_He4_He4_reaclib
       )

    jac[jbe7, jhe3] = (
       -rho*Y[jbe7]*rate_eval.He3_Be7_to_p_p_He4_He4_reaclib
       +rho*Y[jhe4]*rate_eval.He4_He3_to_Be7_reaclib
       )

    jac[jbe7, jhe4] = (
       +rho*Y[jhe3]*rate_eval.He4_He3_to_Be7_reaclib
       )

    jac[jbe7, jli6] = (
       +rho*Y[jp]*rate_eval.p_Li6_to_Be7_reaclib
       +rho*Y[jd]*rate_eval.d_Li6_to_n_Be7_reaclib
       )

    jac[jbe7, jbe7] = (
       -rho*ye(Y)*rate_eval.Be7_to_Li7_reaclib
       -rho*Y[jn]*rate_eval.n_Be7_to_p_Li7_reaclib
       -rho*Y[jn]*rate_eval.n_Be7_to_He4_He4_reaclib
       -rho*Y[jd]*rate_eval.d_Be7_to_p_He4_He4_reaclib
       -rho*Y[jt]*rate_eval.t_Be7_to_n_p_He4_He4_reaclib
       -rho*Y[jhe3]*rate_eval.He3_Be7_to_p_p_He4_He4_reaclib
       )

    return jac
