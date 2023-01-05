import math
import numpy as np

class Statics:
    """Class to hold the static variables between depth and time steps"""
    cumi = 0
    i_end = 0
    lamb = 0
    tp = 0
    pond = 0


def infiltration(time, dt, ppt, k0, cd, m, statics):
    """This format is a python adaption from original Bevin fortran 77 code.
    Variables match the formatting found in that code."""
    f1 = 0.0
    t = time
    if ppt <= 0:
        # no rain.
        statics.cumi = 0
        statics.lamb = 0
        statics.tp = 0
        statics.i_end = 0
        statics.pond = 0

        return 0.0

    if statics.pond == 0:
        if statics.cumi > 0:
            f1 = statics.cumi
            nf = -k0 * m * (cd + f1) / (1 - np.exp(f1 * m))
            if nf < ppt:
                statics.i_end = statics.cumi
                statics.tp = t - dt
                statics.pond = 1
                statics.lamb = 0
        f2 = statics.cumi + ppt * dt
        nf = (-k0 * m * (cd + f2)) / (1 - np.exp(f2 * m))
        if f2 == 0.0 or nf > ppt:
            didt = ppt
            statics.cumi = statics.cumi + didt * dt
            statics.ponding = 0
            return didt

        statics.i_end = statics.cumi + nf * dt
        for i in range(0, 21):
            nf = -k0 * m * (cd + statics.i_end) / (1 - np.exp(statics.i_end * m))
            if nf > ppt:
                f1 = statics.i_end
                statics.i_end = (statics.i_end + f2) / 2.0
                df = statics.i_end - f1
            else:
                f2 = statics.i_end
                statics.i_end = (statics.i_end + f1) / 2.0
                df = statics.i_end - f2
            if abs(df) <= 0.00001:
                break
            if i == 20:
                print("Warning: max iter exceeded at {}".format(t))

        statics.tp = t - dt + (statics.i_end - statics.cumi) / ppt
        if statics.tp > t:
            didt = ppt
            statics.cumi = statics.cumi + didt
            statics.ponding = 0
            return didt

        statics.pond = 1

    if statics.lamb == 0:
        fact = 1
        icd = statics.i_end + cd
        for j in range(1, 11):
            fact = fact * j
            add = (icd * m) ** j / (j * fact)
            statics.lamb = statics.lamb + add
        statics.lamb = math.log(icd) - (math.log(icd) + statics.lamb) / np.exp(cd * m)

        statics.ponding = 1

    statics.i_end = statics.i_end + ppt * (t - statics.tp) / 2.0
    for i in range(0, 21):
        icd = statics.i_end + cd
        add = 0
        fact = 1
        for j in range(1, 11):
            fact = fact * j
            add_0 = math.pow((icd * m), j) / (j * fact)
            add = add + add_0
        f1 = -(math.log(icd) - (math.log(icd) + add) / np.exp(cd * m) - statics.lamb) / (k0 * m) - (
                    t - statics.tp)
        f2 = (np.exp(statics.i_end * m) - 1) / (icd * k0 * m)
        df = -f1 / f2
        statics.i_end = statics.i_end + df
        if abs(df) <= 0.000001:
            break

    if statics.i_end < statics.cumi + ppt:
        didt = (statics.i_end - statics.cumi)
        statics.cumi = statics.i_end
        statics.i_end = statics.i_end / dt
        statics.pond = 1

    else:
        didt = ppt * dt
        statics.cumi = statics.cumi + didt
        statics.pond = 0
    return didt

def static_reset(statics):
    """In the case of no infiltration, statics are reset to 0"""
    statics.cumi = 0
    statics.lamb = 0
    statics.tp = 0
    statics.i_end = 0
    statics.pond = 0