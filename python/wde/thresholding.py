from __future__ import division
import math

def soft_threshold_half(threshold):
    def f(n, j, dn, coeff):
        lvl_factor = math.sqrt(math.sqrt(j + 1) / math.sqrt(n))
        lvl_t = threshold * lvl_factor
        if coeff < 0:
            if -coeff < lvl_t:
                return 0
            else:
                return coeff + lvl_t
        else:
            if coeff < lvl_t:
                return 0
            else:
                return coeff - lvl_t
    return f

def soft_threshold_2(threshold):
    def f(n, j, dn, coeff):
        lvl_factor = (j+1)/n
        lvl_t = threshold * lvl_factor
        if coeff < 0:
            if -coeff < lvl_t:
                return 0
            else:
                return coeff + lvl_t
        else:
            if coeff < lvl_t:
                return 0
            else:
                return coeff - lvl_t
    return f

def soft_threshold(threshold):
    def f(n, j, dn, coeff):
        lvl_factor = math.sqrt(j + 1) / math.sqrt(n)
        lvl_t = threshold * lvl_factor
        if coeff < 0:
            if -coeff < lvl_t:
                return 0
            else:
                return coeff + lvl_t
        else:
            if coeff < lvl_t:
                return 0
            else:
                return coeff - lvl_t
    return f

def hard_threshold(threshold):
    def f(n, j, dn, coeff):
        lvl_factor = math.sqrt(j + 1) / math.sqrt(n)
        lvl_t = threshold * lvl_factor
        if coeff < 0:
            if -coeff < lvl_t:
                return 0
            else:
                return coeff
        else:
            if coeff < lvl_t:
                return 0
            else:
                return coeff
    return f

def block_threshold(dn_threshold):
    def f(n, j, dn, coeff):
        if dn < dn_threshold:
            return 0
        else:
            return coeff
    return f

def soft_block_threshold(th1, th2):
    def f(n, j, dn, coeff):
        if dn == 0:
            return 0
        lvl_t = th1 * math.sqrt(j + 1) / math.sqrt(dn) + th2
        if coeff < 0:
            if -coeff < lvl_t:
                return 0
            else:
                return coeff + lvl_t
        else:
            if coeff < lvl_t:
                return 0
            else:
                return coeff - lvl_t
    return f
