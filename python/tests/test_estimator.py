from __future__ import division
from unittest import TestCase
import sys
import pywt

from wde.estimator import *

def ints(xs):
    return tuple([int(x) for x in xs])

def mkwave(short_name):
    wave = pywt.Wavelet(short_name)
    phi, psi, _ = wave.wavefun(level=8)
    phi_supp, psi_supp = wave_support_info(wave)
    return wave, phi, psi, phi_supp, psi_supp


class AllTests(TestCase):
    def test_all_qx(self):
        dim = 3
        vals = list(all_qx(dim))
        strs = [''.join([str(i) for i in ix2s]) for wx, ix2s in vals]
        self.assertEqual(['000','001','010','011','100','101','110','111'], strs)
        self.assertEqual(range(dim**2-1),[v for v,ix2 in vals])

    def test_gridify_xs(self):
        p1 = [0.25,0.25]
        p2 = [0.25,0.5]
        p3 = [0.75,0.333]
        p4 = [0.75,0.75]
        xs = np.array([p1,p2,p3,p4])
        minx = np.amin(xs,axis=1)
        maxx = np.amax(xs,axis=1)
        grid = gridify_xs(0,2,xs,minx,maxx)
        self.assertTrue((0,0) in grid[0])
        self.assertEquals(4, len(grid[0][(0,0)][0]))
        for zs, where_xs in grid[1].items():
            self.assertEquals(1, len(where_xs[0]))
        for zs, where_xs in grid[2].items():
            for x in xs[where_xs]:
                self.assertEqual(zs, ints(4 * x))

    def test_calc_factor(self):
        factor = calc_factor(100, 3, 1)
        self.assertAlmostEqual(0.23094, factor, 6)
        factor = calc_factor(121, 2, 3)
        self.assertAlmostEqual(0.0969697, factor, 6)

    def test_calculate_nearest_balls(self):
        p1 = [0.25,0.25]
        p2 = [0.25,0.5]
        p3 = [0.75,0.333]
        p4 = [0.75,0.75]
        xs = np.array([p1,p2,p3,p4])
        xs_1 = calculate_nearest_balls(1, xs)
        self.assertAlmostEqualTuple((0.25, 0.25, 0.417, 0.417), tuple(xs_1[:,2].tolist()), 3)
        xs_2 = calculate_nearest_balls(2, xs)
        self.assertAlmostEqualTuple((0.33789479, 0.35143452, 0.33789479, 0.372678), tuple(xs_2[:,2].tolist()), 3)

    def test_calc_coeff(self):
        p1 = [0.25,0.25]
        p2 = [0.25,0.5]
        p3 = [0.75,0.333]
        p4 = [0.75,0.75]
        xs = np.array([p1,p2,p3,p4])
        minx = np.amin(xs,axis=1)
        maxx = np.amax(xs,axis=1)
        wave, phi, psi, phi_supp, psi_supp = mkwave('db1')
        phi = interp1d(np.linspace(*phi_supp, num=len(phi)), phi, fill_value=0.0, bounds_error=False)
        psi = interp1d(np.linspace(*psi_supp, num=len(psi)), psi, fill_value=0.0, bounds_error=False)
        grid = gridify_xs(0,1,xs,minx,maxx)
        xs_1 = calculate_nearest_balls(1, xs)
        xs_2 = calculate_nearest_balls(2, xs)
        # j=0
        j = 0
        jpow = 2 ** j
        # ix2s=(0,0)
        wavef = partial(wave_tensor, (0,0), phi, psi)
        # zs = (0,0) for k=1,2
        v = calc_coeff(wavef, jpow, (0,0), grid[j], xs_1)
        self.assertAlmostEqual(1.3340, v, 3) # TODO verify !!!
        v = calc_coeff(wavef, jpow, (0,0), grid[j], xs_2)
        self.assertAlmostEqual(1.400, v, 3) # TODO verify !!!
        # ix2s=(1,0)
        wavef = partial(wave_tensor, (1,0), phi, psi)
        # zs = (1,0) for k=1,2
        xs_1 = calculate_nearest_balls(1, xs)
        xs_2 = calculate_nearest_balls(2, xs)
        v = calc_coeff(wavef, jpow, (1,0), grid[j], xs_1)
        self.assertAlmostEqual(1.3340, v, 3) # TODO verify !!!
        v = calc_coeff(wavef, jpow, (0,1), grid[j], xs_2)
        self.assertAlmostEqual(1.400, v, 3) # TODO verify !!!
        # j = 1
        j = 1
        jpow = 2 ** j
        # ix2s=(0,1)
        wavef = partial(wave_tensor, (0,0), phi, psi)
        # zs = (0,0) for k=1,2
        xs_1 = calculate_nearest_balls(1, xs)
        xs_2 = calculate_nearest_balls(2, xs)
        v = calc_coeff(wavef, jpow, (0,0), grid[j], xs_1)
        self.assertAlmostEqual(1.3340, v, 3) # TODO verify !!!
        v = calc_coeff(wavef, jpow, (1,1), grid[j], xs_2)
        self.assertAlmostEqual(1.400, v, 3) # TODO verify !!!
        j = 0
        jpow = 2 ** j
        xs_1 = calculate_nearest_balls(1, xs)
        xs_2 = calculate_nearest_balls(2, xs)
        wavef = partial(wave_tensor, (0,0), phi, psi)
        v = calc_coeff(wavef, jpow, (0,0), grid[j], xs_1)
        self.assertAlmostEqual(1.3340, v, 3) # TODO verify !!!
        v = calc_coeff(wavef, jpow, (0,0), grid[j], xs_2)
        self.assertAlmostEqual(1.400, v, 3) # TODO verify !!!

    def assertAlmostEqualTuple(self, t1, t2, num):
        for v1,v2 in zip(t1,t2):
            self.assertAlmostEqual(v1,v2,num)
