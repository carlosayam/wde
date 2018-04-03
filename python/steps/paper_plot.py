"""
Utility library to plot and generate graphs for paper
"""
from __future__ import division
import math
import sys
import os
import csv
import sqlite3
import numpy as np
from datetime import datetime
import pandas as pd
from scipy.interpolate import interp1d

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl

from Cheetah.Template import Template

import scripts2d.utils as u
from wde.estimator import WaveletDensityEstimator
from wde.simple_estimator import SimpleWaveletDensityEstimator
from wde.thresholding import soft_threshold, soft_threshold_2, ultrasoft_threshold, block_threshold

from statsmodels.nonparametric.kernel_density import KDEMultivariate

from steps.common import connect, sample_fname, calc_true_pdf, calc_spwe_ise, grid_points

#
# -- data source
#

def dbname(dist_name, wave_name):
    if wave_name[0:4] == 'sim-':
        wave_name = wave_name[4:]
    return 'data/STEPS/%s/db/%s/data-ise.db' % (dist_name, wave_name)

def sample_name(dist_name, wave_name, fname):
    return 'data/RESP/%s/%s/%s' % (dist_name, wave_name, fname)

def exec_gen(conn, sql, args=()):
    cur = conn.execute(sql, args)
    row = cur.fetchone()
    while row is not None:
        yield row
        row = cur.fetchone()

class TemplateFile(object):
    def __init__(self, fname):
        self._data = {}
        with open(fname, 'r') as fh:
            self._template = Template(fh.read(), searchList=self._data)
    def render(self, data):
        self._data.clear()
        self._data.update(data)
        return str(self._template)
#
# -- generate true plots
#
def generate_plot(dist_code):
    fname = 'data/plots-tex/true-%s.eps' % dist_code
    dist = u.dist_from_code(dist_code)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    X = np.linspace(0.0,1.0, num=75)
    Y = np.linspace(0.0,1.0, num=75)
    XX, YY = np.meshgrid(X, Y)
    Z = dist.pdf((XX, YY))
    # see http://mpastell.com/2013/05/02/matplotlib_colormaps/
    surf = ax.plot_surface(XX, YY, Z, edgecolors='k', linewidth=0.5, cmap=cm.get_cmap('BuGn'))
    #ax.set_zlim(0, 5)
    plt.savefig(fname, pad_inches=0.0, orientation='portrait', frameon=False)
    plt.close(fig)
    return 'true-%s.pdf' % dist_code

def tex_figure(data):
    template = TemplateFile('steps/templates/figures.tex.qtl')
    data = dict(
        label='true_figures',
        figures=data,
        width=0.95/len(data)
    )
    with open('data/plots-tex/figures.tex', 'w') as fh:
        fh.write(template.render(data))

def generate_true_plots():
    data = []
    for code, title in [('mult', 'Gaussian mix 1'), ('mix2','Gaussian mix 2'), ('mix8', '2D Comb')]:
        fname = generate_plot(code)
        data.append(dict(label=code, fname=fname, caption=title))
    tex_figure(data)

def generate_comparison_plot(label, dist):
    fname = 'data/plots-tex/comp-%s.eps' % label
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.view_init(elev=0, azim=-45)
    ax.set_zlim(-2,16)
    X = np.linspace(0.0,1.0, num=256)
    Y = np.linspace(0.0,1.0, num=256)
    XX, YY = np.meshgrid(X, Y)
    Z = dist.pdf((XX, YY))
    # see http://mpastell.com/2013/05/02/matplotlib_colormaps/
    surf = ax.plot_surface(XX, YY, Z, edgecolors='k', linewidth=0.5, cmap=cm.get_cmap('BuGn'))
    #ax.set_zlim(0, 5)
    plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    plt.savefig(fname, pad_inches=0.0, orientation='portrait', frameon=False, bbox_inches='tight')
    #plt.show()
    plt.close(fig)
    return 'comp-%s.pdf' % label

def generate_comparison_plots():
    sql = """select fname
        from results
        where n = ? and k = 1 and j0 = ?
        order by ise asc
        limit 1"""
    dist_code = 'mult'
    n = 4096
    waves = (
        ('db6', 0.314, 'new_wde', 'New estimator'),
        ('sim-db6', 0.425, 'classic_wde', 'Classic estimator')
        )
    dists = []
    dists.append(dict(
        label='comp_true',
        dist=u.dist_from_code(dist_code),
        title='True density'
        ))
    with connect(dist_code) as conn:
        for row in exec_gen(conn, sql, (n, 3)):
            fname = sample_name(dist_code, 'db6', row[0])
            data = u.read_sample(fname)
    print(fname, len(data))
    wde = WaveletDensityEstimator('db6', k = 1, j0 = 3, j1 = 0)
    wde.fit(data)
    dists.append(dict(label='new_wde', dist=wde, title='New estimator'))
    wde = SimpleWaveletDensityEstimator('db6', j0 = 3, j1 = 0)
    wde.fit(data)
    dists.append(dict(label='classic_wde', dist=wde, title='Classic'))
    data = []
    for dist_desc in dists:
        fname = generate_comparison_plot(dist_desc['label'], dist_desc['dist'])
        data.append(dict(label=dist_desc['label'], fname=fname, caption=dist_desc['title']))
    template = TemplateFile('scripts2d/templates/figures-compare.tex.qtl')
    data = dict(
        label='comparison_figures',
        figures=data,
        width=0.95/len(data)
    )
    #return
    with open('data/plots-tex/comparison-figures.tex', 'w') as fh:
        fh.write(template.render(data))

#
# -- generate tables
#
def calc_n_data_wave(dist_code, algorithm, n, wave_code, what):
    max_j = 4 if n <= 1024 or dist_code == 'mul3' else 5
    k = 1 if algorithm == 'SPWE' else 0
    sql = """select case when j0 <= j1 then j1 + 1 else j0 end as j_level, k, count(*), sum(%(what)s), sum(%(what)s * %(what)s)
        from results
        where algorithm = '%(algorithm)s' and wave_code = '%(wave_code)s'
          and n = ? and k = %(k)d and (j1 < j0) and j0 <= %(max_j)d
        group by j_level, k
        order by j_level, k""" % dict(what=what, max_j=max_j, algorithm=algorithm, wave_code=wave_code, k=k)
    data = {}
    with connect(dist_code) as conn:
        for row in exec_gen(conn, sql, (n,)):
            j_level, k, num_p, sum_p, sum2_p = row
            mean_p = sum_p / num_p
            std_p = math.sqrt(sum2_p / num_p - mean_p ** 2)
            data[j_level] = dict(
                j=j_level,
                mise='%.3f' % mean_p,
                first=False,
                best=False,
                _mise=mean_p,
                tick=False
            )
    return data

def get_kde(dist_code, n, what):
    with connect(dist_code) as conn:
        sql = """select count(1), sum(%(what)s), sum(%(what)s * %(what)s)
            from results
            where algorithm = 'KDE' and n = ?""" % dict(what=what)
        rows = list(exec_gen(conn, sql, (n,)))
        num_p, sum_p, sum_p2 = rows[0]
        mean = sum_p / num_p
        std = math.sqrt(sum_p2 / num_p - mean * mean)
        return mean, 0 ##0.6745 * std

def calc_n_data(dist_code, n, what):
    kde, kde_confidence = get_kde(dist_code, n, what)
    if what == 'ise':
        print(dist_code, n, kde, kde_confidence)
    data = calc_n_data_wave(dist_code, 'SPWE', n, 'db10', what)
    data_simple = calc_n_data_wave(dist_code, 'CLWE', n, 'db10', what)
    js = [data[j] for j in sorted(data.keys())]
    js[0]['first'] = True
    best_j = sorted(js, key=lambda v: v['mise'])[0]
    best_j['best'] = True
    if kde - kde_confidence <= best_j['_mise'] <= kde + kde_confidence:
        best_j['tick'] = True
    sorted(data_simple.values(), key=lambda v: v['mise'])[0]['best'] = True
    for data_j in js:
        data_j['simple_mise'] = data_simple[data_j['j']]['mise']
        data_j['simple_best'] = data_simple[data_j['j']]['best']
    return dict(
        js=js,
        len_js=len(js),
        n=n,
        kde='%.3f' % kde,
        kde_confidence='%.3f' % kde_confidence
    )

def calc_table_data(dist_code, title, what):
    ns=[]
    nn = [128,256,512,1024,2048,4096]
    if dist_code != 'mul3':
        nn.append(8192)
    for n in nn:
        ns.append(calc_n_data(dist_code, n, what))
    return dict(
        title=title,
        ns=ns
    )
    
def generate_tables_mise():
    master = TemplateFile('steps/templates/tables-master.tex.qtl')
    master_table = {}
    for what in ['ise', 'hd']:
        template = TemplateFile('steps/templates/tables.tex.qtl')
        data = []
        for code, title in [('mult', 'Gaussian mix (a)'), ('mix2','Gaussian mix (b)'), ('mix8', 'Comb (c)'), ('mul3', 'Gaussian 3D mix') ]: #, ('mul3', 'Gaussian 3D mix (d)')
            data.append(calc_table_data(code, title, what))
        data = dict(
            data=data,
            width=1.90/len(data),
            caption='Comparison for %s' % what.upper(),
            label='mise_%s' % what
            )
        master_table[what] = template.render(data)
    with open('data/plots-tex/tables-master.tex', 'w') as fh:
        data = dict(
            table_ise=master_table['ise'],
            table_hd=master_table['hd'],
            datetime=str(datetime.now())
            )
        temp_str = master.render(data)
        fh.write(temp_str)
    #with open('data/plots-tex/table')

#
#.-- generate contours
#
def generate_threshold_contours(dist_code):
    np.random.seed(1)
    dist = u.dist_from_code('mix2')
    data = dist.rvs(1024)
    dists = []
    dists.append(dict(
        label='comp_true',
        dist=u.dist_from_code(dist_code),
        title='True density'
        ))
    wde = WaveletDensityEstimator('db10', k = 1, j0 = 3, j1 = 4)
    wde.fit(data)
    dists.append(dict(label='new_wde', dist=wde, title='New estimator'))
    wde.fit(data)
    dists.append(dict(label='classic_wde', dist=wde, title='Classic'))
    data = []
    for dist_desc in dists:
        fname = generate_comparison_plot(dist_desc['label'], dist_desc['dist'])
    

def generate_other_statistics():
    pass

def contour_plot_it(dist, data, fname=None):
    fig = plt.figure()
    X = np.linspace(0.0,1.0, num=75)
    Y = np.linspace(0.0,1.0, num=75)
    XX, YY = np.meshgrid(X, Y)
    Z = dist.pdf((XX, YY))
    cs = plt.contour(XX, YY, Z, alpha=0.7, levels=np.linspace(0,16,10))
    plt.clabel(cs, inline=1, fontsize=10)
    #plt.scatter(data[:,0], data[:,1], s=2, alpha=0.0001)
    if fname is not None:
        plt.savefig('data/plots-tex/%s' % fname, pad_inches=0.0, orientation='portrait', frameon=False, bbox_inches='tight')
    plt.clf()
    #plt.show()

def soft_threshold_initial(wde):
    tt1 = max([beta for j in range(wde.j0, wde.j1+1) for beta in wde.get_betas(j)])
    print('Max \\beta', tt1)
    tt1 = 2 * tt1 / (math.sqrt(wde.j1 - wde.j0 + 1) / math.sqrt(wde.n))
    return 0, tt1


def soft_grid(grid_intval):
    tt0, tt1 = grid_intval
    grid_n = 10
    grid_dn = 3
    for i in range(grid_n):
        tt_i = tt0 + (tt1 - tt0) * i / (grid_n - 1)
        yield i, tt_i

def num_threshold_initial(wde):
    tt1 = max(wde.get_nums())
    print('Max nums', tt1)
    return 0, tt1

def num_grid(grid_intval):
    tt0, tt1 = grid_intval
    if tt0 == tt1:
        yield 0, tt0
        raise StopIteration
    grid_n = 16
    for i in range(grid_n):
        tt_i = int(tt0  + (tt1 - tt0) * i / grid_n)
        yield i, tt_i

def threshold_calc(wde, true_pdf, th_fun, th_str, initial_fun, grid_fun):
    # soft thresholding calculation
    grid_intval = initial_fun(wde)
    best_diff = 1
    best_ise = calc_spwe_ise(wde.dim, wde, true_pdf, 64)
    print('Initial ISE at j0=%d, j1=%d' % (wde.j0, wde.j1), best_ise)
    iterations = 0
    while iterations < 3 or best_diff > 0.00001:
        print(grid_intval, '->', best_ise, ',', best_diff)
        new_best = best_ise
        best_i = 0
        grid_vals = list(grid_fun(grid_intval))
        for i, tt_i in grid_vals:
            wde.thresholding = th_fun(tt_i)
            wde.pdf = wde.calc_pdf()
            new_ise = calc_spwe_ise(wde.dim, wde, true_pdf, 64)
            if new_ise < new_best:
                print(i, tt_i, '>', new_ise)
                new_best = new_ise
                best_tt = tt_i
                best_i = i
        best_diff = best_ise - new_best
        best_ise = new_best
        grid_intval = (grid_vals[max(best_i - 2,0)][1], grid_vals[min(best_i + 2, len(grid_vals))][1])
        iterations += 1
    print('Best ISE at j0=%d, j1=%d for %s :' % (wde.j0, wde.j1, th_str), best_ise)
    return best_tt

def generate_plots_threshold(dist_code, wave):
    np.random.seed(1)
    nn = 4096
    print('>>>>', wave, nn)
    j0_ini = 4
    j0 = j0_ini - 3 # 3, 1
    j1 = j0 + 2 #4
    th_fun, th_str = block_threshold, 'block_threshold'

    # ttt = 0.7
    # fig = plt.figure()
    # th = ultrasoft_threshold(ttt)
    # fn = lambda coeff: th(nn, 3, 3, 2, coeff)
    # rr = np.linspace(-int(5 * ttt), int(5 * ttt), 40)
    # yy = [fn(y) for y in rr]
    # plt.plot(rr, yy)
    # plt.plot([min(rr), max(rr)],[min(rr), max(rr)], 'r--', alpha=0.3)
    # plt.plot([min(rr), max(rr)],[min(rr) - ttt, max(rr) - ttt], 'k--', alpha=0.3)
    # plt.plot([min(rr), max(rr)],[min(rr) + ttt, max(rr) + ttt], 'k--', alpha=0.3)
    # plt.show()
    #
    # sys.stdin.read()
    # return

    dist = u.dist_from_code(dist_code)
    t0 = datetime.now()
    fname = sample_fname(dist_code, nn, 1)
    sample = np.genfromtxt(fname, delimiter=',')
    points = grid_points(dist.dim, 64)
    true_pdf_64 = dist.pdf(points)

    # kde = KDEMultivariate(sample, 'c' * sample.shape[1], bw='cv_ml')
    # ise = calc_spwe_ise(sample.shape[1], kde, true_pdf_64, 64)
    # print('Bandwidth', kde.bw)
    # print('KDE ISE', ise)


    wde = WaveletDensityEstimator(wave, k=1, j0=j0_ini, j1=None)
    wde.fit(sample)
    ise = calc_spwe_ise(wde.dim, wde, true_pdf_64, 64)
    elapsed_time = (datetime.now() - t0).total_seconds()
    print('Running', th_str, 'algorithm')
    print('ISE at j0=4', ise, ' (', elapsed_time, ' secs)')

    contour_plot_it(dist, sample, 'threshold-%s-%d-true.eps' % (dist_code, nn))
    contour_plot_it(wde, sample, 'th_%s-%s-%d-wde-raw-j0=%d.eps' % (dist_code, th_str, nn, wde.j0))
    points = grid_points(dist.dim, 64)
    true_pdf_64 = dist.pdf(points)
    wde = WaveletDensityEstimator(wave, k=1, j0=j0, j1=j1)
    wde.fit(sample)
    contour_plot_it(wde, sample, 'th_%s-%s-%d-wde-raw-j=%d,%d.eps' % (dist_code, th_str, nn, wde.j0, wde.j1))
    threshold = threshold_calc(wde, true_pdf_64, th_fun, th_str, num_threshold_initial, num_grid)
    wde.thresholding = th_fun(threshold)
    wde.calc_pdf()
    contour_plot_it(wde, sample, 'th_%s-%s-%d-wde-threshold-j=%d,%d,th=%f.eps' % (dist_code, th_str, nn, wde.j0, wde.j1, threshold))
    elapsed_time = (datetime.now() - t0).total_seconds()
    print('Threshold', threshold, ', secs', elapsed_time)

    #betas = np.array(flatten([wde.get_betas(j) for j in wde.coeffs.keys()]))
    # standard soft_threshold
    #print 'all coeffs #',betas.shape
    return

def generate_all():
    #plt.ioff()
    #generate_true_plots()
    #generate_comparison_plots()
    #generate_tables_mise()
    #generate_plots_threshold('mix8', 'db10')
    #generate_other_statistics()
    print('Done')

if __name__ == "__main__":
    generate_all()