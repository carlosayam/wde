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
import pandas
from scipy.interpolate import interp1d

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl

from Cheetah.Template import Template

import scripts2d.utils as u

#
# -- data source
#

def dbname(dist_name, wave_name):
    return 'data/RESP/%s/%s/data-ise.db' % (dist_name, wave_name)

def exec_gen(conn, sql, args=()):
    cur = conn.execute(sql, args)
    row = cur.fetchone()
    while row is not None:
        yield row
        row = cur.fetchone()

def connect(dist_name, wave_name):
    fname_db = dbname(dist_name, wave_name)
    if not os.path.isfile(fname_db):
        conn = sqlite3.connect(fname_db)
        create_table(conn)
    else:
        conn = sqlite3.connect(fname_db)
    return conn

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
    template = TemplateFile('scripts2d/templates/figures.tex.qtl')
    data = dict(
        label='true_figures',
        figures=data,
        width=0.95/len(data)
    )
    with open('data/plots-tex/figures.tex', 'w') as fh:
        fh.write(template.render(data))

def generate_true_plots():
    data = []
    for code, title in [('beta', 'Beta'), ('mix2','Gaussian mixture'), ('mix3', '3D Comb')]:
        fname = generate_plot(code)
        data.append(dict(label=code, fname=fname, caption=title))
    tex_figure(data)

#
# -- generate tables
#
def calc_n_data_wave(dist_code, n, wave_name):
    sql = """select case when j0 <= j1 then j1 + 1 else j0 end as j_level, k, count(*), sum(ise), sum(ise * ise)
        from results
        where n = ? and k = 1 and ((j1 >= j0 and j1 >= 2) or (j1 < j0))
        group by j_level, k
        order by j_level, k"""
    data = {}
    with connect(dist_code, wave_name) as conn:
        for row in exec_gen(conn, sql, (n,)):
            j_level, k, num_p, sum_p, sum2_p = row
            mean_p = sum_p / num_p
            std_p = math.sqrt(sum2_p / num_p - mean_p ** 2)
            data[j_level] = dict(
                j=j_level,
                mise='%.3f' % mean_p,
                first=False,
                best=False
            )
    return data

def calc_n_data(dist_code, n):
    data = calc_n_data_wave(dist_code, n, 'db6')
    data_simple = calc_n_data_wave(dist_code, n, 'sim-db6')
    js = [data[j] for j in sorted(data.keys())[:-1]]
    js[0]['first'] = True
    sorted(js, key=lambda v: v['mise'])[0]['best'] = True
    sorted(data_simple.values(), key=lambda v: v['mise'])[0]['best'] = True
    for data_j in js:
        data_j['simple_mise'] = data_simple[data_j['j']]['mise']
        data_j['simple_best'] = data_simple[data_j['j']]['best']
    return dict(
        js=js,
        len_js=len(js),
        n=n
    )

def calc_table_data(dist_code, title):
    ns=[]
    for n in [128,256,512,1024,2048,4096]:
        ns.append(calc_n_data(dist_code, n))
    return dict(
        title=title,
        ns=ns
    )

def generate_table_mise():
    template = TemplateFile('scripts2d/templates/tables.tex.qtl')
    data = []
    for code, title in [('beta', 'Beta'), ('mix2','Gaussian mixture'), ('mix3', '2D Comb')]:
        data.append(calc_table_data(code, title))
    data = dict(
        data=data,
        width=0.95/len(data)
        )
    with open('data/plots-tex/tables.tex', 'w') as fh:
        temp_str = template.render(data)
        fh.write(temp_str)

#
#.-- generate contours
#
def generate_mise_contours():
    pass

def generate_other_statistics():
    pass

def generate_all():
    plt.ioff()
    generate_true_plots()
    generate_table_mise()
    generate_mise_contours()
    generate_other_statistics()
    print 'Done'
