import csv
import math
import os
import sqlite3
import sys
import pandas

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl

sys.path.append('.')
import scripts2d.utils as u


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

def generate_plot(dist_code):
    fname = 'data/plots-tex/%s.eps' % dist_code
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
    plt.savefig('data/paper-1/%s' % fname, pad_inches=0.0, orientation='portrait', frameon=False)
    return fname

def tex_figure(fnames):
    pass

def generate_true_plots():
    fname1 = generate_plot('mult')
    fname2 = generate_plot('mix2')
    fname3 = generate_plot('mix3')
    tex_figure([fname1, fname2, fname3])

def generate_all():
    generate_true_plots()
    generate_table_mise()
    generate_mise_contours()
    generate_other_statistics()

if __name__ == "__main__":
    generate_all()