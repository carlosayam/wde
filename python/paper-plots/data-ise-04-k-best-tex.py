import csv
import math
import os
import sqlite3
import sys
import scipy.stats
import pandas

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

def list_all(dist_code, wave_name):
    sql = """select n, case when j0 <= j1 then j1 + 1 else j0 end as j_level, k, count(*) as nump, sum(ise) as ise, sum(ise * ise) as ise2
        from results
        group by n, j_level, k
        order by n, j_level, k"""
    with connect(dist_code, wave_name) as conn:
        df = pandas.read_sql(sql, conn)
    df['mise'] = df.ise / df.nump
    df['std_ise'] = df.ise2 / df.nump - (df.mise ** 2)
    idx = df.groupby(['n'])['mise'].transform(min) == df['mise']
    df2 = df[idx]
    for index, row in df.iterrows():
        n, best_j, best_k, mise, std_ise = row['n'], row['j_level'], row['k'], row['mise'], row['std_ise']

if __name__ == "__main__":
    list_all(sys.argv[1], sys.argv[2])