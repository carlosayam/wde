import csv
import math
import os
import sqlite3
import sys

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

def get_ns(dist_code, wave_name):
    resp = []
    sql = """select distinct n from results order by n"""
    with connect(dist_code, wave_name) as conn:
        for row in exec_gen(conn, sql):
            resp.append(int(row[0]))
    return resp

def list_all(dist_code, wave_name):
    ns = get_ns(dist_code, wave_name)
    sql = """select case when j0 <= j1 then j1 + 1 else j0 end as j_level, k, count(*), sum(ise), sum(ise * ise)
        from results
        where n = ? and k = 1 and ((j1 >= j0 and j1 >= 2) or (j1 < j0))
        group by j_level, k
        order by j_level, k"""
    data = {}
    data_1 = {}
    with connect(dist_code, wave_name) as conn:
        for n in ns:
            for row in exec_gen(conn, sql, (n,)):
                j_level, k, num_p, sum_p, sum2_p = row
                mean_p = sum_p / num_p
                std_p = math.sqrt(sum2_p / num_p - mean_p ** 2)
                if n in data:
                    if data[n]['mise'] > mean_p:
                        data[n] = dict(mise=mean_p, std_ise=std_p, best=(j_level, k))
                else:
                    data[n] = dict(mise=mean_p, std_ise=std_p, best=(j_level, k))
                if k == 1:
                    data_1[(n, j_level)] = mean_p
    for n in ns:
        the_best = data[n]
        print n, the_best['best'], the_best['mise'], the_best['std_ise']

if __name__ == "__main__":
    list_all(sys.argv[1], sys.argv[2])