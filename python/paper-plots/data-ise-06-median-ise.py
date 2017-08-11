import csv
import math
import os
import sqlite3
import sys
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
    ns = [128,256,512,1024,2048,4096]
    sql = """select n, case when j0 <= j1 then j1 + 1 else j0 end as j_level, ise
        from results
        where k = 1 and ((j1 >= j0 and j1 >= 2) or (j1 < j0))
        """
    # and ((j0 <= j1 and j1 > 2) or (j1 < j0))
    with connect(dist_code, wave_name) as conn:
        df = pandas.read_sql(sql, conn)
    df2 = df.groupby([df.n, df.j_level])[['ise']].median()
    print 'Median, k=1'
    print df2
if __name__ == "__main__":
    list_all(sys.argv[1], sys.argv[2])