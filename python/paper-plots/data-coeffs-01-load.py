import csv
import os
import sqlite3
import sys

def dbname(dist_name, wave_name):
    return 'data/RESP/%s/%s/data-coeffs.db' % (dist_name, wave_name)

def csvname(dist_name, wave_name):
    return 'data/RESP/%s/%s/all-coeffs.csv' % (dist_name, wave_name)

def create_table(conn):
    sql = """
    CREATE TABLE IF NOT EXISTS results (
     fname varchar(256) NOT NULL,
     n integer NOT NULL,
     j0 integer NOT NULL,
     j1 integer NOT NULL,
     k integer NOT NULL,
     j integer NOT NULL, 
     qx varchar(8) NOT NULL,
     z1 integer NOT NULL,
     z2 integer NOT NULL,
     coeff real NOT NULL
     )
    """
    conn.execute(sql)
    print('results created')

def connect(dist_name, wave_name):
    fname_db = dbname(dist_name, wave_name)
    if not os.path.isfile(fname_db):
        conn = sqlite3.connect(fname_db)
        create_table(conn)
    else:
        conn = sqlite3.connect(fname_db)
    return conn

def read_rows(fcsv):
    for row in fcsv:
        if len(row) == 0 or len(row[0]) == 0:
            continue
        try:
            # fname, dist_code, wave_code, n, j0, j1, k, j, qx, zs, coeff
            zs = row[9][1:-2].split(',')
            z1, z2 = int(zs[0]), int(zs[1])
            yield row[0], int(row[3]), int(row[4]), int(row[5]), int(row[6]), int(row[7]), row[8], z1, z2, row[10]
        except:
            print(row)
            print('E:',fcsv.line_num, len(row))
            raise

def load_data(dist_name, wave_name):
    fname_csv = csvname(dist_name, wave_name)
    with open(fname_csv, 'r') as f:
        fcsv = csv.reader(f, quotechar='"', delimiter=',')
        with connect(dist_name, wave_name) as conn:
            conn.execute('delete from results')
            headers = next(fcsv)
            print(headers)
            for fname, n, j0, j1, k, j, qx, z1, z2, coeff in read_rows(fcsv):
                try:
                    conn.execute('insert into results (fname, n, j0, j1, k, j, qx, z1, z2, coeff) values (?,?,?,?,?,?,?,?,?,?)', (fname, n, j0, j1, k, j, qx, z1, z2, coeff))
                except sqlite3.Error as e:
                    print(e)
                    print(fname, n)
                    raise
    print('Done')

if __name__ == "__main__":
    load_data(sys.argv[1], sys.argv[2])