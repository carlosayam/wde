import csv
import sqlite3
import argparse
from steps.common import connect


def read_rows(fcsv):
    for row in fcsv:
        if len(row) == 0 or len(row[0]) == 0 or row[0] == 'fname':
            continue
        try:
            # fname, dist_code, wave_code, n, j0, j1, k, ise, hs, elapsed_time
            yield row
        except:
            print('E:',fcsv.line_num, len(row))

def insert_row(conn, algorithm, row):
    fname, wave_code, n, j0, j1, k, ise, hd, etime = \
        row[0], row[2], int(row[3]), int(row[4]), int(row[5]), int(row[6]), float(row[7]), float(row[8]), float(row[9])
    try:
        conn.execute(
            """insert into results 
                (fname, algorithm, wave_code, n, j0,
                 j1, k, ise, hd, etime)
            values (?,?,?,?,?,?,?,?,?,?)""",
                (fname, algorithm, wave_code, n, j0,
                 j1, k, ise, hd, etime)
        )
    except sqlite3.Error as e:
        print(e)
        print(fname, n)
        raise


def load_data(algorithm, dist_name, fnames, delete=False):
    with connect(dist_name, True) as conn:
        if delete:
            print('Deleting data')
            conn.execute('delete from results')
        for fname in fnames:
            print(fname, end='')
            with open(fname, 'r') as f:
                fcsv = csv.reader(f, delimiter=',', skipinitialspace=True)
                for row in read_rows(fcsv):
                    insert_row(conn, algorithm, row)
            print('.')
    print('Done')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="database loading")
    parser.add_argument('-d', '--delete', action='store_true', default=False)
    parser.add_argument('algorithm', help='Algorithm', choices=['SPWE', 'CLWE', 'KDE'])
    parser.add_argument('dist_code', help='code for distribution')
    parser.add_argument('fnames', help='fname to add', nargs='+')
    args = parser.parse_args()

    load_data(args.algorithm, args.dist_code, args.fnames, args.delete)
