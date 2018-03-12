#!/bin/bash
set -ae

for dist in mult mix2 mix3; do
  echo '------'
  echo " $dist "
  echo '------'
  python steps/data-ise-load.py -d SPWE $dist data/STEPS/$dist/spwe/db10/*/*/ise-500-000.csv
  python steps/data-ise-load.py CLWE $dist data/STEPS/$dist/clwe/db10/*/*/ise-500-000.csv
  python steps/data-ise-load.py KDE $dist data/STEPS/$dist/kde_ise_hd/*/kde-512-000.csv
  python steps/data-ise-load.py KDE $dist data/STEPS/$dist/kde_ise_hd/*/kde-128-*.csv
done
