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

for dist in mix5 mix8; do
  echo '------'
  echo " $dist "
  echo '------'
  python steps/data-ise-load.py SPWE $dist data/STEPS/$dist/spwe/db10/*/{00128,00256,00512}/ise-512-000.csv
  python steps/data-ise-load.py SPWE $dist data/STEPS/$dist/spwe/db10/*/{01024,02048,04096,08192}/ise-128-*.csv
  if [ $dist == mix8 ]; then
    python steps/data-ise-load.py CLWE $dist data/STEPS/$dist/clwe/db10/*/*/ise-*.csv
  fi
  python steps/data-ise-load.py KDE $dist data/STEPS/$dist/kde_ise_hd/{00128,00256,00512}/kde-512-000.csv
  python steps/data-ise-load.py KDE $dist data/STEPS/$dist/kde_ise_hd/{01024,02048,04096,08192}/kde-128-*.csv
done

dist='mul3'
echo '------'
echo " $dist "
echo '------'
python steps/data-ise-load.py SPWE $dist data/STEPS/$dist/spwe/db10/*/{00128,00256,00512}/ise-128-*.csv
python steps/data-ise-load.py SPWE $dist data/STEPS/$dist/spwe/db10/*/{01024,02048,04096}/ise-032-*.csv
python steps/data-ise-load.py CLWE $dist data/STEPS/$dist/clwe/db10/*/*/ise-*.csv
python steps/data-ise-load.py KDE $dist data/STEPS/$dist/kde_ise_hd/{00128,00256,00512}/kde-*-000.csv
python steps/data-ise-load.py KDE $dist data/STEPS/$dist/kde_ise_hd/{01024,02048,04096}/kde-*-*.csv
