#!/bin/bash

for fname in data/plots-tex/{th_mix8-soft-8192-wde-threshold-j=4\,4\,th=2.687919,th_mix8-soft-8192-wde-raw-j=4\,4}.eps; do
  echo "$fname -> ${fname%.*}.pdf"
  gs -dNOPAUSE -dBATCH -q -sDEVICE=pdfwrite -dEPSCrop -sOutputFile="${fname%.*}.pdf" -f "$fname"
done