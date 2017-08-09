#!/bin/bash

for fname in data/plots-tex/*.eps; do
  echo "$fname -> ${fname%.*}.pdf"
  gs -dNOPAUSE -dBATCH -q -sDEVICE=pdfwrite -dEPSCrop -sOutputFile="${fname%.*}.pdf" -f "$fname"
done
mv data/plots-tex/*.pdf ~/dev/ME/shape-preserving/images
for fname in data/plots-tex/*.tex; do
  echo "$fname -> includes/"
  egrep -v '^[ ]*$' $fname > ~/dev/ME/shape-preserving/includes/$(basename $fname)
done
