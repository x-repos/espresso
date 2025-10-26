#!/bin/bash
mkdir -p tmp

for i in mg_relax mgh2_relax h2_scf; do
  echo ">>> Running $i ..."
  pw.x -in ${i}.in > ${i}.out
done
