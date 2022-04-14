#!/bin/bash

if [[ $# -ne 7 ]]; then
  echo "Usage: $(basename $0) <sequence fasta prefix> <ntHits BF filenames> <k values> <X parameter> <Y parameter> <num threads> <output file>"
  exit 1
fi

set -e -o pipefail

seqs_basename=$1; shift
bfs=($1); shift
k=($1); shift
X=$1; shift
Y=$1; shift
t=$1; shift
outfile=$1; shift

input_size=$(wc -c ${seqs_basename}.fa | awk '{ print $1 }')

for i in "${!k[@]}"; do
  # Set input sequence for ntEdit
  if test -z "$prev"; then
    input=${seqs_basename}
  else
    input=${prev}
  fi
  ntedit -f ${input}.fa -r ${bfs[$i]} -d5 -i5 -m1 -X$X -Y$Y -b ${input}.k${k[$i]}.X$X.Y$Y -t$t -a1
  prev=${input}.k${k[$i]}.X${X}.Y${Y}_edited
done

# If ntEdit made a mistake and trimmed the hell outta the input sequence, just ignore the output
output_size=$(wc -c ${prev}.fa | awk '{ print $1 }')
size_diff=$((input_size-output_size))
if [[ ${size_diff} -gt 5000 ]]; then
  ln -sf ${seqs_basename}.fa ${outfile}
  echo "run_ntedit.sh: skipped ${seqs_basename}.fa"
else
  ln -sf ${prev}.fa ${outfile}
fi
