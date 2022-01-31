#!/usr/bin/env bash

set -e
#set -o xtrace

popvcf=$1
vcf=$2

if [[ -z $vcf ]]; then
  echo "Usage: $0 <popvcf_binary> <vcf>" >&2
  exit 1
fi

if [[ ${vcf} != "test.vcf" ]]; then
  cp $vcf test.vcf
fi

#echo "cat test.vcf | bgzip -c > test.pipe.vcf.gz"
#time cat test.vcf | bgzip -c > test.pipe.vcf.gz

run ()
(
  echo "$1"
  eval time "$1"
)

#echo "bgzip -k -f test.vcf"
#time bgzip -k -f test.vcf
echo "== Compression times =="
run "bgzip -c -k -f test.vcf > test.vcf.gz"
#run "bgzip -c -k -f -l9 test.vcf > test.vcf.l9.gz"
run "${popvcf} encode test.vcf > test.vcf.popvcf"
#run "${popvcf} encode test.vcf.gz > test.vcf.gz.popvcf"
run "${popvcf} encode test.vcf -Oz > test.vcf.popvcf.gz"
run "${popvcf} encode test.vcf -Oz -l9 > test.vcf.popvcf.l9.gz"
run "spvcf encode --quiet --no-squeeze test.vcf > test.vcf.spvcf"
run "bgzip -c -k -f test.vcf.spvcf > test.vcf.spvcf.gz"
run "bgzip -c -k -f -l9 test.vcf.spvcf > test.vcf.spvcf.l9.gz"
run "${popvcf} encode test.vcf.spvcf -Oz -o test.vcf.popspvcf.gz"
run "${popvcf} encode test.vcf.spvcf -l9 -Oz -o test.vcf.popspvcf.l9.gz"

echo "== Decompression times =="
md5sum test.vcf
run "${popvcf} decode test.vcf.popvcf > test.vcf.popvcf.vcf"
md5sum test.vcf.popvcf.vcf ; rm -f test.vcf.popvcf.vcf
run "${popvcf} decode test.vcf.popvcf.gz > test.vcf.popvcf.gz.vcf"
md5sum test.vcf.popvcf.gz.vcf ; rm -f test.vcf.popvcf.gz.vcf
run "${popvcf} decode test.vcf.popvcf.l9.gz > test.vcf.popvcf.l9.gz.vcf"
md5sum test.vcf.popvcf.l9.gz.vcf ; rm -f test.vcf.popvcf.l9.gz.vcf
#run "bgzip -c -d -k -f test.vcf.gz > test.vcf.gz.vcf"
#md5sum test.vcf.gz.vcf ; rm -f test.vcf.gz.vcf


echo "== Index construction times =="
run "tabix -f test.vcf.gz"
run "tabix -f test.vcf.popvcf.gz"
run "tabix -f test.vcf.popvcf.l9.gz"

echo "== Query times =="
run "tabix test.vcf.gz chr20:17625500-17625600 > /dev/null"
run "tabix test.vcf.popvcf.gz chr20:17625500-17625600 | ${popvcf} decode - > /dev/null"
#run "tabix test.vcf.popvcf.l9.gz chr20:17625500-17625600 | ${popvcf} decode - > /dev/null"

ls -lh test.*gz
ls -l test.*gz

original_size=$(find -L . -name "test.vcf" -printf "%s\n")
find . -name "test.*gz" -printf "%f\t%s\n" | awk -v os="${original_size}" '{print $1"\t"$2"\t"os/$2}'

# cleanup
#echo test.* | tr ' ' '\n' | grep -vP "^test.vcf$" | grep -vP "^test.vcf.gz$" | xargs rm
