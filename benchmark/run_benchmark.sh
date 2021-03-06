#!/usr/bin/env bash

set -e
#set -o xtrace

popvcf=$1
vcf=$2
level="$3"
#level="-l9"

if [[ -z $vcf ]]; then
  echo "Usage: $0 <popvcf_binary> <vcf> [level]" >&2
  echo " i.e. $0 ./popvcf test.vcf \"-l9\""
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
run "bgzip --stdout --force --threads 1 test.vcf ${level} > test.vcf.gz"
run "${popvcf} encode test.vcf -o test.vcf.popvcf"
run "${popvcf} encode test.vcf ${level} -Oz -o test.vcf.popvcf.gz"
run "spvcf encode --quiet --no-squeeze test.vcf > test.vcf.spvcf"
run "spvcf encode --quiet --no-squeeze test.vcf | bgzip -c ${level} > test.vcf.spvcf.gz"

echo "== Decompression times =="
md5sum test.vcf
run "${popvcf} decode test.vcf.popvcf > test.vcf.popvcf.vcf"
md5sum test.vcf.popvcf.vcf ; rm -f test.vcf.popvcf.vcf
run "${popvcf} decode test.vcf.popvcf.gz > test.vcf.popvcf.gz.vcf"
md5sum test.vcf.popvcf.gz.vcf ; rm -f test.vcf.popvcf.gz.vcf
run "bgzip -dc test.vcf.spvcf.gz | spvcf decode --quiet > test.vcf.spvcf.gz.vcf"
md5sum test.vcf.spvcf.gz.vcf ; rm -f test.vcf.spvcf.gz.vcf
run "bgzip -dc test.vcf.gz > test.vcf.gz.vcf"
md5sum test.vcf.gz.vcf ; rm -f test.vcf.gz.vcf

echo "== Index construction times =="
run "tabix -p vcf -f test.vcf.gz"
run "tabix -p vcf -f test.vcf.popvcf.gz"
run "tabix -p vcf -f test.vcf.spvcf.gz"

echo "== Query times =="
region=$(grep -v ^# test.vcf | cut -f1,2 | head -n 20 | tail -n 1 | awk '{print $1":"$2"-"$2+100}')
run "tabix test.vcf.gz ${region} > /dev/null"
run "${popvcf} decode test.vcf.popvcf.gz --region=${region} > /dev/null"
run "spvcf tabix test.vcf.spvcf.gz ${region} | spvcf decode - > /dev/null"

ls -lh test.*gz
ls -l test.*gz

original_size=$(find -L . -name "test.vcf" -printf "%s\n")
find . -name "test.*gz" -printf "%f\t%s\n" | awk -v os="${original_size}" '{print $1"\t"$2"\t"os/$2}'

# cleanup
#echo test.* | tr ' ' '\n' | grep -vP "^test.vcf$" | grep -vP "^test.vcf.gz$" | xargs rm
