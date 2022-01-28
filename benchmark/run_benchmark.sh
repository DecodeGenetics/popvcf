#!/usr/bin/env bash

set -e
#set -o xtrace

popvcf=$1
vcf=$2

if [[ -z $vcf ]]; then
  echo "Usage: $0 <popvcf_binary> <vcf>" >&2
  exit 1
fi

cp $vcf test.vcf

echo "cat test.vcf | bgzip -c > test.pipe.vcf.gz"
time cat test.vcf | bgzip -c > test.pipe.vcf.gz

echo "bgzip -k -f test.vcf"
time bgzip -k test.vcf

echo "${popvcf} encode test.vcf"
time ${popvcf} encode test.vcf > test.popvcf

echo "${popvcf} encode -Iz test.vcf.gz"
time ${popvcf} encode -Iz test.vcf.gz > test.vcfgz.popvcf

echo "${popvcf} encode test.vcf -Oz > test.vcf2popvcf.gz"
time ${popvcf} encode test.vcf -Oz > test.vcf2popvcf.gz

echo "bgzip on popVCF"
cat test.popvcf | time bgzip -c > test.popvcf.gz

#echo "zstd on popVCF"
#cat test.popvcf | time zstd > test.popvcf.zst

echo "popVCF decode"
zcat test.popvcf.gz | time ./popvcf decode | md5sum

echo "expected md5sum"
cat test.vcf | md5sum

echo "bgzip on VCF"
cat test.vcf | time bgzip -c > test.vcf.gz

#ls -lh test.*vcf*
