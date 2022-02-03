#!/usr/bin/env bash

n=100000
echo "##fileformat=VCFv4.2"
echo "##contig=<ID=chr1>"
echo "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype.\">"
echo "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths.\">"
echo "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PHRED-scaled genotype likelihoods.\">"
echo -e -n "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"

awk -v n=${n} 'BEGIN{
  for (i = 1; i <= n; i++){
    printf "\t%08d", i
  }
}'

awk -v n=${n} -v n_alts=3 'BEGIN{
  alt="AC"
  printf "\nchr1\t"n+2
  printf "\t.\tA"
  for (a = 1; a <= n_alts; a++){
    alt=alt"C"
    printf ","alt
  }
  printf "\t0\t.\t.\t.\tGT:AD:PL"
  for (i = 1; i <= n; i++){
    printf "\t0/0:"n_alts*10
    for (a = 1; a <= n_alts; a++){
      printf ","a % 10
    }
    printf ":0"
    for (a = 2; a <= n_alts + 1; a++){
      for (b = 1; b <= a; b++){
        printf ","a+b-2
      }
    }
  }
  printf "\n"}'

awk -v n=${n} -v n_alts=3 'BEGIN{
  alt="AC"
  printf "chr1\t"n+2
  printf "\t.\tA"
  for (a = 1; a <= n_alts; a++){
    alt=alt"C"
    printf ","alt
  }
  printf "\t0\t.\t.\t.\tGT:AD:PL"
  for (i = 1; i <= n; i++){
    printf "\t0/0:"n_alts*10
    for (a = 1; a <= n_alts; a++){
      printf ","a % 10
    }
    printf ":0"
    for (a = 2; a <= n_alts + 1; a++){
      for (b = 1; b <= a; b++){
        printf ","a+b-2
      }
    }
  }
  printf "\n"}'

awk -v n=${n} -v n_alts=3 'BEGIN{
  alt="AC"
  printf "chr1\t"n+2
  printf "\t.\tA"
  for (a = 1; a <= n_alts; a++){
    alt=alt"C"
    printf ","alt
  }
  printf "\t0\t.\t.\t.\tGT:AD:PL"
  for (i = 1; i <= n; i++){
    printf "\t0/0:"n_alts*10
    for (a = 1; a <= n_alts; a++){
      printf ","a % 10
    }
    printf ":0"
    for (a = 2; a <= n_alts + 1; a++){
      for (b = 1; b <= a; b++){
        printf ","a+b-2
      }
    }
  }
  printf "\n"}'


awk -v n=${n} -v n_alts=7 'BEGIN{
  alt="GT"
  printf "chr1\t"n+3
  printf "\t.\tG"
  for (a = 1; a <= n_alts; a++){
    alt=alt"T"
    printf ","alt
  }
  printf "\t0\t.\t.\t.\tGT:AD:PL"
  for (i = 1; i <= n; i++){
    printf "\t0/0:"n_alts*10
    for (a = 1; a <= n_alts; a++){
      printf ","(a+i) % 10
    }
    printf ":0"
    for (a = 2; a <= n_alts + 1; a++){
      for (b = 1; b <= a; b++){
        printf ","a+i+b-2
      }
    }
  }
  printf "\n"}'
