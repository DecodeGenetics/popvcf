## popVCF

popVCF losslessly encodes a multi sample VCF to reduce disk footprint. VCF fields are encoded by pointing to other identical fields in the same row. popVCF provives no benefits on a single sample VCF, and very little benefits until the VCF contains at least 100 samples or so.


### Building
C++17 compiler (Tested on GCC 7+, Clang 10+) is required for building popVCF.

```sh
git clone --recursive <url> popvcf # Clone the repository
cd popvcf
mkdir build-release
cd build-release
cmake ..
make -j3 popvcf
```

### Usage

```sh
cat my.vcf | popvcf encode > my.popvcf
cat my.popvcf | popvcf decode > my.1.vcf
diff my.vcf my.old.vcf # Should be the same

# It is also possible to bgzip and index the output
cat my.vcf | popvcf encode | bgzip > my.popvcf.gz
tabix my.popvcf.gz
zcat my.popvcf.gz | popvcf decode > my.2.vcf
```

### License
MIT