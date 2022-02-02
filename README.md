## popVCF

popVCF losslessly encodes a multi sample VCF to reduce disk footprint. VCF fields are encoded by pointing to other identical fields in the same row. popVCF provives no benefits on a single sample VCF, and very little benefits until the VCF contains at least 100 samples or so.


### Building
Feature complete C++17 compiler is required for building popVCF, i.e. GCC 8/Clang 10 or newer.

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
popvcf encode my.vcf > my.popvcf
popvcf decode my.popvcf > my.new.vcf
diff my.vcf my.new.vcf # Should be the same

# It is also possible to bgzip and index the output
popvcf encode my.vcf -Oz > my.popvcf.gz
tabix my.popvcf.gz
popvcf decode my.popvcf.gz > my.new2.vcf
```

### License
MIT