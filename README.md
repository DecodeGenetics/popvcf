## popVCF

popVCF losslessly encodes a multi sample VCF to reduce disk footprint. VCF fields are encoded by pointing to other exactly identical fields in the same row or in the row above. popVCF compression performance is small on a single sample VCF, but the compression ratio can go up to 40+ on a large population VCF or 5x more compressed than the standard bgzip compression.

Files are encoded with the "popvcf encode" command, and by encoding with the "-Oz" flag you can directly write the output in bgzip format. You can then decode the file back to VCF using the "popvcf decode" command. The decode subcommand can also query a region using option "--region=chrN:A-B".

On a 64 bit linux, you can get the latest static binary from the [Release page](https://github.com/DecodeGenetics/popvcf/releases).

### Usage

```sh
popvcf encode my.vcf > my.popvcf
popvcf decode my.popvcf > my.new.vcf
diff my.vcf my.new.vcf # Should be the same

# It is also possible to bgzip, tabix index and query
popvcf encode my.vcf -Oz > my.popvcf.gz
tabix my.popvcf.gz
popvcf decode my.popvcf.gz > my.new2.vcf
popvcf decode my.popvcf.gz --region=chrN:A-B > my.region.vcf # Random access a region using the tabix index
```

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

### License
MIT