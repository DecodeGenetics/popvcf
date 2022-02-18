## popVCF

popVCF losslessly encodes a multi sample VCF to reduce disk footprint. VCF fields are encoded by pointing to other exactly identical fields in the same row or in the row above. popVCF compression performance is small on a single sample VCF, but the compression ratio can go up to 40+ on a large population VCFs or 5x more compressed than the standard bgzip compression. The compression ratio varies a lot between data sets, see below for benchmarks on several different data sets.

Files are encoded with the "popvcf encode" command, and by encoding with the "-Oz" flag you can directly write the output in bgzip format. You can then decode the file back to VCF using the "popvcf decode" command. The decode subcommand can also query a region using option "--region=chrN:A-B".

On a 64 bit linux, you can get the latest static binary from the [Release page](https://github.com/DecodeGenetics/popvcf/releases).

### Benchmarks

We have benchmarked popVCF against few other compression methods with some large population VCF data. In all experiements, we report wall clock time using /usr/bin/time and used a single CPU thread. The VCF data was read and written to a SSD disk. spVCF was run with the "--no-squeeze" option to prevent any lossy compression. The script run to benchmark is in the benchmark/ directory.  In the WGS benchmarks, we had to exclude genozip and VCFShark as they were unable to compress the data because of repeated runtime errors.

Benchmarked versions: popVCF v1.1.0, spVCF v1.2.0-0-gbecb461, htslib+bcftools v1.14 (with libdeflate), Genozip 13.0.11, VCFShark v1.1.

#### GraphTyper UK biobank WGS-487k individual data

| Method/format | Compression ratio | Compared to bgzip |
| ------------- | ----------------- | ----------------- |
| popVCF+bgzip  |             37.6x |              4.4x |
| spVCF+bgzip   |             17.2x |              2.0x |
| BCF           |             10.5x |              1.2x |
| bgzip (VCF)   |              8.6x |              1.0x |

#### Deep Variant/GLnexus WES-200k individual data

| Method/format | Compression ratio | Compared to bgzip | Compression speed (MB/s) | Decompression speed (MB/s) |
| ------------- | ----------------- | ----------------- | ------------------------ | -------------------------- |
| popVCF+bgzip  |            102.9x |              6.9x |                    194.0 |                      490.7 |
| spVCF+bgzip   |             43.8x |              2.9x |                    129.7 |                      281.5 |
| Genozip       |             35.0x |              2.3x |                     18.0 |                       17.3 |
| VCFShark      |             28.3x |              1.9x |                     22.8 |                       21.7 |
| BCF           |             14.0x |             0.94x |                     62.4 |                      175.2 |
| bgzip (VCF)   |             14.9x |              1.0x |                     91.6 |                      521.3 |

#### GATK UK biobank WGS-150k individual data

| Method/format | Compression ratio | Compared to bgzip | Compression speed (MB/s) | Decompression speed (MB/s) |
| ------------- | ----------------- | ----------------- | ------------------------ | -------------------------- |
| popVCF+bgzip  |             20.1x |              2.8x |                    102.2 |                      295.0 |
| spVCF+bgzip   |             10.0x |              1.4x |                     58.8 |                      165.7 |
| BCF           |              6.7x |             0.94x |                     55.6 |                      174.2 |
| bgzip (VCF)   |              7.1x |              1.0x |                     58.5 |                      474.7 |

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