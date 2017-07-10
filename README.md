xHLA: Fast and accurate HLA typing from short read sequence data
================================================================
The Human Leukocyte Antigen (HLA) gene complex on human chromosome
6 is one of the most polymorphic regions in the human
genome, and contributes in large part to the diversity of the immune
system. Accurate typing of HLA genes with short read sequencing
data has historically been difficult due to the sequence similarity between
the polymorphic alleles.  xHLA iteratively refines the mapping results at
the amino acid level to achieve 99 to 100% 4-digit typing accuracy for both
class I and II HLA genes, taking only about 3 minutes to process a 30X
whole genome BAM file on a desktop computer.

Installation
------------
Compile docker image:

```bash
cd docker
make build
make deploy
```

Usage
-----
Run xHLA caller directly on an indexed BAM file generated using
BWA-mem against hg38 reference without alt contigs:

```bash
docker run -v `pwd`:`pwd` -w `pwd` docker-dev.hli.io/xchao/hla \
    --sample_id test --input_bam_path test.bam \
    --output_path test
```

For other types of BAMs, pre-processing is required. **Please check details
[here](https://github.com/humanlongevity/HLA/wiki/BAMs-compatible-with-xHLA).**

Output is a JSON file that lists 12 HLA alleles, 2 for each of the HLA genes:

```bash
{
 "subject_id": "176444255",
 "creation_time": "2016-05-04T08:25:04Z",
 "report_version": "1.1",
 "report_type": "hla_typing",
 "sample_id": "176444255",
 "hla": {
  "alleles": [
   "A*01:01",
   "A*02:01",
   "B*13:02",
   "B*37:01",
   "C*06:02",
   "C*06:02",
   "DPB1*04:01",
   "DPB1*04:01",
   "DQB1*02:02",
   "DQB1*05:01",
   "DRB1*07:01",
   "DRB1*10:01"
  ]
 }
}
```

Citation
--------
Xie et al. (2017) Fast and accurate HLA typing from short-read next-generation
sequence data with xHLA.
***PNAS***
[doi:10.1073/pnas.1707945114](http://www.pnas.org/content/early/2017/06/27/1707945114)

License
-------
The HLA Typing Software Code (the "Code") is made available by Human
Longevity, Inc. ("HLI") on a non-exclusive, non-sublicensable,
non-transferable basis solely for non-commercial academic research use.
Commercial use of the Code is expressly prohibited.  If you would like to obtain
a license to the Code for commercial use, please contact HLI at
bizdev@humanlongevity.com.  HLI MAKES NO REPRESENTATIONS OR WARRANTIES
WHATSOEVER, EITHER EXPRESS OR IMPLIED, WITH RESPECT TO THE CODE PROVIDED
HEREUNDER. IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR
PURPOSE WITH RESPECT TO CODE ARE EXPRESSLY DISCLAIMED. THE CODE IS FURNISHED
"AS IS" AND "WITH ALL FAULTS" AND DOWNLOADING OR USING THE CODE
IS UNDERTAKEN AT YOUR OWN RISK.  TO THE FULLEST EXTENT ALLOWED BY APPLICABLE
LAW, IN NO EVENT SHALL HLI BE LIABLE, WHETHER IN CONTRACT, TORT, WARRANTY, OR
UNDER ANY STATUTE OR ON ANY OTHER BASIS FOR SPECIAL, INCIDENTAL, INDIRECT,
PUNITIVE, MULTIPLE OR CONSEQUENTIAL DAMAGES SUSTAINED BY YOU OR ANY OTHER PERSON
OR ENTITY ON ACCOUNT OF USE OR POSSESSION OF THE CODE, WHETHER OR NOT
FORESEEABLE AND WHETHER OR NOT HLI HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
DAMAGES, INCLUDING WITHOUT LIMITATION DAMAGES ARISING FROM OR RELATED TO LOSS OF
USE, LOSS OF DATA, DOWNTIME, OR FOR LOSS OF REVENUE, PROFITS, GOODWILL, BUSINESS
OR OTHER FINANCIAL LOSS.
