##ReMatCh - Read Mapping and Check-Coverage

```
rematch.py [-h] [-r [/path/reference.fasta]] [-d [/path/to/workdir]]
                  [-gatk [/path/to/gatk.jar]] [-picard [/path/to/picard]]
                  [-l [/path/to/idenfifiersList.txt]] [-cov [N]] [-qual [N]]
                  [-mul [0.0 - 1.0]] [-threads [N]]
                  [-tax [Streptococcus pneumoniae]] [-xtraSeq [XTRASEQ]]
                  [-bowtieBuild] [-clean] [-rmFastq] [-allplat]
                  [--mergeResults /path/to/workdir]
                  [--sequenceCoverage 0.0 - 1.0]
```
**ReMatCh** is an application which combines a set of bioinformatic tools for
reads mapping against a reference, finds the allelic variants and produces a
consensus sequence. It also allows direct sample download from ENA database to
be used in the analysis.

    optional arguments:
      -h, --help            show this help message and exit
      -cov [N], --minCoverage [N]
                            Minimum coverage depth required for base calling and
                            SNP calling. (default: 10)
      -qual [N], --minQuality [N]
                            Minimum mapping quality for SNP calling (default: 10)
      -mul [0.0 - 1.0], --multipleAlleles [0.0 - 1.0]
                            Minimum reads frequency (confidence) of dominant
                            nucleotide to consider absence of multiple alleles at
                            a given SNP position. (default: 0.75)
      -threads [N]          Number of threads used to run bowtie2 (default: 1)
      -tax [Streptococcus pneumoniae]
                            Name taxon to download sequences. Results will be
                            stored in /path/to/idenfifiersList.txt (default: None)
      -xtraSeq [XTRASEQ]    For trimming extra sequence lenght 5' and 3' (default:
                            0)
      -bowtieBuild          Run build bowtie (default: False)
      -clean                Clean intermediate files produced by the application
                            (.bam, .vcf, index files, coverage file) (default:
                            False)
      -rmFastq              Remove fastq files after the analysis (default: False)
      -allplat              Use all platforms. By default, only Illumina runs are
                            used (default: False)

    required arguments:
      -r [/path/reference.fasta]
                            Path for the reference sequence (default: None)
      -d [/path/to/workdir], --workdir [/path/to/workdir]
                            Working directory. Downloaded files will be stored
                            here, but it can also already contain folders with
                            fastq files. Results will be stored here. (default:
                            None)
      -gatk [/path/to/gatk.jar]
                            Path for the Genome Analysis Toolkit jar file
                            (default: None)
      -picard [/path/to/picard]
                            Path for Picard (default: None)
      -l [/path/to/idenfifiersList.txt]
                            Path to a list with ids to run. IDs can be ENA run
                            accession numbers for download or direcory names where
                            fastqs are stored in --workdir. Run accession numbers
                            retrieved from ENA using -tax will be stored here.
                            (default: None)
    
    merge results arguments. To be used after ReMatCh run:
      --mergeResults /path/to/workdir
                            Merge all ReMatCh results available at --workdir.
                            Option to be used alone or with --sequenceCoverage.
                            (default: None)
      --sequenceCoverage 0.0 - 1.0
                            Minimum sequence length to consider the gene to be
                            present. This is a relative measure. To be used with
                            --mergeResults (default: 0.8)

 

 

##Output

**ReMatch** outputs two distinct files for each sample: *sample_mappingCheck.tab* and *sample_sequences.fasta*

*sample_sequences.fasta* - Multi-fasta file storing sample consensus sequences.
 *sample_mappingCheck.tab* - Tabular file with sequences statistics.

***sample_mappingCheck.tab*** information:

**Sequence** - Target sequence name.
**Duplication** - Sequence nucleotide frequency with more than 1.25x mean coverage depth.
**Indel** - Sequence nucleotide frequency with less than 0.1x mean coverage depth.
**Coverage** - Sequence nucleotide frequency with less than minimum coverage required for base calling (0 - all nucleotide present; 1 - sequence absent).
**LowSNPsQualityScore** - Presence of SNPs with low mapping quality score (< -qual).
**SNPCoverage** - Presence of SNPs with low coverage depth (< -cov).
**SNPMultipleAlleles** - Presence of SNPs with multiple alleles.
**meanSequenceCoverage** - Mean coverage depth of the entire sequence.
**meanTrimmedSequenceCoverage** - Mean coverage depth of the sequence without the trimmed ends.

**ReMatCh** *--mergeResults* also outputs a tabular file with gene presence/absence information for samples existing in the *--workdir*, with the name *mergedResults.tab*.

**meanTrimmedSequenceCoverage** are reported for each gene (in columns) for each sample (in rows) unless multiple alleles are found ('**Mul.Alleles**' reported). If the sequence length is lower than *--sequenceCoverage*, gene is consider to be absent ('**Absent**' reported).  




##Dependencies
Bowtie2 http://bowtie-bio.sourceforge.net/bowtie2/index.shtml - Reads mapping
SAMtools http://www.htslib.org/ - .sam and .bam manipulation. SNP calling.
BEDtools http://bedtools.readthedocs.org/en/latest/ - Coverage retrieval.
Genome Analysis Toolkit (GATK) https://www.broadinstitute.org/gatk/ - Consensus sequence creation.
Picard http://broadinstitute.github.io/picard/ - Reference sequence indexing. 
Python Module: NumPy http://www.numpy.org/ - Matrix manipulation.

##Authors

Miguel Machado (https://github.com/miguelpmachado)
Bruno Gonçalves (https://github.com/bfrgoncalves)
Mickael Silva (https://github.com/mickaelsilva)

> Written with [StackEdit](https://stackedit.io/).