ReMatCh
=======
**ReMatCh** - *Reads Mapping and Check-Coverage*

ReMatCh is an application which combines a set of bioinformatics tools for reads mapping against a reference, finds the allelic variants and produces a consensus sequence. It also allows direct sample download from ENA database to be used in the analysis.

----------

ReMatCh **antibiotics_book_chapter** branch
-------------------------------------------

This ReMatCh branch is for use illustration in "Epidemiological Surveillance and Typing Methods to Track Antibiotic Resistant Strains Using High Throughput Sequencing" chapter from *Antibiotics - Methods and Protocols* book (Methods Molecular Biology Series, Springer, 2016).

Therefore, **It will not be subject to further updates!**

<https://github.com/miguelpmachado/ReMatCh/tree/antibiotics_book_chapter>

----------


Dependencies
------------
 - [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) >= v2.2.6
 - [Picard](http://broadinstitute.github.io/picard/) = v1.2
 - [Samtools](http://www.htslib.org/) = v1.2
 - [Bcftools](http://www.htslib.org/) = v1.2
 - [Bedtools](http://bedtools.readthedocs.io/en/latest/) >= v2.22
 - [GATK](https://software.broadinstitute.org/gatk/) = v1.2
 - *Java JDK* >= v1.8


Installation
------------
    git clone -b antibiotics_book_chapter https://github.com/miguelpmachado/ReMatCh.git


Usage
-----
    usage: rematch.py [-h] [--version] {ReMatCh,mergeResults} ...


**rematch.py ReMatCh** module

*rematch.py ReMatCh* module will execute ReMatCh main scripts: it will run and analyse each sample


    usage: rematch.py ReMatCh [-h] -r /path/reference.fasta -d /path/to/workdir -l
                              [/path/to/identifiersList.txt] [-cov [N]]
                              [-qual [N]] [-mul [0.0 - 1.0]] [-j [N]]
                              [--tax "Streptococcus pneumoniae"]
                              [--xtraSeq [XTRASEQ]]
                              [--asperaKey /path/to/asperaweb_id_dsa.openssh]
                              [-bowtieBuild] [-clean] [-rmFastq] [-allplat]
                              [--useOmicsDataType GENOMIC,TRANSCRIPTOMIC]

    optional arguments:
      -h, --help            show this help message and exit

    ReMatCh required options:
      -r /path/reference.fasta, --reference /path/reference.fasta
                            Path for the reference sequence (default: None)
      -d /path/to/workdir, --workdir /path/to/workdir
                            Working directory. Downloaded files will be stored
                            here under sampleID/fastq/, but it can also already
                            contain folders with fastq files. Results will be
                            stored here. (default: None)
      -l [/path/to/identifiersList.txt]
                            Path to a list with ids to run. IDs can be ENA run
                            accession numbers for download or directory names
                            where fastq files are stored (inside a folder named
                            fastq) within --workdir directory. Run accession
                            numbers retrieved from ENA using --tax option will
                            be stored here. (default: None)

    ReMatCh optional options:
      -cov [N], --minCoverage [N]
                            Minimum coverage depth required for base calling and
                            SNP calling. (default: 10)
      -qual [N], --minQuality [N]
                            Minimum mapping quality for SNP calling (default: 10)
      -mul [0.0 - 1.0], --multipleAlleles [0.0 - 1.0]
                            Minimum reads frequency (confidence) of dominant
                            nucleotide to consider absence of multiple alleles at
                            a given SNP position. (default: 0.75)
      -j [N], --threads [N]
                            Number of threads used to run bowtie2 (default: 1)
      --tax "Streptococcus pneumoniae"
                            Name taxon to download sequences. Results will be
                            stored in /path/to/idenfifiersList.txt (default:
                            [None])
      --xtraSeq [XTRASEQ]   For trimming extra sequence lenght 5' and 3' (default:
                            0)
      --asperaKey /path/to/asperaweb_id_dsa.openssh
                            Tells ReMatCh to download run files from ENA using
                            Aspera Connect. The path to Private-key file
                            asperaweb_id_dsa.openssh normaly found in
                            ~/.aspera/connect/etc/asperaweb_id_dsa.openssh needs
                            to be provided. (default: [None])
      -bowtieBuild          Run build bowtie (default: False)
      -clean                Clean intermediate files produced by the application
                            (.bam, .vcf, index files, coverage file) (default:
                            False)
      -rmFastq              Remove fastq files after the analysis (default: False)
      -allplat              Use all platforms. By default, only Illumina runs are
                            used (default: False)
      --useOmicsDataType GENOMIC,TRANSCRIPTOMIC
                            Tells ReMatCh to analyse these OMICS data type.
                            Possible choises are
                            GENOMIC,TRANSCRIPTOMIC,METAGENOMIC,SYNTHETIC,OTHER,ALL
                            (default: ['ALL'])


**rematch.py mergeResults** module

*rematch.py mergeResults* module will execute ReMatCh main scripts: it will run and analyse each sample


    usage: rematch.py mergeResults [-h] --mrWorkdir /path/to/workdir
                                   [--mrSequenceCoverage 0.0 - 1.0]
                                   [--mrOutdir /path/to/Results/Outdir/]

    optional arguments:
      -h, --help            show this help message and exit

    mergeResults required options:
      --mrWorkdir /path/to/workdir
                            Merge all ReMatCh results available at --mrWorkdir.
                            (default: None)

    mergeResults optional options:
      --mrSequenceCoverage 0.0 - 1.0
                            Minimum sequence length to consider the gene to be
                            present. This is a relative measure. (default: [0.8])
      --mrOutdir /path/to/Results/Outdir/
                            Specify a different directory from --mrWorkdir to
                            output the merged results (otherwise merged results
                            will be stored in --mrWorkdir). (default: [None])


Output
------


**rematch.py ReMatCh** module

ReMatch main module outputs two distinct files for each sample: *sample_mappingCheck.tab* and *sample_sequences.fasta*


*sample_sequences.fasta* - Multi-fasta file containing sample consensus sequences.

*sample_mappingCheck.tab* - Tabular file with sequences statistics.


*sample_mappingCheck.tab* information:

 - *Sequence* - Target sequence name.
 - *Duplication* - Sequence nucleotide frequency with more than 1.25x mean coverage depth.
 - *Indel* - Sequence nucleotide frequency with less than 0.1x mean coverage depth.
 - *Coverage* - Sequence nucleotide frequency with less than minimum coverage required for base calling (0 - all nucleotide present; 1 -  sequence absent).
 - *LowSNPsQualityScore* - Presence of SNPs with low mapping quality score (< -qual).
 - *SNPCoverage* - Presence of SNPs with low coverage depth (< -cov).
 - *SNPMultipleAlleles* - Presence of SNPs with multiple alleles.
 - *meanSequenceCoverage* - Mean coverage depth of the entire sequence.
 - *meanTrimmedSequenceCoverage* - Mean coverage depth of the sequence without the trimmed ends.


**rematch.py mergeResults** module

Outputs a tabular file (*mergedResults.tab*) with gene presence/absence information for samples existing in ReMatCh main module *--workdir* directory


*mergedResults.tab* information:

The *rematch.py mergeResults* module will create the mergedResults.tab file inside ​~/methods_protocols/antibiotic_resistance/rematch_run_AR/merged_results/ folder that will report which genes are present in the different samples. In the case of genes being present (genes with equal to or more than 85 % of nucleotides with 10 reads coverage minimum, *--mrSequenceCoverage* option), the script will provide the mean sequence coverage, otherwise will report “Absent” for genes not present, or "Mul_Allele" for those genes that might have multiple alleles (depending on ReMatCh *-mul* option).


Authors
-------

 - Miguel Machado (https://github.com/miguelpmachado)
 - Bruno Gonçalves (https://github.com/bfrgoncalves)
 - Mickael Silva (https://github.com/mickaelsilva)

> Written with [StackEdit](https://stackedit.io/).
