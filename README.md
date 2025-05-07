# AncestrySnpsToPcaDemo

## Inroduction
Many people are interested in the breed origins of their dogs. One will often hear people say things along the lines of: "My dog is definitely part dingo", or "I'm pretty sure Rex is a wolf hybrid", or "Fritz is a Carolina dog". Given preternatural curiosity and a non-trivial chunk of money, one can now get a better handle on the genetic ancestry of your canine pet by submitting a DNA sample via [Ancestry for dogs](https://petdna.ancestry.com/). Ancestry uses a genotyping array designed for dogs, which generates genotypes at hundreds of thousands of variable sites in the dog genome. The company then performs ancestry inference, comparing your dog's genotypes against those for a large panel of bered dogs, from which the percentages of ancestry from various breeds in your dog can be estimated.

This repository demonsrates how one would analyze the data yourself, if say you and a bunch of friends got back results for your dogs and wanted to compare them. There are some peculiar features of the genotyping array data in general, and the data provided by Ancestry in particular, and I show you how to deal with these issues so that you can do two things:

* merge all of the Ancestry-generated dog genotypes from a bunch of samples
* merge these with a panel of publicly available canid samples.

## The publicly available data
In this particular demo, we will merge the dog data with genotypes obtained from whole-genome sequencing for 6 canids: dingo, Bajenji, Chinese wolf, Croatian wolf, Israeli wolf, and Golden Jackal. These data were first published alongside a study investigating the timing and geographic origins of dogs [Freedman et al. 2014, *PLoS Genetics*](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004016).
All samples were sequenced to relatively high coverage (approximately 20x each), meaning that, on average, enough sequencing reads overlap any particular site in the genome such that heterozygous sites can be genotyped with high confidence. 

## Data analysis: Ancestry samples
### PLINK files: 1st steps
PLINK is a software package for manipulating genotype array data. Unsurprisingly, PLINK also refes to a particular file format, or more specifically, a set of files that store genotypes that can be used with PLINK and other tools that can take PLINK format files as input. There are text-reader readible PLINK files (with one *ped,*fam, and *map file for each sample) and there are binary versions ("bim and *bed). As is typicaly in bioinformatics, the same prefix can be used for different software and different data analysis contexts, e.g. "bed" files also refer to flat text files that represent genomc intervals. If you decide to stay in biology and as a result, do a bit of data crunching, you will get used to this! Or, you won't and may wind up being the owner of a vegan bakery instead.

#### Counting the number of chromosomes
Annoyingly, genotype arrays can have a number of chromosomes, or have chromosome names that the PLINK software does not like. One can easily write *ped, *fam and *map files without every opening PLINK. But try and get PLINK to read your files, and, the headaches will begin. So, what we do is, first count the number of chromosomes. Because all of our example dog files from Ancestry have the same structure, we'll use the genotypes from Logan (the dog).

```bash
awk '{print $1}' Logan.map | awk '!a[$0]++'  |wc -l > number_chromosomes.txt
``

This is just a command line parsing of the map file (which contains the chromosome names of each genotyped site in the genotype array, and counts the unique values. If you open number_chromosomes.txt, it will tell you there are 82 chromosomes in the dog genome. The version of the genome is CanFam3.1 and there are **NOT** 82 chromosomes in that genome. There are a bunch of chromosomes, and a bunch of other shorter scaffolds that aren't chromosome-level scaffolds.

#### Create binary version of plink files
All we do here is, in the shell, create an array called samples, add sample names to it, then use PLINK to stick, one by name, the sample names in a search for the ped and map files and thenm output binary versions. The important thing here is that we explcitly set the number of chromosoems to 82 and use the *--allow-extra-chr* switch to tell PLINK not to comlain about the weird number of chromosomes.

```bash
samples=()
for f in *.ped; do samples+=("${f%.ped}");done
for i in $samples; do plink  --ped ${i}.ped --map ${i}.map --make-bed --chr-set 82  --allow-extra-chr --out ${i}_binary;done
```

#### Create a list of dog sample binary PLINK files 
Why do we need to do this? Becaue in order to merge the individual sample PLINK files into one that contains the genotypes of all of the Ancestry dogs, we need to supply a file that lists these files.

```bash
ls *.ped | sed 's/.ped$//' | while read i; do
  awk -v val="$i" 'BEGIN { print val"_binary.bed "val"_binary.bim "val"_binary.fam" }' |grep -v Bambino>> merge_list.txt
done
```

`grep -v` means we are excluding the dog named Bambino. Why? Because of PLINK weirdness, we need to supply the files for a particular sample (dog) and then a list of the other samples you want to merge with it. 

#### Merge dog binaries into multi-sample PLINK files
Once again, we need to specify arguments so that PLINK doesn't complain about the non-standard number of chromosomes, merging all the files corresponding to samples in merge_list.txt to those for Bambino:

```bash 
plink --bfile Bambino_binary --chr-set 82 --merge-list merge_list.txt --make-bed --out dogs_merged --allow-extra-chr
``` 

#### Convert merged dog files to vcf format
[vcf](https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format) format is the standard format for representing genotypes calculated using genome sequencing data, whether that be whole-genome sequencing or some form of "reduced representation" , such as when sequencing is done for targed regions of a genome, e.g. protein-coding genes. It contains a number of header fields that are "commented out" with "#" characters, that typically describe the chromosomes names in the file, and what various codes mean. The data part of the file describes which variants are observed and in which samples at a given genomic position, e.g.:


<img src="img/vcf.png" width="100%" height="100%"/>


