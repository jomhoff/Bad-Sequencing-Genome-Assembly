## **Background**

Ok. 
So,
this workflow is made to check the read quality, assemble the genome, and check the quality of the assembly for a _Plestiodon fasciatus_ genome. Although this genome was sequenced with PacBio3 Revio sequencing, the results we very poor, likely due to preservation in NAP Buffer.
Here are the stats:
- Yield (Gb): 4.98
- Read Length (kb): 8.71

Check the .README for details regarding extraction QC, library QC, and read length distribution.

## **Raw Read Quality Assessment with FastQC -- Adapted from [Amanda Markee](https://github.com/amandamarkee/actias-luna-genome.git)**

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a quality assessment tool for next generation sequencing, often used to assess raw reed quality and highlight problem areas in the form of visualizations (see results below).

Ideally, you would create a FastQC directory for the output data and script as so:
```
mkdir fastqc_raw
```

I didn't make another directory, rather, I copied the following script into the existing directory that contained my raw hifi reads:
```
#!/bin/sh
#SBATCH --job-name fastqc_fasciatus
#SBATCH --nodes=1
#SBATCH --mem=100gb
#SBATCH --tasks-per-node=20 # Number of cores per node
#SBATCH --time=100:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jhoffman1@amnh.org
#SBATCH --output=slurm-%j-%x.out
#conda init
source ~/.bash_profile
conda activate fasciatus_ass

fastqc m84082_240409_161410_s4.hifi_reads.bc2049.fastq
~
```

The FastQC output files include an .html file, which will contain visualizations for the following results:
- Basic Statistics
- Per base sequence quality
- Per sequence quality scores (Phred scores)
- Per base sequence content
- Per sequence GC content
- Per base N content
- Sequence Length Distribution
- Sequence Duplication Levels
- Overrepresented sequences
- Adapter content

Check my FastQC results here: [_Plestiodon fasciatus_](fasciatus_fastqc.html)


## **09/19/2022; Genome Assembly with hifiasm -- Once Again Adapted from THE [Amanda Markee](https://github.com/amandamarkee/actias-luna-genome.git)**

[hifiasm](https://hifiasm.readthedocs.io/en/latest/) is a fast and easy haplotype-resolved de novo assembly software for PacBio HiFi reads
 - hifiasm documentation explaining input parameters: https://hifiasm.readthedocs.io/en/latest/pa-assembly.html
 - hifiasm documentation explaining output files: https://hifiasm.readthedocs.io/en/latest/interpreting-output.html

Originally, this source script for hifiasm utilized the aggressive duplicate purging options in Hifiasm (option -l 2). By default, hifiasm purges haplotig duplications. Normally, this would be a good approach to take. For some cases, such as with inbred or homozygous genomes, it is useful to specify 0 haplotig duplications. Since this genome is already struggling for size and completeness, I decided to assemble it twice, once with fairly aggressive purging (-l 2) and one with no purging (-l 0) to see if I can salvage more data.

With purging:
```
#!/bin/sh
#SBATCH --job-name hoff_hifiasm
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --tasks-per-node=1 # Number of cores per node
#SBATCH --time=30:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jhoffman1@amnh.org
#SBATCH --output=slurm-%j-%x.outt
#conda init
source ~/.bash_profile
conda activate fasciatus_ass

hifiasm -o hoff_hifi_assembly.asm -l 2 -t 32 /home/jhoffman1/mendel-nas1/fasciatus_genome/6354_import-dataset/hifi_reads/m84082_240409_161410_s4.hifi_reads.bc2049.fastq

```
Without purging:
```
#!/bin/sh
#SBATCH --job-name hoff_hifiasm
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --tasks-per-node=1 # Number of cores per node
#SBATCH --time=30:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jhoffman1@amnh.org
#SBATCH --output=slurm-%j-%x.out
#conda init
source ~/.bash_profile
conda activate fasciatus_ass

hifiasm -o hoff_hifi_assembly.asm -l 0 -t 32 /home/jhoffman1/mendel-nas1/fasciatus_genome/6354_import-dataset/hifi_reads/m84082_240409_161410_s4.hifi_reads.bc2049.fastq

```
<br />

## **10/03/2022; Genome Assembly Quality Assessment with assemblystats.py -- Again, Adapted From the ~~Late~~ Great [Amanda Markee](https://github.com/amandamarkee/actias-luna-genome.git)**

- After assembly with hifiasm, we can assess assembly quality using the [assemblystats.py script](https://github.com/MikeTrizna/assembly_stats/tree/0.1.4) created by Mike Trizna.
- The version of assemblystats.py used here was modified by Paul Frandsen (Brigham Young University).

First, I copied this script into my working directory, and called it assemblystats.py

```
#!/usr/bin/env python

import numpy as np
from itertools import groupby
import json
import sys


def fasta_iter(fasta_file):
    """Takes a FASTA file, and produces a generator of Header and Sequences.
    This is a memory-efficient way of analyzing a FASTA files -- without
    reading the entire file into memory.

    Parameters
    ----------
    fasta_file : str
        The file location of the FASTA file

    Returns
    -------
    header: str
        The string contained in the header portion of the sequence record
        (everything after the '>')
    seq: str
        The sequence portion of the sequence record
    """

    fh = open(fasta_file)
    fa_iter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in fa_iter:
        # drop the ">"
        header = next(header)[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.upper().strip() for s in next(fa_iter))
        yield header, seq


def read_genome(fasta_file):
    """Takes a FASTA file, and produces 2 lists of sequence lengths. It also
    calculates the GC Content, since this is the only statistic that is not
    calculated based on sequence lengths.

    Parameters
    ----------
    fasta_file : str
        The file location of the FASTA file

    Returns
    -------
    contig_lens: list
        A list of lengths of all contigs in the genome.
    scaffold_lens: list
        A list of lengths of all scaffolds in the genome.
    gc_cont: float
        The percentage of total basepairs in the genome that are either G or C.
    """

    gc = 0
    total_len = 0
    contig_lens = []
    scaffold_lens = []
    for _, seq in fasta_iter(fasta_file):
        scaffold_lens.append(len(seq))
        if "NN" in seq:
            contig_list = seq.split("NN")
        else:
            contig_list = [seq]
        for contig in contig_list:
            if len(contig):
                gc += contig.count('G') + contig.count('C')
                total_len += len(contig)
                contig_lens.append(len(contig))
    gc_cont = (gc / total_len) * 100
    return contig_lens, scaffold_lens, gc_cont


def calculate_stats(seq_lens, gc_cont):
    stats = {}
    seq_array = np.array(seq_lens)
    stats['sequence_count'] = seq_array.size
    stats['gc_content'] = gc_cont
    sorted_lens = seq_array[np.argsort(-seq_array)]
    stats['longest'] = int(sorted_lens[0])
    stats['shortest'] = int(sorted_lens[-1])
    stats['median'] = np.median(sorted_lens)
    stats['mean'] = np.mean(sorted_lens)
    stats['total_bps'] = int(np.sum(sorted_lens))
    csum = np.cumsum(sorted_lens)
    for level in [10, 20, 30, 40, 50]:
        nx = int(stats['total_bps'] * (level / 100))
        csumn = min(csum[csum >= nx])
        l_level = int(np.where(csum == csumn)[0])
        n_level = int(sorted_lens[l_level])

        stats['L' + str(level)] = l_level
        stats['N' + str(level)] = n_level
    return stats


if __name__ == "__main__":
    infilename = sys.argv[1]
    contig_lens, scaffold_lens, gc_cont = read_genome(infilename)
    contig_stats = calculate_stats(contig_lens, gc_cont)
    scaffold_stats = calculate_stats(scaffold_lens, gc_cont)
    stat_output = {'Contig Stats': contig_stats,
                   'Scaffold Stats': scaffold_stats}
    print(json.dumps(stat_output, indent=2, sort_keys=True))
```

Next, I changed permissions as follows to allow execution permissions.
```
chmod +x assemblystats.py
```

Then, I produced a FASTA file from the initial GFA output files from the hifiasm assembly output for haplotype 1, haplotype 2, and total. I used the primary contig file, as indicated with asm.bp.p_ctg.fa (p_ctg meaning primary contig)
```
awk '/^S/{print ">"$2;print $3}' hoff_hifi_assembly.asm.bp.hap1.p_ctg.gfa > hoff_hifi_assembly.hap1-l0.ctg.fa
awk '/^S/{print ">"$2;print $3}' hoff_hifi_assembly.asm.bp.hap2.p_ctg.gfa > hoff_hifi_assembly.hap2-l0.ctg.fa
awk '/^S/{print ">"$2;print $3}' hoff_hifi_assembly.asm.bp.p_ctg.gfa > hoff_hifi_assembly.total-l0.ctg.fa
```

Lastly, I ran the assemblystats.py script on the newly generated fasta file of the fga in the format of scriptfilepath/scirptname.py nameofassembly.fa and save as a txt file 
```
./assemblystats.py ./hoff_hifi_assembly.hap1-l0.ctg.fa >> hoff_hifi_ass.stats-h1l0.txt
./assemblystats.py ./hoff_hifi_assembly.hap2-l0.ctg.fa >> hoff_hifi_ass.stats-h2l0.txt
./assemblystats.py ./hoff_hifi_assembly.total-l0.ctg.fa >> hoff_hifi_ass.stats-Tl0.txt
```
Save results as a text file as shown.
```
./assemblystats.py aclu_hifi_assembly_06-14-2022.asm.bp.p_ctg.fa >> aclu_hifi_assembly_06-14-2022.txt
```
Let's compare the total assembly contig statistics for the two assemblies:

With Aggressive Purging:
  > Contig Stats: 
  >  L10: 866,
  >  L20: 2135,
  >  L30: 3714,
  >  L40: 5589,
  >  L50: 7780,
  >  N10: 87238,
  >  N20: 67486,
  >  N30: 55699,
  >  N40: 47290,
  >  N50: 40767,
  >  gc_content: 44.955823538395784,
  >  longest: 286376,
  >  mean: 36067.164390756305,
  >  median: 29920.0,
  >  sequence_count: 26656,
  >  shortest: 8950,
  >  total_bps: 961406334

With No Purging:
  > Contig Stats: 
  >  L10: 919,
  >  L20: 2274,
  >  L30: 3967,
  >  L40: 5977,
  >  L50: 8326,
  >  N10: 83702,
  >  N20: 64740,
  >  N30: 53535,
  >  N40: 45424,
  >  N50: 39086,
  >  gc_content: 44.93135646459214,
  >  longest: 286376,
  >  mean: 34883.374845738865,
  >  median: 28997.0,
  >  sequence_count: 28361,
  >  shortest: 7839,
  >  total_bps: 989327394

As we can see, the No Purging run resulted in a slightly better N50 value. The rest of the statistics are similar, with the No Purging run having slightly shorter contigs including that likely were purged in the other run.

Moving forwards, I will be using the No Purging assembly. 

## **Checking Assembly Completeness with BUSCO**

[BUSCO](https://busco.ezlab.org/busco_userguide.html) is a program that estimates genome completeness based on evolutionarily-informed expectations of gene content of near-universal single-copy orthologs.

Since the worker nodes of the AMNH's computational clusters don't have access to the internet, it is necesarry to install BUSCO locally:
```
#clone repository
git clone https://gitlab.com/ezlab/busco.git
cd busco/

#pip install
python -m pip install .
```

Make sure you have all the [dependencies](https://busco.ezlab.org/busco_userguide.html#editing-busco-run-configuration) installed for the type of BUSCO run you are planning on running.

Then, I ran this shell file:
```
#!/bin/sh
#SBATCH --job-name Busco_Genomepfas
#SBATCH --nodes=20
#SBATCH --mem=100gb
#SBATCH --time=144:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jhoffman1@amnh.org

source ~/.bash_profile
conda activate fasciatus_ass
pfas="/home/jhoffman1/mendel-nas1/fasciatus_genome/6354_import-dataset/hifi_reads/hoff_hifi_assembly.total-l0.ctg.fa"
busco -m genome -i $pfas -o pfasBUSCO_saur -l sauropsida_odb10 -f --metaeuk --offline --download_path /home/jhoffman1/mendel-nas1/fasciatus_genome
```


## **Ragtag onto a Chromosome Level Genome**

In order to make some use out of this long-read genome, let's map it onto the closest chromosome-level reference genome (Tiliqua scincoides), a mere ~85 million yesars divereged. 

For this, we will be using [Ragtag](https://github.com/malonge/RagTag), which is a collection of software tools for scaffolding and improving modern genome assemblies.

First, let's copy the assembled genome to our Ragtag directory
```
cp ./fasciatus_genome/6354_import-dataset/hifi_reads/hoff_hifi_assembly.total-l0.ctg.fa ./ragtag
```

ok cool, now we download the Tiliqua genome into the same folder. For me, I uploaded it from local using [Cyberduck](https://cyberduck.io/).

Install ragtag
```
#create conda environment
conda create -n ragtag
conda activate ragtag

#install ragtag
conda install -c bioconda ragtag
```

Shell script for running ragtag
```
#!/bin/sh
#SBATCH --job-name ragtag
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --time=500:00:00
#SBATCH --mail-type=ALL
#SBATCH --tasks-per-node=40
#SBATCH --mail-user=jhoffman1@amnh.org
#SBATCH --output=slurm-%j-%x.out

source ~/.bash_profile
conda activate ragtag
 ragtag.py scaffold GCA_035046505.1_rTilSci1.hap2_genomic.fna hoff_hifi_assembly.total-l0.ctg.fa -r -g 2 -m 100000000 -o pseudochrom_pfas/
```







