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

I didn't make another direcotry, rather, I copied the following script into the existing directory that contained my raw hifi reads:
```
#!/bin/sh
#SBATCH --job-name fastqc_fasciatus
#SBATCH --nodes=1
#SBATCH --mem=40gb
#SBATCH --tasks-per-node=5 # Number of cores per node
#SBATCH --time=40:00:00
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

Check it out here: [_Plestiodon fasciatus_](fasciatus_fastqc.html)
