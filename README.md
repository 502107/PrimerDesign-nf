# PrimerDesign Nextflow Pipeline
[![Generic badge](https://img.shields.io/badge/Requires-nextflow-<COLOR>.svg)](https://www.nextflow.io/)
[![Generic badge](https://img.shields.io/badge/Requires-primer3-<COLOR>.svg)](https://anaconda.org/bioconda/primer3-py)
[![Generic badge](https://img.shields.io/badge/Requires-clustalw-<COLOR>.svg)](https://anaconda.org/bioconda/clustalw)
[![Generic badge](https://img.shields.io/badge/Requires-BLAST+-<COLOR>.svg)](https://anaconda.org/bioconda/blast)


- [Running the pipeline](#running-the-pipeline)  
- [Documentation](#documentation)  
- [Still to come](#still-to-come)  

## Running the pipeline

When working on the local machine, just run:
```
./main.nf -profile local
```
> Or, `./main.nf -profile local -resume`

When working on the HPC, run:
```
./main.nf -profile slurm
```

> [!NOTE]
> - You can install the dependencies (i.e. conda/mamba and environment), as well as download missing reference genomes for `assets/refs/` by `./nf-install.sh`
> - Currently we check either the forward or reverse primer, by changing the direction variable of the _PrimerDesign_ process inside `main.nf`, to either _forward_ or _reverse_.

## Documentation

The pipeline currently works on Pgt210 and Pst104e annotated genes. In a future implementation this will be updated to also consider PgtCRL75, although not actively working to update this.
There are 4 main processes in this workflow, outlined below.
  > The pipeline could be made faster by only considering one of the rusts, since nextflow is well-suited for relay processes.
  > This is not an option in the current implementation, because of post-sorting of the outputs inside succeeding processes.
  > Essentially, we need to wait for each process to finish before moving on to the next one.

### Extracting Genes
> Process: ExtractValidGenes

The workflow begins by considering the gene names inside p[g/s]t_genes.txt.

1. Searching for the gene:
    - The location of the genes in the fasta file of the corresponding genome, is found from the annotation (.gff) file.
    - The reference isolate is sorted and the gene extracted from the location specified in the annotation.

2. Padding the gene:
   - Add 500 bp padding regardless of gene size.

3. Verifying diversity:
   - We do a quick search to see if the sequence is repeated anywhere else. If yes, we add 500 bp padding to the already padded gene, and search again.
  
### BLAST
> Process: BlastSearch

Here, we do a simple BLAST search of the extracted genes against all the reference isolates[^1].

### Find Conserved Sites
> Process: FindConservedSites

Using the BLAST outputs from the previous process, here we filter the results for each of the reference isolates,
discounting insignificant hits, using the following conditions:
1. If pident (% of identical positions) is 100%, or
2. If pident >= 93 and the length of the sequence >= 1000, then
3. Check if there are any other hits from the BLAST search under the same sseqid, if there are and if they come from the same gene, then concatenate them in one sequence, with hyphens between the gap of the two hits.
4. Check that all the blast hits have at most 500b size difference from the original - if there's more then drop the hit. If all the blast results have >500b difference in size, then just consider everything, but most likely the primer search will fail.

Following the filtering, the regions that were present in the BLAST results of all the isolates are concatinated in a fasta file, with SNPs in the sequence replaced with hyphens.

### Designing Primers
> Process: DesignPrimers

Using the conserved fasta files created in the previous step, here we use Primer3 to generate 10,000 primer pairs for each of the genes.
> The optimal conditions are pretty much the same as the defaults[^2], with certain changes, such as the length of the product being in the range of 1,000 - 3,000 bp.

We need each of our primer pairs to only amplify one gene, and work across all the reference isolates. We also need to ensure that the primer pair doesn't amplify any region in the host genome.
In order to filter out these primers and find the best hits for each of the genes, the following steps are applied:
1. First, the primer pairs are sorted with respect to how close they are to the actual gene sequence (i.e. without padding). At the moment, we are not considering any primers that go into the sequence.
2. In order to check for regions the _forward_ primer and its reverse complement may align with, we search the length of the combined reference isolates with a moving window, for any region that is similar to the primer by at least 4 mismatches (i.e. if a region has 4 SNPs or less, compared to the primer sequence). This primer would therefore be disregarded, and the next closest to the actual sequence checked. Here, we also check for regions that have 4 SNPs or less, if at least one SNP exists in the last 3 base pairs (or first three for the reverse complement of the right primer) -- if they do we then accept that hit.
3. Obviously, this condition would always return a hit, given we have the gene sequence in the fasta files of the reference isolates. So, here we also check the blast results again, filtered as discussed in the previous process, and we obtain the sseqid (subject sequence ID -- i.e. reference isolate scaffold ID), start and end positions of the hits. These regions are replaced by ambiguities (N), so during the mismatch search these regions are not checked.
4. The same principle of checking for mismatches, is also applied for the host genome.
5. The valid genes are thus saved to csv files.

## Still to come
- [X] Some of the genes we're interested in, go above the 3 kbp threshold we have inexplicitly applied here. In a future update we need to handle genes whose lengths are >3.5 kbp with multiple overlapping primer pairs.
- [ ] Add option to handle Pgt genes from CRL75.
- [X] Streamline installation for users of Mac and Linux.
- [X] `assets/` directory, with the different isolates is not included due to file size restrictions. Working on a cloud link to download outside of github. This directory looks like this:
> assets/  
> └── refs  
>     ├── host  
>     │   └── GCF_018294505.1.fna  
>     ├── pgt210  
>     │   ├── GCA_008522505.1_Pgt_210_genomic.fna  
>     │   └── pgt210_genomic.gff  
>     ├── pgt_all  
>     │   ├── GCA_000149925.1_ASM14992v1_genomic.fna  
>     │   ├── GCA_002762355.2_KSU_Pgt_99KS76A_2.0_genomic.fna  
>     │   ├── GCA_008520325.1_Pgt_Ug99_genomic.fna  
>     │   ├── GCA_008522505.1_Pgt_210_genomic.fna  
>     │   └── GCA_903797515.1_Puccinia_graminis_f._sp._tritici_UK-01_genomic.fna  
>     ├── pgtcrl  
>     │   ├── GCA_000149925.1_ASM14992v1_genomic.fna  
>     │   └── crl75_genomic.gff  
>     ├── pst104  
>     │   ├── pst104e.fna  
>     │   └── pst104e.gff  
>     └── pst_all  
>         ├── GCA_011750755.1_ASM1175075v1_genomic.fna  
>         ├── GCA_021901695.1_Pst134E36_v1_pri_genomic.fna  
>         ├── GCA_025169535.1_ASM2516953v1_genomic.fna  
>         ├── GCA_025169555.1_ASM2516955v1_genomic.fna  
>         ├── GCA_025869495.1_ASM2586949v1_genomic.fna  
>         └── pst104e.fna  


[^1]: Pgt reference isolates: GCA_903797515.1_Puccinia_graminis_f._sp._tritici_UK-01, GCA_008522505.1_Pgt_210, GCA_008520325.1_Pgt_Ug99, GCA_002762355.2_KSU_Pgt_99KS76A_2.0, GCA_000149925.1_ASM14992v1.  
Pst reference isolates: GCA_025869495.1_ASM2586949v1, GCA_025169555.1_ASM2516955v1, GCA_011750755.1_ASM1175075v1, GCA_025169535.1_ASM2516953v1, GCA_021901695.1_Pst134E36_v1_pri, pst104e.
[^2]: [Primer3 Manual](https://primer3.org/manual.html)
