#!/bin/bash

pckg="biopython nextflow python=3.11 blast clustalw primer3-py"

if command -v mamba &> /dev/null; then
	if mamba env list | grep -q nextflow-env; then
                :
        else
                mamba create -n nextflow-env -y -c bioconda -c conda-forge $pckg
        fi
elif [[ "$OSTYPE" == "linux-gnu"* || "$OSTYPE" == "darwin"* ]]; then
        echo "$OSTYPE"
        if [[ "$OSTYPE" == "linux-gnu"* ]]; then
                sudo apt-get install bzip2
        elif [[ "$OSTYPE" == "darwin"* ]]; then
                /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
                brew install bzip2
        fi
        wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
        ./bin/micromamba shell init -s bash -p ~/micromamba
        source ~/.bashrc
        sleep 1
        micromamba create -n base -c conda-forge -y
        micromamba activate base
        micromamba install mamba -c conda-forge -y
        mamba create -n nextflow-env -y -c bioconda -c conda-forge $pckg
        mamba init
else
        echo "Unsupported operating system"
        exit 1
fi

# Need to update shebangs in main.nf and python scripts

nf_path="#!$HOME/mambaforge/envs/nextflow/bin/nextflow"
python_path="#!$HOME/mambaforge/envs/nextflow/bin/python"

sed -i-tmp "1s|#.*|${nf_path}|" main.nf 
rm main.nf-tmp

sed -i-tmp "1s|#.*|${python_path}|" bin/*py

# Download the contents of the assets directory

mkdir -p assets/refs/host
mkdir -p assets/refs/pgt210
mkdir -p assets/refs/pst104
mkdir -p assets/refs/pgt_all
mkdir -p assets/refs/pst_all

# Host reference
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/294/505/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna.gz -P assets/refs/host/
# Pst references
# For pst104e need to login here: https://genome.jgi.doe.gov/portal/Pucstr1/Pucstr1.download.html; so i'll include it using git lfs
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/021/901/695/GCA_021901695.1_Pst134E36_v1_pri/GCA_021901695.1_Pst134E36_v1_pri_genomic.fna.gz -P assets/refs/pst_all/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/025/169/535/GCA_025169535.1_ASM2516953v1/GCA_025169535.1_ASM2516953v1_genomic.fna.gz -P assets/refs/pst_all/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/750/755/GCA_011750755.1_ASM1175075v1/GCA_011750755.1_ASM1175075v1_genomic.fna.gz -P assets/refs/pst_all/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/025/169/555/GCA_025169555.1_ASM2516955v1/GCA_025169555.1_ASM2516955v1_genomic.fna.gz -P assets/refs/pst_all/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/025/869/495/GCA_025869495.1_ASM2586949v1/GCA_025869495.1_ASM2586949v1_genomic.fna.gz -P assets/refs/pst_all/
# Pgt references
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/522/505/GCA_008522505.1_Pgt_210/GCA_008522505.1_Pgt_210_genomic.fna.gz -P assets/refs/pgt210/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/522/505/GCA_008522505.1_Pgt_210/GCA_008522505.1_Pgt_210_genomic.gff.gz -P assets/refs/pgt210/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/903/797/515/GCA_903797515.1_Puccinia_graminis_f._sp._tritici_UK-01/GCA_903797515.1_Puccinia_graminis_f._sp._tritici_UK-01_genomic.fna.gz -P assets/refs/pgt_all/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/762/355/GCA_002762355.2_KSU_Pgt_99KS76A_2.0/GCA_002762355.2_KSU_Pgt_99KS76A_2.0_genomic.fna.gz -P assets/refs/pgt_all/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/008/520/325/GCA_008520325.1_Pgt_Ug99/GCA_008520325.1_Pgt_Ug99_genomic.fna.gz -P assets/refs/pgt_all/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/149/925/GCA_000149925.1_ASM14992v1/GCA_000149925.1_ASM14992v1_genomic.fna.gz -P assets/refs/pgt_all/

gunzip assets/refs/host/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna.gz
mv assets/refs/host/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna assets/refs/host/GCF_018294505.1.fna

gunzip assets/refs/pgt210/GCA_008522505.1_Pgt_210_genomic.fna.gz
gunzip assets/refs/pst104/pst104e.fna.gz
cp assets/refs/pgt210/GCA_008522505.1_Pgt_210_genomic.fna assets/refs/pgt_all/.
cp assets/refs/pst104/pst104e.fna assets/refs/pst_all/.

gunzip assets/refs/pst_all/GCA_021901695.1_Pst134E36_v1_pri_genomic.fna.gz
gunzip assets/refs/pst_all/GCA_025169535.1_ASM2516953v1_genomic.fna.gz
gunzip assets/refs/pst_all/GCA_011750755.1_ASM1175075v1_genomic.fna.gz
gunzip assets/refs/pst_all/GCA_025169555.1_ASM2516955v1_genomic.fna.gz
gunzip assets/refs/pst_all/GCA_025869495.1_ASM2586949v1_genomic.fna.gz

gunzip assets/refs/pgt_all/GCA_903797515.1_Puccinia_graminis_f._sp._tritici_UK-01_genomic.fna.gz
gunzip assets/refs/pgt_all/GCA_002762355.2_KSU_Pgt_99KS76A_2.0_genomic.fna.gz
gunzip assets/refs/pgt_all/GCA_008520325.1_Pgt_Ug99_genomic.fna.gz
gunzip assets/refs/pgt_all/GCA_000149925.1_ASM14992v1_genomic.fna.gz
