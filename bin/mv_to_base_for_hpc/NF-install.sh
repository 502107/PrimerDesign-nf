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

## Need to update shebangs in main.nf and python scripts

