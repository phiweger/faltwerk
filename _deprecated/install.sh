

# !wget https://mmseqs.com/foldseek/foldseek-osx-universal.tar.gz; tar xvzf foldseek-osx-universal.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH
# !foldseek
# foldseek Version: 7d0c07f89a8f254789342dc041726cfc1c7223e1

ENV_NAME=faltwerk
SHELL_NAME=bash

conda init $SHELL_NAME

conda env create -y -n $ENV_NAME --file install/requirements.yml
conda activate $NAME

pip install -r install/more_requirements.txt

# On my laptop this is not necessary, but it might be:
# > In order to use with JupyterLab you must install the JupyterLab
# extension -- https://pypi.org/project/py3Dmol/
# jupyter labextension install jupyterlab_3dmol

pip install -e .

conda deactivate
conda clean -a



