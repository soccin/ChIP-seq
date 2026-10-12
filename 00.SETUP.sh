#
# Do not try to run needs to be sourced, from the ChIP-seq directory
#   . 00.SETUP.sh
#
# Builds the venv for MACS2. MACS2 2.2.9.1 builds under the python 3.10
# on the IRIS PATH with numpy pinned to 1.23.0 (the recipe ATAC-seq uses
# on IRIS).
#

python3 -m venv venv

. venv/bin/activate
pip install --upgrade pip
pip install numpy==1.23.0
pip install MACS2==2.2.9.1

deactivate
