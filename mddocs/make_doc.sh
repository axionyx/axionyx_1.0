source ~/anaconda2/etc/profile.d/conda.sh
conda activate py36
module load gcc
source ../../enable.sh
which python
doxygen doc.cfg
doxybook -i /tmp/temporary_doc_axionyx/xml -o .
