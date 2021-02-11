export PYTHONPATH=$PYTHONPATH:$(readlink -f ../midas/MIDAS)
export PATH=$PATH:$(readlink -f ../midas/MIDAS)/scripts
export MIDAS_DB=$(readlink -f ../midas/midas_db_v1.2)
