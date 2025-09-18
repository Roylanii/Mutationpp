#mutationpp.bashrc
cd $( dirname $BASH_SOURCE)
export MPP_DIRECTORY="$( pwd )"
cd -
#export MPP_DIRECTORY=/home/zhangjc/WORKDIR/Mutationpp
export MPP_DATA_DIRECTORY=$MPP_DIRECTORY/data
export PATH=$MPP_DIRECTORY/install/bin:$PATH
export LD_LIBRARY_PATH=$MPP_DIRECTORY/install/lib:$LD_LIBRARY_PATH
