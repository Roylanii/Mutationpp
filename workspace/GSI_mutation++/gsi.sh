echo compile
g++ -g -I./ -I$MPP_DIRECTORY/install/include/mutation++ -c input.cpp -o build/input.o
g++ -g -I./ -I$MPP_DIRECTORY/install/include/mutation++ -c gsitest.cpp -o build/gsitest.o
echo link
g++ -g -L$MPP_DIRECTORY/install/lib build/gsitest.o build/input.o -lmutation++ -std=c++11 -o xx
echo exec
./xx --workdir /home/zhangjc/code/ablation/GSI_mutation++/data --mixture seb_oxidation_NASA9_ChemNonEq1T --pressure 101325 --temperature 3000
