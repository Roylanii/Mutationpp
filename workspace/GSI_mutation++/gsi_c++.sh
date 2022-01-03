echo compile
g++ -g -I./ -I$MPP_DIRECTORY/install/include/mutation++ -c src/sub/input.cpp -o obj_c/input.o
g++ -g -I./src/sub -I$MPP_DIRECTORY/install/include/mutation++ -c src/main/GSI_mutation++.cpp -o obj_c/GSI_mutation++.o
echo link
g++ -g -L$MPP_DIRECTORY/install/lib obj_c/GSI_mutation++.o obj_c/input.o -lmutation++ -std=c++11 -o gsi_c++
echo obj_c_exect
#./xx --workdir /home/zhangjc/code/ablation/GSI_mutation++/data --mixture seb_oxidation_NASA9_ChemNonEq1T --pressure 101325 --temperature 3000
