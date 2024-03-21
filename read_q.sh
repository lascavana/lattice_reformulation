# compile code
rm -r build
mkdir build
cd build
cmake .. -DSCIP_DIR=$SCIPOPTDIR
make
cd ..

# mkdir Q
# for i in {1..30}
# do
#   ./build/reformulate benchmarks/nostruct_s/instance_${i}.lp ahl settingsfile.set
#   mv Q.txt Q/Q${i}.txt
# done

mkdir Q
mkdir Aext
declare -a arr=("ej.mps" "neos859080.mps" "p2m2p1m1p0n100.mps" "ponderthis0517-inf.mps" "pb-market-split8-70-4.mps" "gt2.mps" "p0201.mps" "stein45inf.mps" "stein15inf.mps")
for i in "${arr[@]}"
do
    ./build/reformulate benchmarks/MIPLIB/$i ahl settingsfile.set
    mv Aext.txt Aext/${i}.txt
    mv Q.txt Q/${i}.txt
done