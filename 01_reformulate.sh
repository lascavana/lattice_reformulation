# compile code
rm -r build
mkdir build
cd build
cmake .. -DSCIP_DIR=$SCIPOPTDIR
make
cd ..

# for i in {1..50}
# do
#   ./build/reformulate benchmarks/marksplit/instance_${i}.lp settingsfile.set
# done

# for i in {1..5}
# do
#   ./build/reformulate benchmarks/cuww/instance_${i}.lp settingsfile.set
# done

# for i in {1..10}
# do
#   ./build/reformulate benchmarks/prob_struct/prob${i}.lp settingsfile.set
# done

# for i in {11..20}
# do
#   ./build/reformulate benchmarks/prob_nostruct/prob${i}.lp settingsfile.set
# done

# for i in {1..50}
# do
#   ./build/reformulate benchmarks/mknapsack/instance_${i}_100_6.lp settingsfile.set
# done

# for i in {1..50}
# do
#   ./build/reformulate benchmarks/msplit_feasibility/instance_${i}.lp settingsfile.set
# done

./build/reformulate benchmarks/cuww/instance_1.lp kc settingsfile.set