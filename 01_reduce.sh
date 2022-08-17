# compile code
rm -r build
mkdir build
cd build
cmake .. -DSCIP_DIR=/Users/lvscavuzzomont/opt/scipoptsuite-8.0.0
make
cd ..

# for i in {1..50}
# do
#   ./build/reformulate benchmarks/marksplit/instance_${i}.lp settingsfile.set
# done

for i in {1..5}
do
  ./build/reformulate benchmarks/cuww/instance_${i}.lp settingsfile.set
done

# for i in {1..5}
# do
#   ./build/reformulate benchmarks/prob/prob${i}.lp settingsfile.set
# done

# for i in {1..50}
# do
#   ./build/reformulate benchmarks/mknapsack/instance_${i}_100_6.lp settingsfile.set
# done
