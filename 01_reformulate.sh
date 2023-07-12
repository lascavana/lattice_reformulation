# compile code
rm -r build
mkdir build
cd build
cmake .. -DSCIP_DIR=$SCIPOPTDIR
make
cd ..

# struct_s
for i in {1..30}
do
  ./build/reformulate benchmarks/struct_s/instance_${i}_n10_struct.lp ahl settingsfile.set
  ./build/reformulate benchmarks/struct_s/instance_${i}_n10_struct.lp ahl_diag settingsfile.set
  ./build/reformulate benchmarks/struct_s/instance_${i}_n10_struct.lp ahl_poor settingsfile.set
  ./build/reformulate benchmarks/struct_s/instance_${i}_n10_struct.lp kp settingsfile.set
done

# struct_b
for i in {1..30}
do
  ./build/reformulate benchmarks/struct_b/instance_${i}_n100_struct.lp ahl settingsfile.set
  ./build/reformulate benchmarks/struct_b/instance_${i}_n100_struct.lp ahl_diag settingsfile.set
  ./build/reformulate benchmarks/struct_b/instance_${i}_n100_struct.lp ahl_poor settingsfile.set
  ./build/reformulate benchmarks/struct_b/instance_${i}_n100_struct.lp kp settingsfile.set
done

# nostruct_s
for i in {1..30}
do
  ./build/reformulate benchmarks/nostruct_s/instance_${i}_n10_nostruct.lp ahl settingsfile.set
  ./build/reformulate benchmarks/nostruct_s/instance_${i}_n10_nostruct.lp ahl_diag settingsfile.set
  ./build/reformulate benchmarks/nostruct_s/instance_${i}_n10_nostruct.lp ahl_poor settingsfile.set
  ./build/reformulate benchmarks/nostruct_s/instance_${i}_n10_nostruct.lp kp settingsfile.set
done

# nostruct_b
for i in {1..30}
do
  ./build/reformulate benchmarks/nostruct_s/instance_${i}_n100_nostruct.lp ahl settingsfile.set
  ./build/reformulate benchmarks/nostruct_s/instance_${i}_n100_nostruct.lp ahl_diag settingsfile.set
  ./build/reformulate benchmarks/nostruct_s/instance_${i}_n100_nostruct.lp ahl_poor settingsfile.set
  ./build/reformulate benchmarks/nostruct_s/instance_${i}_n100_nostruct.lp kp settingsfile.set
done

# Market split
for i in {1..30}
do
  ./build/reformulate benchmarks/MS/instance_${i}.lp ahl settingsfile.set
  ./build/reformulate benchmarks/MS/instance_${i}.lp ahl_poor settingsfile.set
  ./build/reformulate benchmarks/MS/instance_${i}.lp kp settingsfile.set
done

# Generalized assignment 
for i in {1..30}
do
  ./build/reformulate benchmarks/GAP/instance_${i}_100_6.lp ahl settingsfile.set
  ./build/reformulate benchmarks/GAP/instance_${i}_100_6.lp ahl_poor settingsfile.set
  ./build/reformulate benchmarks/GAP/instance_${i}_100_6.lp kp settingsfile.set
done

# Combinatorial auctions
for i in {1..30}
do
  ./build/reformulate benchmarks/CA/instance_${i}.lp ahl settingsfile.set
  ./build/reformulate benchmarks/CA/instance_${i}.lp ahl_poor settingsfile.set
  ./build/reformulate benchmarks/CA/instance_${i}.lp kp settingsfile.set
done



