# compile code
rm -r build
mkdir build
cd build
cmake .. -DSCIP_DIR=$SCIPOPTDIR
make
cd ..

rm -r logs
mkdir logs

# struct_s
for i in {1..30}
do
  ./build/reformulate benchmarks/struct_s/instance_${i}_n10_struct.lp ahl settingsfile.set >> logs/struct_s.txt
  ./build/reformulate benchmarks/struct_s/instance_${i}_n10_struct.lp ahl_diag settingsfile.set >> logs/struct_s.txt
  ./build/reformulate benchmarks/struct_s/instance_${i}_n10_struct.lp ahl_poor settingsfile.set >> logs/struct_s.txt
  ./build/reformulate benchmarks/struct_s/instance_${i}_n10_struct.lp kp settingsfile.set >> logs/struct_s.txt
done

# struct_b
for i in {1..30}
do
  ./build/reformulate benchmarks/struct_b/instance_${i}_n100_struct.lp ahl settingsfile.set >> logs/struct_b.txt
  ./build/reformulate benchmarks/struct_b/instance_${i}_n100_struct.lp ahl_diag settingsfile.set >> logs/struct_b.txt
  ./build/reformulate benchmarks/struct_b/instance_${i}_n100_struct.lp ahl_poor settingsfile.set >> logs/struct_b.txt
  ./build/reformulate benchmarks/struct_b/instance_${i}_n100_struct.lp kp settingsfile.set >> logs/struct_b.txt
done

# nostruct_s
for i in {1..30}
do
  ./build/reformulate benchmarks/nostruct_s/instance_${i}_n10_nostruct.lp ahl settingsfile.set >> logs/nostruct_s.txt
  ./build/reformulate benchmarks/nostruct_s/instance_${i}_n10_nostruct.lp ahl_diag settingsfile.set >> logs/nostruct_s.txt
  ./build/reformulate benchmarks/nostruct_s/instance_${i}_n10_nostruct.lp ahl_poor settingsfile.set >> logs/nostruct_s.txt
  ./build/reformulate benchmarks/nostruct_s/instance_${i}_n10_nostruct.lp kp settingsfile.set >> logs/nostruct_s.txt
done

# nostruct_b
for i in {1..30}
do
  ./build/reformulate benchmarks/nostruct_s/instance_${i}_n100_nostruct.lp ahl settingsfile.set >> logs/nostruct_b.txt
  ./build/reformulate benchmarks/nostruct_s/instance_${i}_n100_nostruct.lp ahl_diag settingsfile.set >> logs/nostruct_b.txt
  ./build/reformulate benchmarks/nostruct_s/instance_${i}_n100_nostruct.lp ahl_poor settingsfile.set >> logs/nostruct_b.txt
  ./build/reformulate benchmarks/nostruct_s/instance_${i}_n100_nostruct.lp kp settingsfile.set >> logs/nostruct_b.txt
done

# Market split
for i in {1..30}
do
  ./build/reformulate benchmarks/MS/instance_${i}.lp ahl settingsfile.set >> logs/ms.txt
  ./build/reformulate benchmarks/MS/instance_${i}.lp ahl_poor settingsfile.set >> logs/ms.txt
  ./build/reformulate benchmarks/MS/instance_${i}.lp kp settingsfile.set >> logs/ms.txt
done

# Generalized assignment 
for i in {1..30}
do
  ./build/reformulate benchmarks/GAP/instance_${i}_100_6.lp ahl settingsfile.set >> logs/gap.txt
  ./build/reformulate benchmarks/GAP/instance_${i}_100_6.lp ahl_poor settingsfile.set >> logs/gap.txt
  ./build/reformulate benchmarks/GAP/instance_${i}_100_6.lp kp settingsfile.set >> logs/gap.txt
done

# Combinatorial auctions
for i in {1..30}
do
  ./build/reformulate benchmarks/CA/instance_${i}.lp ahl settingsfile.set >> logs/ca.txt
  ./build/reformulate benchmarks/CA/instance_${i}.lp ahl_poor settingsfile.set >> logs/ca.txt
  ./build/reformulate benchmarks/CA/instance_${i}.lp kp settingsfile.set >> logs/ca.txt
done



