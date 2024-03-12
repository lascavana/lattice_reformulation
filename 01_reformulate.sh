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
  ./build/reformulate benchmarks/struct_s/instance_${i}.lp ahl settingsfile.set >> logs/struct_s_ahl.txt
  ./build/reformulate benchmarks/struct_s/instance_${i}.lp ahl_diag settingsfile.set >> logs/struct_s_ahl_diag.txt
  ./build/reformulate benchmarks/struct_s/instance_${i}.lp ahl_poor settingsfile.set >> logs/struct_s_ahl_poor.txt
  ./build/reformulate benchmarks/struct_s/instance_${i}.lp kp settingsfile.set >> logs/struct_s_kp.txt
done

# struct_b
for i in {1..30}
do
  ./build/reformulate benchmarks/struct_b/instance_${i}.lp ahl settingsfile.set >> logs/struct_b_ahl.txt
  ./build/reformulate benchmarks/struct_b/instance_${i}.lp ahl_diag settingsfile.set >> logs/struct_b_ahl_diag.txt
  ./build/reformulate benchmarks/struct_b/instance_${i}.lp ahl_poor settingsfile.set >> logs/struct_b_ahl_poor.txt
  ./build/reformulate benchmarks/struct_b/instance_${i}.lp kp settingsfile.set >> logs/struct_b_kp.txt
done

# nostruct_s
for i in {1..30}
do
  ./build/reformulate benchmarks/nostruct_s/instance_${i}.lp ahl settingsfile.set >> logs/nostruct_s_ahl.txt
  ./build/reformulate benchmarks/nostruct_s/instance_${i}.lp ahl_diag settingsfile.set >> logs/nostruct_s_ahl_diag.txt
  ./build/reformulate benchmarks/nostruct_s/instance_${i}.lp ahl_poor settingsfile.set >> logs/nostruct_s_ahl_poor.txt
  ./build/reformulate benchmarks/nostruct_s/instance_${i}.lp kp settingsfile.set >> logs/nostruct_s_kp.txt
done

# nostruct_b
for i in {1..30}
do
  ./build/reformulate benchmarks/nostruct_s/instance_${i}.lp ahl settingsfile.set >> logs/nostruct_b_ahl.txt
  ./build/reformulate benchmarks/nostruct_s/instance_${i}.lp ahl_diag settingsfile.set >> logs/nostruct_b_ahl_diag.txt
  ./build/reformulate benchmarks/nostruct_s/instance_${i}.lp ahl_poor settingsfile.set >> logs/nostruct_b_ahl_poor.txt
  ./build/reformulate benchmarks/nostruct_s/instance_${i}.lp kp settingsfile.set >> logs/nostruct_b_kp.txt
done

# Market split
for i in {1..30}
do
  ./build/reformulate benchmarks/MS/instance_${i}.lp ahl settingsfile.set >> logs/ms_ahl.txt
  ./build/reformulate benchmarks/MS/instance_${i}.lp ahl_poor settingsfile.set >> logs/ms_ahl_poor.txt
  ./build/reformulate benchmarks/MS/instance_${i}.lp kp settingsfile.set >> logs/ms_kp.txt
done

# Generalized assignment 
for i in {1..30}
do
  ./build/reformulate benchmarks/GAP/instance_${i}.lp ahl settingsfile.set >> logs/gap_ahl.txt
  ./build/reformulate benchmarks/GAP/instance_${i}.lp ahl_poor settingsfile.set >> logs/gap_ahl_poor.txt
  ./build/reformulate benchmarks/GAP/instance_${i}.lp kp settingsfile.set >> logs/gap_kp.txt
done

# Combinatorial auctions
for i in {1..30}
do
  ./build/reformulate benchmarks/CA/instance_${i}.lp ahl settingsfile.set >> logs/ca_ahl.txt
  ./build/reformulate benchmarks/CA/instance_${i}.lp ahl_poor settingsfile.set >> logs/ca_ahl_poor.txt
  ./build/reformulate benchmarks/CA/instance_${i}.lp kp settingsfile.set >> logs/ca_kp.txt
done

