# compile code
rm -r build
mkdir build
cd build
cmake .. -DSCIP_DIR=$SCIPOPTDIR
make
cd ..

rm -r logs
mkdir logs

# MIPLIB
for i in benchmarks/MIPLIB/*.mps; do
    echo $i
    ./build/reformulate $i ahl settingsfile.set
    ./build/reformulate $i kp settingsfile.set
done