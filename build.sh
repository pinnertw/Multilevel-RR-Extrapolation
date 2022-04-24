mkdir obj
mkdir executable

echo "Compile object files..."
for i in structural_params models estimator
do
    g++ -O3 -Isrc -std=c++11 -c src/$i.cpp -o obj/$i.o
done
echo "Done!"

echo "Compile main files..."
for i in BS_call BS_lookback BS_barrier compound
do
    g++ -O3 -Isrc -std=c++11 obj/* simulations/$i.cpp -o executable/$i.out
done
echo "Done!"
echo "Installation done!"
