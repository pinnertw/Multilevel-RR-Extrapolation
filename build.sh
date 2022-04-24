g++ -O3 -std=c++11 -Isrc -c src/structural_params.cpp -o obj/structural_params.o
g++ -O3 -std=c++11 -Isrc -c src/models.cpp -o obj/model.o
g++ -O3 -std=c++11 -Isrc -c src/estimator.cpp -o obj/estimator.o

g++ -O3 -std=c++11 -Isrc obj/* main.cpp -o test.out
g++ -O3 -std=c++11 -Isrc obj/* BS_call.cpp -o call.out
g++ -O3 -std=c++11 -Isrc obj/* BS_lookback.cpp -o lookback.out
g++ -O3 -std=c++11 -Isrc obj/* BS_barrier.cpp -o barrier.out
g++ -O3 -std=c++11 -Isrc obj/* compound.cpp -o compound.out


# To comment a line, use #
./test.out
#./call.out
#./lookback.out
#./barrier.out
#./compound.out
