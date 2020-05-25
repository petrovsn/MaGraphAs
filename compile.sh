#/usr/bash

cd src
cd libs

g++ -O2 -c -std=c++11 Node.cpp -I../libs


cd ..
g++  -E main.cpp -o GGT_10KG1.cpp -I libs/ -std=c++11
g++   GGT_10KG1.cpp libs/Node.o -o ../MaGraphAs -pthread -lm -L libs -I libs -std=c++11
cd ..
#./GGT_10KG1

