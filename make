#! /bin/bash

g++ -o Plot -Ofast -m64 -march=native -mtune=native -std=gnu++20 main.cpp BMP/EasyBMP.cpp && echo Made! && echo
