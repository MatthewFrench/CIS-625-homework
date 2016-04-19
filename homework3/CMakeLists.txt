cmake_minimum_required(VERSION 3.5)
project(homework3)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    contig_generation.h
    kmer_hash.h
    packingDNAseq.h
    pgen.upc
    serial.c
    sort.cpp)

add_executable(homework3 ${SOURCE_FILES})