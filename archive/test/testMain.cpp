#include <iostream>
#include <chrono>

#include "test.hpp"


int main(int argc, char const *argv[])
{
    // test tabix
    if (argc < 1)
        return 0;

    std::string fileName(argv[1]);
    auto tic = std::chrono::high_resolution_clock::now();
    if (test_tabix(fileName) != 0)
    {
        std::cerr << "TEST_TABIX: failed.";
        return 1;
    }
    auto toc = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = toc - tic;
    std::cout << "TEST_TABIX: passed. " << elapsed.count() << " s.\n";
    
    // test argument parser

    return 0;
}
