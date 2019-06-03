#include <iostream>
#include <string>
#include <htslib/bgzf.h>

#include "version.hpp"


int extractumi(int argc, char const *argv[])
{
    int counter = 0;
    kstring_t hdr = {0, 0, NULL};
    kstring_t seq = {0, 0, NULL};
    kstring_t opt = {0, 0, NULL};
    kstring_t qual = {0, 0, NULL};
    

    BGZF *bgzf_fp = bgzf_open(argv[1], "r");
    if (bgzf_fp == nullptr) {
        std::cerr << "extractumi :: error, could not open file: " << argv[1] << std::endl;
        return 1;
    }

    auto tic = std::chrono::high_resolution_clock::now();

    while (bgzf_getline(bgzf_fp, '\n', &hdr) > 0) {
        bgzf_getline(bgzf_fp, '\n', &seq);
        bgzf_getline(bgzf_fp, '\n', &opt);
        bgzf_getline(bgzf_fp, '\n', &qual);

        std::string header(hdr.s);
        std::string sequence(seq.s);
        std::string quality(qual.s);

        int pos = header.find_first_of(' ');

        std::cout << header << std::endl;
        std::cout << header.substr(0, pos) << std::endl;
        std::cout << sequence << std::endl;
        
        break;

        //std::cout << hdr.s << std::endl;
        //std::cout << seq.s << std::endl;
        //std::cout << opt.s << std::endl;
        //std::cout << qual.s << std::endl;


        counter++;

        //if (counter == 1)
        //    break;
    }
    
    auto toc = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = toc - tic;
    std::cout << "ribotools extractumi " << elapsed.count() << " s.\n";
    std::cout << "counter: " << counter << std::endl;

    free(hdr.s);
    free(seq.s);
    free(opt.s);
    free(qual.s);
    bgzf_close(bgzf_fp);
    
    return 0;
}