#include <iostream>
#include <string>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/tbx.h>

#include <chrono>

#include "version.hpp"
#include "bedrecord.hpp"

inline uint32_t bam_calend(const bam1_core_t *core, const uint32_t *cigar)
{
    return core->pos + (core->n_cigar ? bam_cigar2rlen(core->n_cigar, cigar) : 1);
}


int countasite(int arg, char const *argv[])
{
    auto tic = std::chrono::high_resolution_clock::now();
    std::string fileBed(argv[1]);
    std::string fileBam(argv[2]);

    // open tabix file
    htsFile *fhBed = hts_open(fileBed.c_str(), "r");
    tbx_t *fhTabix = tbx_index_load(fileBed.c_str());
    hts_itr_t *iter = nullptr;

    // loop over bam file
    samFile *fhBam = sam_open(fileBam.c_str(), "r");
    bam_hdr_t *hdrBam = sam_hdr_read(fhBam);
    bam1_t *aln = bam_init1();
    kstring_t buffer = {0, 0, nullptr};
    
    auto bed = BedRecord();

    int counter = 0;
    int counterReused = 0;
    int counterNew = 0;
    int counterIntergenic = 0;
    while (sam_read1(fhBam, hdrBam, aln) >= 0) {
        counter++;

        std::string readChrom = std::string(hdrBam->target_name[aln->core.tid]);
        uint32_t readStart = aln->core.pos;
        uint32_t readEnd = bam_calend(&aln->core, bam_get_cigar(aln));
        int32_t readLength = aln->core.l_qseq;

        if (bed.overlap(readChrom, readStart, readEnd)) {
            counterReused++;
        }
        else {
            int i = 0;
            iter = tbx_itr_queryi(fhTabix, aln->core.tid, readStart, readStart+1);
            while (tbx_itr_next(fhBed, fhTabix, iter, &buffer) >= 0) {
                auto ss = std::stringstream(buffer.s);
                ss >> bed;
                i++;
            }
            

            if (i == 0) {
                counterIntergenic++;
            }
            else {
                counterNew++;
            }
        }        

        //std::cout << i << std::endl;
        if (counter > 2000000)
            break;
    }


    tbx_itr_destroy(iter);
    tbx_destroy(fhTabix);
    hts_close(fhBed);
    free(buffer.s);
    bam_destroy1(aln);
    bam_hdr_destroy(hdrBam);
    sam_close(fhBam);

    std::cout << "Counts: " << counter << std::endl;
    std::cout << "Reused: " << counterReused << std::endl;
    std::cout << "New: " << counterNew << std::endl;
    std::cout << "Intergenic: " << counterIntergenic << std::endl;

    auto toc = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = toc - tic;
    std::cout << "ribotools countasite " << elapsed.count() << " s.\n";
    return 0;
}
