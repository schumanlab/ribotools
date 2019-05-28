#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>

#include <htslib/hts.h>
#include <htslib/faidx.h>
#include <htslib/tbx.h>
#include <htslib/kstring.h>


#include "bedrecord.hpp"


struct AAIndex {
    int index;
    char code;
};

void loadCodonMap(std::unordered_map<std::string, AAIndex> &codonMap, const std::string &fileCodonTable);


int main(int argc, char *argv[])
{
    auto fileCodonMap = std::string(argv[1]);
    auto fileBed = std::string(argv[2]);
    auto fileFasta = std::string(argv[3]);
    auto fileGbed = std::string(argv[4]);
    
    // Load Codon Map
    std::unordered_map<std::string, AAIndex> codonMap;
    loadCodonMap(codonMap, fileCodonMap);

    faidx_t *fai = fai_load(fileFasta.c_str());
    if (!fai) {
        std::cerr << "Error: failed to load fasta reference" << std::endl;
        return -1;
    }

    htsFile *fhGbed = hts_open(fileGbed.c_str(), "r");
    if (!fhGbed) {
        std::cerr << "Error: failed to load Gbed file" << std::endl;
        return -1;
    }

    tbx_t *fhGbedIdx = tbx_index_load(fileGbed.c_str());
    if (!fhGbedIdx) {
        std::cerr << "Error: failed to load Gbed index" << std::endl;
        return -1;
    }

    hts_itr_t *iterGbed = nullptr;

    std::ifstream fh;
    std::string line;
    
    // constructor
    fh.open(fileBed);

    while (std::getline(fh, line)) {
        auto bed = BedRecord();
        std::istringstream isline(line);
        isline >> bed;
        bed.parseExons();
        char *rnaSeq = faidx_fetch_seq(fai, bed.name.c_str(), 0, bed.span, &bed.span);
        
        int queryStart = bed.cdsStart + 60;
        int queryEnd = bed.cdsEnd - 60;
        int querySpan = queryEnd - queryStart;
        double *rnaDepth = (double *)std::calloc(querySpan, sizeof(double));
        int queryTid = tbx_name2id(fhGbedIdx, bed.transcript().c_str());
        kstring_t buffer;
        iterGbed = tbx_itr_queryi(fhGbedIdx, queryTid, queryStart, queryEnd);
        
        double totalDepth = 0.0;
        while (tbx_itr_next(fhGbed, fhGbedIdx, iterGbed, &buffer) >= 0) {
            
            std::string gbedLine = std::string(buffer.s);
            std::istringstream isbuffer(gbedLine);
            std::string chrom;
            int chromStart;
            int chromEnd;
            double chromDepth;
            isbuffer >> chrom;
            isbuffer >> chromStart;
            isbuffer >> chromEnd;
            isbuffer >> chromDepth;

            for (int k = chromStart; k < chromEnd; k++) {
                int index = k - queryStart;
                if ((0 <= index) && (index < querySpan)) {
                    rnaDepth[index] = chromDepth;
                    totalDepth += chromDepth;
                }
            }
        }

        double averageDepth = totalDepth / querySpan;
        
        for (int k = 0; k < querySpan; k+=3) {
            double codonCoverage = rnaDepth[k] + rnaDepth[k+1] + rnaDepth[k+2];
            codonCoverage = (codonCoverage - (static_cast<int>(codonCoverage) % 3)) / 3.0;
            double codonDepth = codonCoverage / averageDepth;
            char codonSeq[4];
            memset(codonSeq, '\0', 4);
            std::strncpy(codonSeq, &rnaSeq[queryStart + k], 3);
            codonSeq[3] = '\0';
            
            if (codonSeq[0] == '\0')
                continue;

            if (codonSeq[1] == '\0')
                continue;

            if (codonSeq[2] == '\0')
                continue;

            if (std::strchr(codonSeq, 'N'))
                continue;

            if ((codonDepth > 0) && (codonDepth <= 10))
                std::cout << codonSeq << "\t" << codonDepth << std::endl;
        }
        
        
        if (rnaSeq)
            free(rnaSeq);
        
        if (rnaDepth)
            free(rnaDepth);
    }

    // destructor
    if (iterGbed != nullptr)
        tbx_itr_destroy(iterGbed);
    
    fh.close();
    fai_destroy(fai);
    hts_close(fhGbed);
    tbx_destroy(fhGbedIdx);
    
    return 0;
}


void loadCodonMap(std::unordered_map<std::string, AAIndex> &codonMap, const std::string &fileCodonTable)
{
    std::ifstream fh;
    int index = 0;
    std::string line;

    fh.open(fileCodonTable);
    while (std::getline(fh, line)) {
        std::istringstream isline(line);
        std::string aaKey;
        char aaValue;
        isline >> aaKey;
        isline >> aaValue;
        codonMap[aaKey] = {index, aaValue};
        index++;
    }
    fh.close();
}
