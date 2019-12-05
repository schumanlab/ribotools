#include <iostream>
#include <map>
#include <algorithm>
#include <numeric>

#include "argumentparser.h"
#include "bamio.h"
#include "bedio.h"
#include "version.h"


struct PSiteOffset {
    int32_t readLength, offsetPrev, offsetNext = 0;
};

bool operator<(PSiteOffset const& lhs, PSiteOffset const& rhs) {
        if (lhs.readLength != rhs.readLength)
            return lhs.readLength < rhs.readLength;
        else if (lhs.offsetPrev != rhs.offsetPrev)
            return lhs.offsetPrev < rhs.offsetPrev;
        else
            return lhs.offsetNext < rhs.offsetNext;
}



struct PSiteOffsetComparator {

    // opt into being transparent comparator
    using is_transparent = void;

    bool operator() (PSiteOffset const& lhs, PSiteOffset const& rhs) const {
        return lhs < rhs;
    }

    bool operator() (int lhs, PSiteOffset const& rhs) const {
        return lhs <= rhs.readLength;
    }

    bool operator() (PSiteOffset const& lhs, int rhs) const {
        return lhs.readLength <= rhs;
    }
};


int closestNumberInFrame(int n, int m = 3) {

    // quotient
    int q = n / m;

    // left number
    int n_left = m * q;

    // right number
    int n_right = (n * m) > 0 ? (m * (q + 1)) : (m * (q - 1));

    // choose closest
    if (std::abs(n - n_left) < std::abs(n - n_right))
        return n_left;

    return n_right;
}



void accumulatePSiteOffset(std::map<PSiteOffset, uint32_t, PSiteOffsetComparator> &map_pSiteOffset, BedIO &hBed, BamAuxiliary &hBam) {

    int readCount = 0;
    int bedCount = 0;

    // loop over each bed record
    while (hBed.next()) {

        hBam.query(hBed.bed().name(1), hBed.bed().orfStart(), hBed.bed().orfEnd());
        int32_t thickStart = hBed.bed().orfStart();

        while (hBam.next()) {

            int32_t readLength = hBam.readLength();
            int32_t readStart = hBam.readStart();
            int32_t readEnd = hBam.readEnd();

            int32_t shiftStart = readStart - thickStart;
            int32_t shiftEnd = readEnd - thickStart;

            // correct for reads not spannig start codon
            if (thickStart < readStart) {
                int steps = (readEnd - thickStart) / readLength;
                int anchorStart = closestNumberInFrame(steps * readLength);
                shiftStart -= anchorStart;
                shiftEnd -= anchorStart;
            }

            PSiteOffset key = {readLength, shiftStart, shiftEnd};
            map_pSiteOffset[key]++;
            readCount++;
        }

        bedCount++;
        //if(bedCount == 1000) break;
    }

    std::cerr << "used reads: " << readCount << std::endl;
}


void bestPSiteOffset(std::map<int, int> &map_pSiteBest, std::map<PSiteOffset, uint32_t, PSiteOffsetComparator> &map_pSiteOffset, bool flagPrintOffset) {

    const std::vector<int> offset = {12,13,11};

    if (flagPrintOffset)
        std::cout << "#read.length\tpsite.offset\treads.frame0\treads.frame1\treads.frame2" << std::endl;

    // calculate best frame per length
    for (int l = 1; l <= 100; ++l) {

        std::map<PSiteOffset, uint32_t, PSiteOffsetComparator>::const_iterator itBegin = map_pSiteOffset.upper_bound(l);
        std::map<PSiteOffset, uint32_t, PSiteOffsetComparator>::const_iterator itEnd = map_pSiteOffset.lower_bound(l);

        // skip empty lengths
        if (itBegin == itEnd) continue;

        // accumulate counts per frame
        std::vector<int> counts(3, 0);
        for (std::map<PSiteOffset, uint32_t>::const_iterator it = itBegin; it != itEnd; ++it) {
            std::size_t frame = static_cast<std::size_t>(std::abs(it->first.offsetPrev) % 3);
            if (frame < 3)
                counts[frame] += it->second;
        }

        // best frame per length
        auto itMax = std::max_element(counts.begin(), counts.end());
        std::rotate(counts.begin(), itMax, counts.end());
        int bestOffset = offset.at(static_cast<std::size_t>(itMax - counts.begin()));

        map_pSiteBest[l] = bestOffset;

        if (flagPrintOffset)
            std::cout << l << "\t"
                      << bestOffset << "\t"
                      << counts.at(0) << "\t"
                      << counts.at(1) << "\t"
                      << counts.at(2)
                      << std::endl;
    }
}


void pileUpRegion(BamAuxiliary &hBam,
                  const std::string &chrom,
                  int chromStart,
                  int chromEnd,
                  const std::map<int, int> &map_pSiteBest,
                  std::vector<int> &hist) {

    hBam.query(chrom, chromStart, chromEnd);
    while (hBam.next()) {

        int32_t readLength = hBam.readLength();
        int32_t readStart = hBam.readStart();


        if (map_pSiteBest.find(readLength) == map_pSiteBest.end())
            continue;

        int32_t readPSite = readStart + map_pSiteBest.find(readLength)->second;

        int offset = readPSite - chromStart;
        if (offset < 0) continue;

        std::size_t idx = static_cast<std::size_t>(offset);
        if (idx < hist.size())
            hist.at(idx)++;

    }

}



void pileUpHistogram(const std::map<int, int> &map_pSiteBest,
                     std::vector<int> &histStart,
                     std::vector<int> &histCenter,
                     std::vector<int> &histEnd,
                     BedIO &hBed,
                     BamAuxiliary &hBam) {

    int bedCount = 0;

    // rewind bed file
    hBed.rewind();

    // loop over each bed record
    while (hBed.next()) {

        // start histogram
        pileUpRegion(hBam,
                     hBed.bed().name(1),
                     hBed.bed().orfStart() - 25,
                     hBed.bed().orfStart() + 75,
                     map_pSiteBest,
                     histStart);

        // center histogram
        int orfCenter = hBed.bed().orfStart() + closestNumberInFrame(hBed.bed().orfSpan() / 2);
        pileUpRegion(hBam,
                     hBed.bed().name(1),
                     orfCenter - 50,
                     orfCenter + 50,
                     map_pSiteBest,
                     histCenter);

        // end histogram
        pileUpRegion(hBam,
                     hBed.bed().name(1),
                     hBed.bed().orfEnd() - 75,
                     hBed.bed().orfEnd() + 25,
                     map_pSiteBest,
                     histEnd);

        bedCount++;
        //if(bedCount == 1000) break;
    }

}



int main_periodicity(int argc, const char *argv[])
{
    std::string fileBed;
    std::string fileBam;
    bool flagPrintOffset;

    auto p = ArgumentParser("periodicity", std::string(VERSION), "creates P-site periodicity histogram");
    p.addArgumentRequired("annotation").setKeyShort("-a").setKeyLong("--bed").setHelp("BED file containing transcript annotation");
    p.addArgumentRequired("alignment").setKeyShort("-b").setKeyLong("--bam").setHelp("BAM file for alignment");
    p.addArgumentFlag("offset").setKeyShort("-o").setKeyLong("--offset").setHelp("prints PSite offset per read length");
    p.addArgumentFlag("help").setKeyShort("-h").setKeyLong("--help").setHelp("prints help message");
    p.addArgumentFlag("version").setKeyShort("-v").setKeyLong("--version").setHelp("prints major.minor.build version");


    try {
        p.parse(argc, argv);
        fileBed = p.get<std::string>("annotation");
        fileBam = p.get<std::string>("alignment");
        flagPrintOffset = p.get<bool>("offset");
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }



    auto hBed = BedIO(fileBed);
    if (!hBed.isOpen()) {
        std::cerr << "ribotools::" + hBed.error() << std::endl;
        return EXIT_FAILURE;
    }


    auto hBam = BamAuxiliary(fileBam, 255);
    if (!hBam.isOpen()) {
        std::cerr << "ribotools::" + hBam.what() << std::endl;
        return EXIT_FAILURE;
    }


    std::map<PSiteOffset, uint32_t, PSiteOffsetComparator> map_pSiteOffset;
    std::map<int, int> map_pSiteBest;
    std::vector<int> histStart(100, 0);
    std::vector<int> histCenter(100, 0);
    std::vector<int> histEnd(100, 0);

    accumulatePSiteOffset(map_pSiteOffset, hBed, hBam);
    bestPSiteOffset(map_pSiteBest, map_pSiteOffset, flagPrintOffset);
    pileUpHistogram(map_pSiteBest, histStart, histCenter, histEnd, hBed, hBam);

    std::cout << "#index\thist.start\thist.center\thist.end" << std::endl;
    for (std::size_t i = 0; i < 100; ++i) {
        std::cout << i << "\t" << histStart[i] << "\t" << histCenter[i] << "\t" << histEnd[i] << std::endl;
    }


    return 0;
}



