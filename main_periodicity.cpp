#include <iostream>
#include <map>

#include "argumentparser.h"
#include "bamio.h"
#include "bedio.h"
#include "version.h"


struct PSiteShift {
    int32_t readLength, shiftPrev, shiftNext = 0;
};

bool operator<(PSiteShift const& lhs, PSiteShift const& rhs) {
        if (lhs.readLength != rhs.readLength)
            return lhs.readLength < rhs.readLength;
        else if (lhs.shiftPrev != rhs.shiftPrev)
            return lhs.shiftPrev < rhs.shiftPrev;
        else
            return lhs.shiftNext < rhs.shiftNext;
}



struct PSiteShiftComparator {

    // opt into being transparent comparator
    using is_transparent = void;

    bool operator() (PSiteShift const& lhs, PSiteShift const& rhs) const {
        return lhs < rhs;
    }

    bool operator() (int lhs, PSiteShift const& rhs) const {
        return lhs < rhs.readLength;
    }

    bool operator() (PSiteShift const& lhs, int rhs) const {
        return lhs.readLength < rhs;
    }
};



/*
struct PSiteShift {
    int32_t readLength, shiftFPrime, shiftTPrime = 0;

    bool operator==(const PSiteShift &qry) const {
        return (readLength == qry.readLength) && (shiftFPrime == qry.shiftFPrime) && (shiftTPrime == qry.shiftTPrime);
    }

    bool operator<(const PSiteShift &qry) const {
        if (readLength != qry.readLength)
            return readLength < qry.readLength;
        else if (shiftFPrime != qry.shiftFPrime)
            return shiftFPrime < qry.shiftFPrime;
        else
            return shiftTPrime < qry.shiftTPrime;
    }
};
*/


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



void calculatePSiteShift(std::map<PSiteShift, uint32_t, PSiteShiftComparator> &map_pSiteShift, BedIO &hBed, BamAuxiliary &hBam) {

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

            PSiteShift key = {readLength, shiftStart, shiftEnd};
            map_pSiteShift[key]++;
            readCount++;
        }

        bedCount++;
        if(bedCount == 500) break;
    }

    std::cerr << "used reads: " << readCount << std::endl;
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


    std::map<PSiteShift, uint32_t, PSiteShiftComparator> map_pSiteShift;
    calculatePSiteShift(map_pSiteShift, hBed, hBam);

    // print options
    if (flagPrintOffset) {
        //for (const auto &[key, value] : map_pSiteShift) {
        //    std::cout << key.readLength << "\t" << key.shiftPrev << "\t" << key.shiftNext << "\t" << value << std::endl;
        //}

        std::map<PSiteShift, uint32_t, PSiteShiftComparator>::const_iterator itBegin = map_pSiteShift.upper_bound(28);
        std::map<PSiteShift, uint32_t, PSiteShiftComparator>::const_iterator itEnd = map_pSiteShift.lower_bound(32);


        for(std::map<PSiteShift, uint32_t, PSiteShiftComparator>::const_iterator it = itBegin;
            it != itEnd; ++it) {
            std::cout << it->first.readLength << "\t"
                      << it->first.shiftPrev << "\t"
                      << it->first.shiftNext << "\t"
                      << it->second << std::endl;
        }





    }


    return 0;
}



