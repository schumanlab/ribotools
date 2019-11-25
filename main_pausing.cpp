#include <iostream>
#include <memory>
#include <numeric>

//#include <fstream>

//#include <cmath>
//#include <htslib/faidx.h>

#include "argumentparser.h"
//#include "bamhandle.h"
//#include "bedrecord.h"
#include "bedio.h"
#include "bamio.h"
#include "seqio.h"

#include "aminoacidtable.h"
#include "version.h"

bool cmpf(float A, float B, float epsilon = 0.000005f)
{
    return (std::fabs(A - B) < epsilon);
}


std::vector<float> reduceToCodons(const std::vector<int> &nucleotides) {
    std::vector<float> codons(nucleotides.size()/3, 0);
    const float codon_step = 3.0;

    for (std::size_t i = 0; i < codons.size(); ++i) {
        codons.at(i) = (nucleotides.at(3*i) + nucleotides.at(3*i+1) + nucleotides.at(3*i+2)) / codon_step;
    }
    return codons;
}



void calculatePauseScore(std::vector<float>::iterator ov,
                                       std::vector<float>::iterator ovEnd,
                                       const std::vector<float> &coverage,
                                       int basic, int flank) {

    // calculate basic score
    float backgroundBasicSum = 0.0;
    int backgroundBasicCount = std::min(basic, static_cast<int>(coverage.size()));


    for (std::vector<float>::const_iterator iv = coverage.begin(); iv != coverage.begin() + backgroundBasicCount; ++iv)
        backgroundBasicSum += *iv;
    float backgroundBasic = backgroundBasicSum / backgroundBasicCount;

    // calculate sliding score
    std::vector<float>::const_iterator ivPrev = coverage.begin();
    std::vector<float>::const_iterator ivNext = coverage.begin() + flank;

    // initialize previous background
    float backgroundPrevSum = 0.0;
    int backgroundPrevCount = 0;

    // initialize next background
    float backgroundNextSum = 0.0;
    int backgroundNextCount = flank;
    for (std::vector<float>::const_iterator iv = coverage.begin();
         iv != coverage.begin() + flank; ++iv)
        backgroundNextSum += *iv;

    // loop coverage
    for (std::vector<float>::const_iterator iv = coverage.begin();
         iv != coverage.end(); ++iv) {

        // background next
        backgroundNextSum -= *iv;
        if (ivNext < coverage.end()) {
            backgroundNextSum += *ivNext;
            ++ivNext;
        }
        else {
            backgroundNextCount--;
        }
        float backgroundNext = (backgroundNextCount > 0) ? (backgroundNextSum / backgroundNextCount) : backgroundNextSum;

        // background previous
        float backgroundPrev = (backgroundPrevCount > 0) ? (backgroundPrevSum / backgroundPrevCount) : backgroundPrevSum;
        backgroundPrevSum += *iv;
        backgroundPrevCount++;
        if ((iv - ivPrev) == flank) {
            backgroundPrevSum -= *ivPrev;
            ++ivPrev;
            backgroundPrevCount--;
        }

        // calculate z-score
        float background = std::max(backgroundBasic, std::max(backgroundPrev, backgroundNext));
        if (ov != ovEnd)
        *ov = (background > 0) ? (*iv - background) / std::sqrt(background) : 0.0;
        ++ov;
    }
}


void calculatePauseScore(std::vector<float>::iterator ov,
                         std::vector<float>::iterator ovEnd,
                         const std::vector<float> &coverage,
                         int window = 11) {

    // window is  odd -> centered around current position
    // window is even -> centered around current and previous position
    int windowPrev = window / 2;
    int windowNext = windowPrev;
    if (window % 2 == 0)
        windowNext -= 1;

    std::vector<float>::const_iterator ivPrev = coverage.begin();
    std::vector<float>::const_iterator ivNext = coverage.begin() + windowNext;

    // intialize sum
    float windowSum = 0.0f;
    for (std::vector<float>::const_iterator iv = ivPrev; iv != (ivNext + 1); ++iv)
        windowSum += *iv;

    for (std::vector<float>::const_iterator iv = coverage.begin();
         iv != coverage.end(); ++iv) {

        // assign output
        if (ov != ovEnd) {
            float background = windowSum / (ivNext - ivPrev + 1); // sliding window average
            *ov = (background > 0.0f) ? (*iv - background) / std::sqrt(background) : 0.0f;
        }
        ++ov;

        // update next step
        if (ivNext != coverage.end() - 1) {
            ++ivNext;
            windowSum += *ivNext;
        }

        // update previous step
        if ((iv - ivPrev) == windowPrev) {
            windowSum -= *ivPrev;
            ++ivPrev;
        }


    }


}




int main_pausing(const int argc, const char *argv[])
{
    std::string fileBed;
    std::string fileFasta;
    std::vector<std::string> filesBam;
    //std::vector<std::shared_ptr<BamHandle>> handlesBam;
    int backgroundWindow_basic;
    int backgroundWindow_flank;
    int skipCodons;

    auto p = ArgumentParser("pausing", std::string(VERSION), "calculates z-score pausing score per codon");
    p.addArgumentRequired("annotation").setKeyShort("-a").setKeyLong("--bed").setHelp("BED file containing transcript annotation");
    //p.addArgumentRequired("sequence").setKeyShort("-f").setKeyLong("--fasta").setHelp("FASTA file containing transcript sequence");
    p.addArgumentOptional("basic").setKeyShort("-b").setKeyLong("--basic").setDefaultValue<int>(150).setHelp("number of codons to calculate background score");
    p.addArgumentOptional("flank").setKeyShort("-f").setKeyLong("--flank").setDefaultValue<int>(5).setHelp("number of flanking codons to calculate background score");
    p.addArgumentOptional("skip").setKeyShort("-s").setKeyLong("--skip").setDefaultValue<int>(10).setHelp("number of codons to skip after/before start/stop codon");
    p.addArgumentPositional("alignment").setCount(-1).setHelp("list of RiboSeq BAM files");
    p.addArgumentFlag("help").setKeyShort("-h").setKeyLong("--help").setHelp("prints help message");
    p.addArgumentFlag("version").setKeyShort("-v").setKeyLong("--version").setHelp("prints major.minor.build version");

    try {
        p.parse(argc, argv);
        fileBed = p.get<std::string>("annotation");
        //fileFasta = p.get<std::string>("sequence");
        backgroundWindow_basic = p.get<int>("basic");
        backgroundWindow_flank = p.get<int>("flank");
        skipCodons = p.get<int>("skip");
        filesBam = p.get<std::vector<std::string>>("alignment");
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }


    auto hBed = BedIO(fileBed);
    if (!hBed.isOpen()) {
        std::cerr << "ribotools::" + hBed.error() << std::endl;
        return EXIT_FAILURE;
    }

    /*
    auto hFas = SeqIO(fileFasta);
    if (!hFas.isOpen()) {
        std::cerr << "ribotools::" + hFas.error() << std::endl;
        return EXIT_FAILURE;
    }
    */

    std::vector<std::shared_ptr<BamIO>> hBams;
    for (auto fileName : filesBam) {
        auto handle = std::make_shared<BamIO>(fileName);
        if (!handle->isOpen()) {
            std::cerr << "ribotools::" + handle->error() << std::endl;
            return EXIT_FAILURE;
        }
        hBams.emplace_back(handle);
    }


    // print header
    std::cout << "#gene\tcodon\tregion";
    for (auto handle : hBams)
        std::cout << "\t" << handle->name();
    std::cout << std::endl;

    while (hBed.next()) {


        //std::cout << hBed.bed().name(1) << ":" << hBed.bed().orfStart() << "-" << hBed.bed().orfEnd() << std::endl;


        // pause score per BAM
        int handleCounter = 0;
        int spanInCodons = hBed.bed().orfSpan() / 3;
        std::vector<float> pauseScore(static_cast<std::size_t>(spanInCodons) * hBams.size(), 0.0f); // buffer container

        bool useFlag = false;

        for (auto handle : hBams) {

            // count reads
            handle->query(hBed.bed().name(1), hBed.bed().orfStart(), hBed.bed().orfEnd());
            int readCount = handle->count();
            float readsPerCodon = static_cast<float>(readCount) / spanInCodons;
            if (readsPerCodon > 0.1f) {

                // activate usage
                useFlag = true;

                // depth
                std::vector<int> orfDepth = handle->depth(hBed.bed().name(1), hBed.bed().orfStart(), hBed.bed().orfEnd());

                // reduce to codons
                std::vector<float> orfCodons = reduceToCodons(orfDepth);

                // calculate pause score
                std::vector<float>::iterator itBegin = pauseScore.begin() + spanInCodons * handleCounter;
                std::vector<float>::iterator itEnd = pauseScore.begin() + spanInCodons * (handleCounter + 1);
                //calculatePauseScore(itBegin, itEnd, orfCodons);
                for (std::vector<float>::const_iterator iv = orfCodons.begin(); iv != orfCodons.end(); ++iv) {
                    if (itBegin != itEnd)
                        *itBegin = *iv;
                    ++itBegin;
                }


            }

            handleCounter++;
        }


        //for(std::vector<float>::const_iterator iv = pauseScore.begin();
        //    iv != pauseScore.end(); ++iv)
        //    std::cout << *iv << std::endl;

        if (!useFlag) continue;

        // output matrix
        std::vector<float>::iterator it = pauseScore.begin();
        int regionTag = 0;
        for (int i = 0; i < spanInCodons; ++i) {

            std::stringstream outBuffer;

            // update region tag : 0 - initiation, 1 - elongation, 2 - termination
            if (i <= skipCodons)
                regionTag = 0;
            else if ((skipCodons < i) && (i < (spanInCodons - skipCodons)))
                regionTag = 1;
            else
                regionTag = 2;

            outBuffer << hBed.bed().name(2) << "\t" << i << "\t" << regionTag;
            int handleCounter = 0;
            int isZero = 0;
            int isNegative = 0;
            for (auto handle : hBams) {
                float value = *(it + handleCounter * spanInCodons);
                if (cmpf(value, 0.0f)) isZero++;
                if (value < 0.0f) isNegative++;

                outBuffer << "\t" << value;
                handleCounter++;
            }
            outBuffer << std::endl;
            ++it;

            //if ((isNegative < handleCounter/2) && (isZero < handleCounter/2))
                std::cout << outBuffer.str();
        }



        // retrieve ORF sequence
        //hFas.fetch(hBed.bed().name(), hBed.bed().orfStart(), hBed.bed().orfEnd());
        //std::cout << hFas.sequence() << std::endl;


    }





    return 0;
}



/*
    // open BAM handles
    for (const auto &fileName : filesBam) {
        std::shared_ptr<BamHandle> handle = std::make_shared<BamHandle>(fileName, 255, 0);
        handlesBam.emplace_back(handle);
    }

    // open FASTA file
    faidx_t *fhFai = fai_load(fileFasta.c_str());
    if (!fhFai) {
        std::cerr << "ribotools::pausing::error, failed to load FASTA file " << fileFasta << std::endl;
        return 1;
    }

    // open BED file
    std::ifstream fhBed;
    fhBed.open(fileBed);
    if (!fhBed.is_open()) {
        std::cerr << "ribotools::pausing::error, failed to open BED file " << fileBed << std::endl;
        return 1;
    }

    // print header
    //std::cout << "#name\tcodon\tcount\tmin.pausing.score\tmax.pausing.score" << std::endl;

    // loop over BED record
    std::string line;
    //int line_counter = 0;
    while (std::getline(fhBed, line)) {

        // read bed record
        auto bed = BedRecord();
        std::stringstream iss(line);
        iss >> bed;

        int lengthCodons = bed.cdsSpan / 3;

        if (lengthCodons <= 2*skipCodons + 2*backgroundWindow_flank + 1) continue;

        // accumulate coverage from all files
        std::vector<int> codons(static_cast<std::size_t>(lengthCodons), 0);
        for (auto handle : handlesBam)
            handle->calculateSiteCoverage(codons, bed.transcript, 0, bed.span, bed.cdsStart, true);

        // fiter based on coverage and empty codons
        double background_average = static_cast<double>(std::accumulate(codons.begin(), codons.end(), 0)) / lengthCodons;
        if (background_average < 0.1) continue;

        // calculate basic background
        int backgroundWindow_offset = std::min(backgroundWindow_basic, lengthCodons - 2*skipCodons);
        double background_basic = static_cast<double>(std::accumulate(codons.begin() + skipCodons,
                                                                      codons.begin() + skipCodons + backgroundWindow_offset,
                                                                      0)) / backgroundWindow_offset;


        // retrieve sequence
        char *sequence = faidx_fetch_seq(fhFai, bed.name.c_str(), 0, bed.span, &bed.span);
        auto aa = AminoAcidTable();
        int count_paused = 0;

        // count paused codons
        std::vector<int>::const_iterator it;
        std::vector<int>::const_iterator it_prev;
        std::vector<int>::const_iterator it_next;
        for (it = codons.begin() + skipCodons; it != (codons.end() - skipCodons); ++it) {

            // calculate previous background
            double background_prev = 0.0;
            it_prev = it - backgroundWindow_flank;
            if (codons.begin() <= it_prev)
                background_prev = static_cast<double>(std::accumulate(it_prev, it, background_prev)) / backgroundWindow_flank;

            // calculate next background
            double background_next = 0.0;
            it_next = it + backgroundWindow_flank + 1;
            if ((it + 1) < codons.end() && (it_next <= codons.end()))
                background_next = static_cast<double>(std::accumulate(it + 1, it_next, background_next)) / backgroundWindow_flank;

            // best background
            double background = std::max(background_prev, background_prev);
            background = std::max(background, background_basic);

            // calculate zscore
            double zscore = 0.0;
            if (background > 0.1)
                zscore = (*it - background) / std::sqrt(background);

            if (zscore >= 10.0)
                count_paused++;
            // current codon code
            int idx_codon = static_cast<int>(it - codons.begin());
            int idx_nucleotide = idx_codon * 3 + bed.cdsStart;
            char codon[4];
            std::strncpy(codon, &sequence[idx_nucleotide], 3);
            codon[3] = '\0';

            if (zscore != 0.0)
                aa.addPausingScore(std::string(codon), zscore);



            //if (zscore != 0.0)
                //std::cout << bed.gene << "\t" << (it - codons.begin()) << "\t" << codon << "\t" << background << "\t" << zscore << std::endl;
        }

        // output of AminoAcid map
        aa.log(bed.name);
        //std::cout << bed.gene << "\t" << bed.cdsSpan / 3 << "\t" << count_paused << std::endl;

        if (sequence)
            free(sequence);


    }


    // destructors
    fhBed.close();
    fai_destroy(fhFai);
 */



/*
 *  VERSION 1.0
 *
 *
 */
    /*

    if (parser.find("--bam")) {
        std::string fileNameNext;
        while (parser.next(fileNameNext)) {
            BamHandle *handle = new BamHandle(fileNameNext, 255, 0);
            handlesBam.push_back(handle);
        }
    }
    else {
        std::cerr << "ribotools::pausing::error, provide BAM file." << std::endl;
        return 1;
    }


    // open BED file
    std::ifstream fhBed;
    fhBed.open(fileBed);
    if (!fhBed.is_open()) {
        std::cerr << "ribotools::pausing::error, failed to open BED file " << fileBed << std::endl;
        return 1;
    }

    // open FASTA file
    faidx_t *fhFai = fai_load(fileFasta.c_str());
    if (!fhFai) {
        std::cerr << "ribotools::pausing::error, failed to load fasta reference " << fileFasta << std::endl;
        return 1;
    }

    // output header
    std::cout << "# name\tcodon\tn.codons\tn.asites\tsum.zscore" << std::endl;

    // loop over BED record
    std::string line;
    int line_counter = 0;
    while (std::getline(fhBed, line)) {

        // read bed record
        auto bed = BedRecord();
        std::istringstream iss(line);
        iss >> bed;

        int lengthCodons = bed.cdsSpan / 3;

        // check if codon length is too small
        if (lengthCodons <= 2*skipCodons + 2*backgroundWindow_flank + 1) continue;

        // calculate A-site codon coverage
        std::vector<int> codons(static_cast<size_t>(lengthCodons), 0);
        for (auto handle : handlesBam)
            handle->calculateSiteCoverage(codons, bed.transcript, 0, bed.span, bed.cdsStart, true);

        // count empty codons
        //int emptyCodons = static_cast<int>(std::count(codons.begin(), codons.end(), 0));
        //double emptyRatio = static_cast<double>(emptyCodons) / lengthCodons;

        // get average coverage
        double background_average = static_cast<double>(std::accumulate(codons.begin(), codons.end(), 0)) / lengthCodons;

        // fiter based on coverage and empty codons
        //if (background_average < 0.1 || emptyRatio > 0.75) continue;
        if (background_average < 0.1) continue;


        // calculate background
        int backgrounWindow_offset = std::min(backgroundWindow_basic, lengthCodons - 2*skipCodons);
        double background_basic = static_cast<double>(std::accumulate(codons.begin() + skipCodons, codons.begin() + backgrounWindow_offset, 0)) / backgrounWindow_offset;

        // retrieve sequence
        char *sequence = faidx_fetch_seq(fhFai, bed.name.c_str(), 0, bed.span, &bed.span);
        std::vector<int>::const_iterator it;
        std::vector<int>::const_iterator it_prev;
        std::vector<int>::const_iterator it_next;
        auto aa = AminoAcids();

        for (it = codons.begin() + skipCodons; it != (codons.end() - skipCodons); ++it)
        {
            // skip empty codons
            if (*it == 0) continue;

            // calculate previous background
            double background_prev = 0.0;
            it_prev = it - backgroundWindow_flank;
            if (codons.begin() <= it_prev)
                background_prev = static_cast<double>(std::accumulate(it_prev, it, background_prev)) / backgroundWindow_flank;

            // calculate next background
            double background_next = 0.0;
            it_next = it + backgroundWindow_flank + 1;
            if ((it + 1) < codons.end() && (it_next <= codons.end()))
                background_next = static_cast<double>(std::accumulate(it+1, it_next, background_next)) / backgroundWindow_flank;

            // background value
            double background = std::max(background_next, background_prev);
            background = std::max(background, background_basic);

            // skip if background is empty
            if (background == 0.0) continue;

            // zscore
            double zscore = (*it - background) / std::sqrt(background);

            // current codon code
            int idx_codon = static_cast<int>(it - codons.begin());
            int idx_nucleotide = idx_codon * 3 + bed.cdsStart;
            char codon[4];
            std::strncpy(codon, &sequence[idx_nucleotide], 3);
            codon[3] = '\0';

            aa.addTime(std::string(codon), zscore, static_cast<double>(*it));

            // debug output
            //std::cout << idx_codon << "\t" << codon << "\t" << *it << "\t" << background << "\t" << zscore << std::endl;
        }

        // output of AminoAcid map
        aa.log(bed.name);

        if (sequence)
            free(sequence);

        line_counter++;
    }

    std::cerr << "used genes: " << line_counter << std::endl;

    // destructors
    fhBed.close();
    fai_destroy(fhFai);
    for (auto handle : handlesBam)
        delete handle;


    return 0;
    */

