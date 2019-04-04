#include "coverage.h"

Coverage::Coverage(char *fileName) :
    m_buffer({0, 0, 0}),
    m_iterator(nullptr),
    m_handleFile(nullptr),
    m_handleIndex(nullptr)
{
    m_handleFile = hts_open(fileName, "r");
    if (!m_handleFile) {
        std::cerr << "Error::Coverage::Constructor: could not open coverage file " << fileName << std::endl;
        return;
    }

    m_handleIndex = tbx_index_load(fileName);
    if (!m_handleIndex) {
        std::cerr << "Error::Coverage::Constructor: could not load coverage index " << fileName << std::endl;
        return;
    }

    std::cout << "Coverage::Constructor: file and index handles are created " << fileName << std::endl;
}

Coverage::~Coverage()
{
    std::cout << "Coverage::Destructor" << std::endl;

    if (m_iterator != nullptr)
        tbx_itr_destroy(m_iterator);
    
    if (m_handleFile != nullptr)
        hts_close(m_handleFile);

    if (m_handleIndex != nullptr)
        tbx_destroy(m_handleIndex);
    
    free(ks_release(&m_buffer));
}


void Coverage::query(char *chrom, int32_t chromStart, int32_t chromEnd)
{
    int32_t tid = tbx_name2id(m_handleIndex, chrom);
    if (tid == -1) {
        std::cerr << "Coverage::Query: the sequence not present in this file " << chrom << std::endl;
        return;
    }
    
    m_iterator = tbx_itr_queryi(m_handleIndex, tid, chromStart, chromEnd + 1);
    if (!m_iterator) {
        std::cerr << "Coverage::Query: failed to query index " << chrom << ":" << chromStart << "-" << chromEnd << std::endl;
        return;
    }
    
    int32_t counter = 0;
    while(tbx_itr_next(m_handleFile, m_handleIndex, m_iterator, &m_buffer) >= 0) {
        counter++;
    }
    
    std::cout << "Coverage::Query: " << chrom << ":" << chromStart << "-" << chromEnd << " " << counter << std::endl;    
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	
    // Get the command string
    char cmd[64];
	if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
		mexErrMsgTxt("First input should be a command string less than 64 characters long.");
        
    // New
    if (!strcmp("new", cmd)) {
        // Check parameters
        if (nlhs != 1)
            mexErrMsgTxt("New: One output expected.");


        char fileName[128];
        mxGetString(prhs[1], fileName, sizeof(fileName));

        // Return a handle to a new C++ instance
        plhs[0] = convertPtr2Mat<Coverage>(new Coverage(fileName));
        return;
    }
    
    // Check there is a second input, which should be the class instance handle
    if (nrhs < 2)
		mexErrMsgTxt("Second input should be a class instance handle.");
    
    // Delete
    if (!strcmp("delete", cmd)) {
        // Destroy the C++ object
        destroyObject<Coverage>(prhs[1]);
        // Warn if other commands were ignored
        if (nlhs != 0 || nrhs != 2)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        return;
    }
    
    // Get the class instance pointer from the second input
    Coverage* obj = convertMat2Ptr<Coverage>(prhs[1]);
    
    // Call the various class methods
    // query  
    if (!strcmp("query", cmd)) {
        // Check parameters
        //if (nlhs < 1 || nrhs < 6)
        //    mexErrMsgTxt("Query: Unexpected arguments.");

        // Get query parameters
        char chrom[128];
        mxGetString(prhs[2], chrom, sizeof(chrom));
        int32_t chromStart = (int32_t)*mxGetPr(prhs[3]);
        int32_t chromEnd = (int32_t)*mxGetPr(prhs[4]);
        //int32_t *track;
        //plhs[0] = mxCreateNumericMatrix(1, trackSize, mxINT32_CLASS, mxREAL);
        //track = (int32_t*)mxGetPr(plhs[0]);
        //track[0] = 1;
        //track[1] = 2;

        // Call the method
        //obj->query(track, chromStart, trackSize);
        obj->query(chrom, chromStart, chromEnd);
        
        // return result
        mexPrintf("chrom = %s\n", chrom);
        mexPrintf("chromStart = %d\n", chromStart);
        mexPrintf("chromEnd = %d\n", chromEnd);
        
        return;
    }
    
    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}