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
    std::cout << "Coverage::destructor" << std::endl;

    if (m_iterator != nullptr)
        tbx_itr_destroy(m_iterator);
    
    if (m_handleFile != nullptr)
        hts_close(m_handleFile);

    if (m_handleIndex != nullptr)
        tbx_destroy(m_handleIndex);
    
    free(ks_release(&m_buffer));
}


void Coverage::query(int32_t *track, int32_t chromStart, int32_t trackSize)
{
    std::cout << "Coverage::query:: " << chromStart << std::endl;

    for (int k = 0; k < trackSize; ++k) {
        track[k] = k;
    }
        
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
        if (nlhs < 1 || nrhs < 3)
            mexErrMsgTxt("Query: Unexpected arguments.");

        // Get query parameters
        int32_t chromStart = (int32_t)*mxGetPr(prhs[2]);
        int32_t trackSize = (int32_t)*mxGetPr(prhs[3]);
        int32_t *track;
        plhs[0] = mxCreateNumericMatrix(1, trackSize, mxINT32_CLASS, mxREAL);
        track = (int32_t*)mxGetPr(plhs[0]);
        track[0] = 1;
        track[1] = 2;

        // Call the method
        obj->query(track, chromStart, trackSize);

        // return result
        mexPrintf("Value01 = %d\n", track[0]);
        return;
    }
    
    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}