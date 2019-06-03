#ifndef TEST_H
#define TEST_H

#include <iostream>
#include <htslib/tbx.h>
#include <htslib/kstring.h>
#include <htslib/bgzf.h>

int test_tabix(const std::string &fileName);

#endif /* TEST_H */