#include <iostream>

#define VERSION 1.01

void version();
void usage();
int main_codonfc(int argc, const char *argv[]);


int main(int argc, const char *argv[])
{
    if (argc < 2) {usage(); return 1;}

    if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {usage(); return 1;}
    
    else if (strcmp(argv[1], "-v") == 0 || strcmp(argv[1], "--version") == 0) {version(); return 1;}

    else if (strcmp(argv[1], "codonfc") == 0) {main_codonfc(argc++, argv++);}

    else {usage(); return 1;}
    
    return 0;
}


void version()
{
    std::cout << "Ribotools " << VERSION << std::endl;
}


void usage()
{
    std::cout << "Ribotools help" << std::endl;
}