#ifndef ARGPARSER_H
#define ARGPARSER_H

#include <iostream>

class ArgParser
{
public:
    ArgParser();
    ~ArgParser();

    void parse(int argc, char const *argv[]);

private:
    //std::vector<std::string> *m_keys;
    //std::vector<std::string> *m_values;
    //std::vector<int> *m_index;


};

#endif /* ARGPARSER_H */