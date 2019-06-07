#ifndef PARSERARGV_H
#define PARSERARGV_H

#include <string>
#include <vector>
#include <algorithm>

class ParserArgv
{
public:
    explicit ParserArgv(const int argc, const char *argv[]);
    bool find(const std::string &argument);
    bool next(std::string &argument);

private:
    std::vector<std::string> m_tokens;
    std::vector<std::string>::const_iterator m_iter;
};

#endif