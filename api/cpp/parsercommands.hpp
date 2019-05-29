#ifndef PARSERCOMMANDS_H
#define PARSERCOMMANDS_H

#include <string>
#include <vector>
#include <algorithm>

class ParserCommands
{
public:
    explicit ParserCommands(const int &argc, const char **argv);

    bool findOption(const std::string &option);
    bool nextArgument(std::string &argument);
    
private:
    std::vector<std::string> m_tokens;
    std::vector<std::string>::const_iterator m_iter;
};

#endif /* PARSERCOMMANDS_H */