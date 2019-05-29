#include "parsercommands.hpp"

ParserCommands::ParserCommands(const int &argc, const char **argv)
{
    for (int i = 1; i < argc; ++i)
        m_tokens.push_back(std::string(argv[i]));
    
    m_iter = m_tokens.end();
}

bool ParserCommands::findOption(const std::string &option)
{
    m_iter = std::find(m_tokens.begin(), m_tokens.end(), option);
    return m_iter != m_tokens.end();
}

bool ParserCommands::nextArgument(std::string &argument)
{   
    ++m_iter;
    if (m_iter == m_tokens.end()) {
        argument = "";
        return false;
    }
        
    
    argument = *m_iter;
    if (argument.rfind("-", 0) == 0) {
        argument = "";
        return false;
    }
        
    return true;
}
