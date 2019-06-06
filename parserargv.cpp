#include "parserargv.hpp"

ParserArgv::ParserArgv(const int argc, const char *argv[])
{
    for (int i = 1; i < argc; ++i)
        m_tokens.push_back(std::string(argv[i]));
    
    m_iter = m_tokens.end();
}


bool ParserArgv::find(const std::string &argument)
{
    m_iter = std::find(m_tokens.begin(), m_tokens.end(), argument);
    return m_iter != m_tokens.end();
}


bool ParserArgv::next(std::string &argument)
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