#ifndef ARGUMENTPARSER_H
#define ARGUMENTPARSER_H

#include<iostream>
#include<vector>
#include <any>
#include <sstream>
#include <utility>
#include <unordered_map>
#include <stdexcept>

//#include <list>
//#include <map>

namespace is_stl_container_helper {
template <typename T> struct is_stl_container:std::false_type{};
template <typename... Args> struct is_stl_container<std::vector<Args...>>:std::true_type{};

// expand stl container helper to other container types
//template <typename... Args> struct is_stl_container<std::list<Args...>>:std::true_type{};
//template <typename... Args> struct is_stl_container<std::map<Args...>>:std::true_type{};
}

template <typename T> struct is_stl_container {
    static constexpr bool const value = is_stl_container_helper::is_stl_container<std::decay_t<T>>::value;
};


/**
 * @brief The Argument class
 *
 * @details describes single argument, defined by user,
 * provides methods to get<value> or get<container>.
 * values are kept as std::string and converted upon request
 *
 */
class Argument
{
public:
    friend class ArgumentParser;
    enum class Rank {POSITIONAL, COMMAND, FLAG, REQUIRED, OPTIONAL};

    explicit Argument(const bool used = false,
                      const Rank rank = Rank::OPTIONAL,
                      const int count = 1,
                      const std::string name = "",
                      const std::string keyShort = "",
                      const std::string keyLong = "",
                      const std::string help = "") :
        m_used(used),
        m_rank(rank),
        m_count(count),
        m_name(name),
        m_keyShort(keyShort),
        m_keyLong(keyLong),
        m_help(help),
        m_valueDefaultText(""),
        m_valueDefault(std::any()),
        m_values(std::vector<std::string>())
    {}

    Argument& setUsed(const bool used) {m_used = std::move(used); return *this;}
    bool used() const {return m_used;}

    Argument& setRank(const Argument::Rank rank) {m_rank = std::move(rank); return *this;}
    Argument::Rank rank() const {return m_rank;}

    Argument& setCount(const int count) {m_count = std::move(count); return *this;}
    int count() const {return m_count;}

    Argument& setName(const std::string name) {m_name = std::move(name); return *this;}
    std::string name() const {return m_name;}

    Argument& setKeyShort(const std::string key) {m_keyShort = std::move(key); return *this;}
    std::string keyShort() const {return m_keyShort;}

    Argument& setKeyLong(const std::string key) {m_keyLong = std::move(key); return *this;}
    std::string keyLong() const {return m_keyLong;}

    Argument& setHelp(const std::string help) {m_help = std::move(help); return *this;}
    std::string help() const {return m_help;}

    template<typename T>
    Argument& setDefaultValue(const T value) {

        // assign default value
        m_valueDefault = std::move(value);

        // convert to string
        std::stringstream ss;
        ss << std::boolalpha << value;
        m_valueDefaultText = ss.str();

        return *this;
    }

    template<typename T>
    T getDefaultValue() const {

        if (m_valueDefault.has_value()) {
            return std::any_cast<T>(m_valueDefault);
        }

        throw std::logic_error("error, default value is missing.");
    }


    // set method for value
    Argument& setValue(const std::string value) {m_values.emplace_back(std::move(value)); return *this;}

    // get method for value
    template<typename T>
    std::enable_if_t<!is_stl_container<T>::value, T>
    get() const {

        if (m_values.empty())
            return getDefaultValue<T>();

        return getValue<T>(m_values.front());
    }

    // get method for container
    template<typename Cont>
    std::enable_if_t<is_stl_container<Cont>::value, Cont>
    get() const {
        using ContainerValueType = typename Cont::value_type;
        Cont container = {};

        if (m_values.empty()) {
            container.push_back(getDefaultValue<ContainerValueType>());
        }
        else {
            for (auto const &value: m_values) {
                container.push_back(getValue<ContainerValueType>(value));
            }
        }

        return container;
    }

    // format help message
    std::string printHelp(const std::size_t offset) const {
        auto margin = std::string(offset - printTag().size() + 4, ' ');
        return printTag() + margin + printInfo();
    }

private:
    bool m_used;
    Rank m_rank;
    int m_count;
    std::string m_name;
    std::string m_keyShort;
    std::string m_keyLong;
    std::string m_help;
    std::string m_valueDefaultText;
    std::any m_valueDefault;
    std::vector<std::string> m_values;

    // convert string value to basic type
    template<typename T>
    T getValue(const std::string &value) const {

        if (value.empty())
            return getDefaultValue<T>();

        T valueType;
        std::istringstream iss(value);
        iss >> std::noskipws >> valueType;
        if (!iss.eof() || iss.fail())
            throw std::logic_error("error, failed to convert value " + value);

        return valueType;
    }


    template<>
    std::string getValue(const std::string &value) const {
        if (value.empty())
            return getDefaultValue<std::string>();
        return value;
    }


    template<>
    bool getValue(const std::string &value) const {

        if (value.empty())
            return getDefaultValue<bool>();

        // convert to lower case
        std::string valueLC = "";
        for (char const &c: value)
            valueLC += static_cast<char>(std::tolower(c));

        bool valueType = false;
        std::istringstream iss(valueLC);
        iss >> std::noskipws >> valueType;
        if (iss.fail()) {
            iss.clear();
            iss >> std::boolalpha >> valueType;
        }

        if (!iss.eof() || iss.fail())
            throw std::logic_error("error, failed to convert value " + value);

        return valueType;
    }

    // print components
    std::string printDefaultValue() const {
        return m_valueDefaultText;
    }


    std::string printTag() const {

        if (keyShort().empty() && keyLong().empty())
            return name();

        auto tag = std::string("");
        if (!keyShort().empty())
            tag += keyShort();

        if (!keyLong().empty()) {
            if (!tag.empty())
                tag += ", ";
            tag += keyLong();
        }

        return tag;
    }


    std::string printInfo() const {

        auto info = std::string("");

        if ((rank() == Rank::COMMAND) || (rank() == Rank::POSITIONAL)) {
            info += help();
        }
        else {
            info += name();
            if (!help().empty())
                info += ", " + help();
        }


        if (!printDefaultValue().empty() && ((rank() == Rank::REQUIRED) || (rank() == Rank::OPTIONAL))) {
            if (!info.empty())
                info += ", ";
            info += "default [" + printDefaultValue() + "]";
        }


        if (rank() == Rank::REQUIRED) {
            if (!info.empty())
                info += ", ";
            info += "(required)";
        }

        return info;
    }

};



/**
 * @brief The ArgumentParser class
 *
 * @details maintains a container of all defined arguments,
 * provide a get<type> method to retrieve values based on
 * argument key or name. search in linear, but list of arguments if often short
 *
 */
class ArgumentParser
{
public:
    explicit ArgumentParser(const std::string name = "",
                            const std::string version = "0.0.1",
                            const std::string description = ""
                            ) :
        m_name(name),
        m_version(version),
        m_description(description)
    {}

    void setName(const std::string name) {m_name = std::move(name);}
    std::string name() const {return m_name;}

    void setVersion(const std::string version) {m_version = std::move(version);}
    std::string version() const {return m_version;}

    void setDescription(const std::string description) {m_description = std::move(description);}
    std::string description() const {return m_description;}

    Argument& addArgumentPositional(const std::string name) {
        return addArgument(name, Argument::Rank::POSITIONAL, m_args_positional);
    }

    Argument& addArgumentRequired(const std::string name) {
        return addArgument(name, Argument::Rank::REQUIRED, m_args_keyvalue);
    }

    Argument& addArgumentOptional(const std::string name) {
        return addArgument(name, Argument::Rank::OPTIONAL, m_args_keyvalue);
    }

    Argument& addArgumentFlag(const std::string name) {
        Argument& arg = addArgument(name, Argument::Rank::FLAG, m_args_keyvalue);
        arg.setDefaultValue<bool>(false).setCount(0);
        return arg;
    }

    Argument& addArgumentCommand(const std::string name) {
        Argument& arg = addArgument(name, Argument::Rank::COMMAND, m_args_command);
        arg.setDefaultValue<bool>(false).setCount(0);
        return arg;
    }


    /**
     * @brief parse
     * @details parse argc and argv to container of arguments
     * @param argc
     * @param argv
     */
    void parse(const int argc, const char *argv[]) {

        // build arguments map
        buildArgumentsMap(m_args_command);
        buildArgumentsMap(m_args_keyvalue);

        // assign name
        if (name().empty())
            parseName(std::string(argv[0]));

        // commands
        parseCommands(argc, argv);

        // parse
        parseArguments(argc, argv);

        // check for default flags
        if (checkDefaultFlag("-h", "--help"))
            throw std::logic_error(printHelp());

        if (checkDefaultFlag("-v", "--version"))
            throw std::logic_error(printVersion());

        // validate
        validateArguments(m_args_keyvalue);
        validateArguments(m_args_positional);
    }



    template<typename T>
    T get(const std::string &name) {

        std::shared_ptr<Argument> arg = findArgument(name);

        if (arg == nullptr)
            throw std::logic_error("error, argument with name=" + name + " is missing.\n\n" + printHelp());

        return arg->get<T>();
    }



    std::string printHelp() const {

        std::stringstream oss;

        oss << name() << std::endl;
        oss << printVersion() << std::endl;

        // Description
        if (!description().empty())
            oss << description() << std::endl;

        oss << std::endl;

        // Usage
        oss << printUsage() << std::endl << std::endl;

        // calculate margin
        std::size_t margin = 0;
        getArgumentMargin(margin, m_args_command);
        getArgumentMargin(margin, m_args_positional);
        getArgumentMargin(margin, m_args_keyvalue);


        // Commands
        if (!m_args_command.empty()) {
            oss << "[commands]" << std::endl;
            for (const auto &p_arg : m_args_command) {
                oss << "    " << p_arg->printHelp(margin) << std::endl;
            }
            oss << std::endl;
        }


        // Positional arguments
        if (!m_args_positional.empty()) {
            oss << "[arguments]" << std::endl;
            for (const auto &p_arg : m_args_positional) {
                oss << "    " << p_arg->printHelp(margin) << std::endl;
            }
            oss << std::endl;
        }

        // Required / Optional arguments
        if (!m_args_keyvalue.empty()) {
            oss << "[options]" << std::endl;
            for (const auto &p_arg : m_args_keyvalue) {
                oss << "    " << p_arg->printHelp(margin) << std::endl;
            }
            oss << std::endl;
        }

        return oss.str();
    }



    std::string printVersion() const {
        return "version " + m_version;
    }


    std::string printUsage() const {

        auto usage = std::string("");

        // add commands
        if (!m_args_command.empty())
            usage += " <command>";

        // add flags / optional / required
        if (!usage.empty())
            usage += " ";
        usage += printArgumentTag(Argument::Rank::FLAG);

        if (!usage.empty())
            usage += " ";
        usage += printArgumentTag(Argument::Rank::OPTIONAL);

        if (!usage.empty())
            usage += " ";
        usage += printArgumentTag(Argument::Rank::REQUIRED);


        // add positional
        for(const auto &p_arg : m_args_positional) {
            usage += " <" + p_arg->name() + ">";
        }

        return "usage: " + name() + usage;
    }

private:
    std::string m_name;
    std::string m_version;
    std::string m_description;
    std::vector<std::shared_ptr<Argument>> m_args_positional;
    std::vector<std::shared_ptr<Argument>> m_args_keyvalue;
    std::vector<std::shared_ptr<Argument>> m_args_command;
    std::unordered_map<std::string, std::shared_ptr<Argument>> m_args_map;


    Argument& addArgument(const std::string name,
                          const Argument::Rank rank,
                          std::vector<std::shared_ptr<Argument>> &args) {
        std::shared_ptr<Argument> arg = std::make_shared<Argument>();
        arg->setName(name);
        arg->setRank(rank);
        args.emplace_back(arg);
        return *arg;
    }


    // search arguments
    void buildArgumentsMap(std::vector<std::shared_ptr<Argument>> args) {
        for(std::vector<std::shared_ptr<Argument>>::iterator it = args.begin();
            it != args.end(); ++it) {

            if (!(*it)->name().empty())
                m_args_map[(*it)->name()] = (*it);

            if (!(*it)->keyShort().empty())
                m_args_map[(*it)->keyShort()] =  (*it);

            if (!(*it)->keyLong().empty())
                m_args_map[(*it)->keyLong()] = (*it);
        }
    }

    std::shared_ptr<Argument> findArgument(const std::string &key) {
        std::unordered_map<std::string, std::shared_ptr<Argument>>::iterator it = m_args_map.find(key);
        if (it == m_args_map.end())
            return nullptr;
        return it->second;
    }



    // validators
    bool hasHyphen(const std::string &opt, std::size_t index) {
        return opt.at(index) == '-';
    }


    bool hasEqualSiqn(const std::string &opt) {
        std::size_t found = opt.find('=');
        return found != std::string::npos;
    }


    bool hasLetter(const std::string &opt, std::size_t index) {
        return (('A' <= opt.at(index)) &&
                (opt.at(index) <= 'Z')) ||
               (('a' <= opt.at(index)) &&
                (opt.at(index) <= 'z'));
    }


    bool isShortKey(const std::string &opt) {
        if (opt.size() < 2)
            return false;
        return hasHyphen(opt, 0) && hasLetter(opt, 1);
    }


    bool isLongKey(const std::string &opt) {
        if (opt.size() < 4)
            return false;
        return hasHyphen(opt, 0) &&
               hasHyphen(opt, 1) &&
               hasLetter(opt, 2) &&
               hasLetter(opt, 3);
    }


    bool isCompositeFlag(const std::string &opt) {
        if(opt.size() < 3)
            return false;
        return isShortKey(opt) && !hasEqualSiqn(opt);
    }

    // parsers
    void parseName(const std::string &fullname) {

        // default name
        std::size_t posNameStart = 0;

        // remove path
        std::size_t posSlash = fullname.find_last_of("/");
        if (posSlash != std::string::npos)
            posNameStart = posSlash + 1;

        std::string name = fullname.substr(posNameStart);
        setName(name);
    }


    void parseKeyValuePair(const std::string &opt) {

        std::size_t found = opt.find('=');

        if (found == std::string::npos)
            throw std::invalid_argument("error, invalid key=value pair " + opt + "\n\n" + printHelp());

        auto key = opt.substr(0, found);
        auto value = opt.substr(found + 1, opt.size() - found);

        std::shared_ptr<Argument> arg = findArgument(key);
        if (arg == nullptr)
            throw std::invalid_argument("error, invalid argument " + key + "\n\n" + printHelp());

        arg->setValue(value).setUsed(true);
    }


    int parseKeyValueList(const std::string &opt, const int k, const int argc, const char *argv[]) {

        // find current argument
        std::shared_ptr<Argument> arg = findArgument(opt);
        if (arg == nullptr)
            throw std::invalid_argument("error, invalid argument " + opt + "\n\n" + printHelp());

        // assign as many values as required
        int steps = arg->count();
        if (steps == 0) {
            arg->setValue("1"); // used as a flag
        }
        else {

            // if steps are negative, switch to greedy mode
            // consume all values up to next key or last argument
            if (steps < 0)
                steps = argc - k;

            int i = k;
            while (i < (steps + k)) {
                auto value = std::string(argv[i]);

                if (isShortKey(value) || isLongKey(value))
                    break;

                arg->setValue(value);
                i++;
            }
            steps = i - k;
        }


        arg->setUsed(true);

        return steps;
    }


    void parseCompositeFlag(const std::string &opt) {

        for (std::size_t i = 1; i < opt.size(); ++i) {
            auto flag = std::string() + "-" + opt.at(i);
            std::shared_ptr<Argument> arg = findArgument(flag);
            if (arg == nullptr)
                throw std::invalid_argument("error, invalid argument " + opt + "\n\n" + printHelp());
            arg->setValue("1").setUsed(true);
        }

    }


    void parseCommands(const int argc, const char *argv[]) {

        if (m_args_command.empty())
            return;

        if (argc == 1)
            throw std::logic_error("error, command is required.\n\n" + printHelp());

        auto command = std::string(argv[1]);
        std::shared_ptr<Argument> arg = findArgument(command);

        if (arg == nullptr)
            throw std::invalid_argument("error, invalid command " + command + "\n\n" + printHelp());

        arg->setValue("1").setUsed(true);
    }


    void parseArguments(const int argc, const char *argv[]) {

        int k = 1;
        if (!m_args_command.empty())
            k = 2;

        std::vector<std::shared_ptr<Argument>>::iterator it = m_args_positional.begin();
        while (k < argc) {

            // current option
            auto opt = std::string(argv[k]);

            // update iterator
            k++;

            // parse long key --key value / --key=value
            if (isLongKey(opt)) {

                if (hasEqualSiqn(opt)) {
                    parseKeyValuePair(opt);
                    continue;
                }
                else {
                    int steps = parseKeyValueList(opt, k, argc, argv);
                    k += steps;
                    continue;
                }
            }

            // parse shor key -k value / -k=value
            // parse composite flags -abc {no value}
            if (isShortKey(opt)) {

                if (hasEqualSiqn(opt)) {
                    parseKeyValuePair(opt);
                    continue;
                }
                else {

                    if (isCompositeFlag(opt)) {
                        parseCompositeFlag(opt);
                        continue;
                    }
                    else {
                        int steps = parseKeyValueList(opt, k, argc, argv);
                        k += steps;
                        continue;
                    }
                }
            }

            // assign positional argument
            if (it != m_args_positional.end()) {
                (*it)->setValue(opt).setUsed(true);
                ++it;
                continue;
            }

            // if reached throw unknown option
            throw std::invalid_argument("error, unknown argument " + opt + "\n\n" + printHelp());
        }
    }


    void validateArguments(const std::vector<std::shared_ptr<Argument>> &args) {

        for(std::vector<std::shared_ptr<Argument>>::const_iterator it = args.begin();
            it != args.end(); ++it) {
            if (((*it)->rank() == Argument::Rank::REQUIRED) &&
                ((*it)->used() == false))
                throw std::invalid_argument("error, required argument " + (*it)->name() + " is unused.\n\n" + printHelp());
        }
    }


    bool checkDefaultFlag(const std::string &keyShort, const std::string &keyLong) {

        // check key Short
        std::shared_ptr<Argument> arg = findArgument(keyShort);
        if (arg != nullptr)
            return arg->used();

        // check key Long
        arg = findArgument(keyLong);
        if (arg != nullptr)
            return arg->used();

        return false;
    }


    std::string printArgumentTag(const Argument::Rank rank) const {

        auto argumentTag = std::string("");

        for (auto const &p_arg : m_args_keyvalue) {

            auto tag = std::string("");
            if (!p_arg->keyShort().empty())
                tag += p_arg->keyShort();

            if (!p_arg->keyLong().empty()) {
                if (!tag.empty())
                    tag += "|";
                tag += p_arg->keyLong();
            }

            auto value = std::string("");
            if (p_arg->count() > 0)
                value = " <value>";

            if (p_arg->count() < 0)
                value = " <value>...";

            if ((p_arg->rank() == rank) && (!tag.empty())) {
                if (!argumentTag.empty())
                    argumentTag += " ";
                argumentTag += tag +  value;
            }

        }

        if (!argumentTag.empty())
            argumentTag = "[" + argumentTag + "]";

        return argumentTag;
    }


    void getArgumentMargin(std::size_t &margin, const std::vector<std::shared_ptr<Argument>> &args) const {

        for (const auto &ptr : args) {
            std::size_t value = ptr->printTag().size();
            if (margin < value)
                margin = value;
        }

    }

};
#endif // ARGUMENTPARSER_H
