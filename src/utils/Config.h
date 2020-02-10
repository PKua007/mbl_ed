//--------------------------------------------------------------------------------------------
// Simple configuration file parser
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#ifndef _CONFIG_H
#define _CONFIG_H


#include <map>
#include <stdexcept>
#include <iosfwd>
#include <vector>
#include <algorithm>


/**
 * @brief An exception thrown if there was a problem parsing config in Config.
 */
class ConfigParseException : public std::runtime_error
{
public:
    explicit ConfigParseException(const std::string & _what) : std::runtime_error(_what)
    {}
};


/**
 * @brief An exception thrown if there was an access to nonexistent field in Config.
 */
class ConfigNoFieldException : public std::runtime_error
{
public:
    explicit ConfigNoFieldException(const std::string & _what) : std::runtime_error(_what)
    {}
};


/**
 * @brief A key=value config file parser (INI format) with sections support.
 *
 * File format:
 * \code
 * key1=value1
 * # standalone comment
 * key2=value2 # end-line comment;
 *
 * # empty lines are omitted
 * key3 = value3 # whitespace is trimmed
 * # key3 = value3 - duplicate fields are forbidden
 * \endcode
 */
class Config
{
private:
    std::map<std::string, std::string>  fieldMap;
    std::vector<std::string>            keys;
    std::vector<std::string>            rootSections;

    struct Field {
        std::string key;
        std::string value;
    };

    static void stripComment(std::string &line);
    static bool isItSectionEntry(std::string &line, std::size_t lineNum);
    static Field splitField(const std::string &line, char delim, std::size_t lineNum,
                            const std::string &currentSection);

    void buildRootSections();

    Config() = default;

public:

    /**
     * @brief Parses given stream.
     *
     * Format:
     * \code
     * key1=value1
     * # standalone comment
     * ; semicolon is also supported
     * key2=value2 # end-line comment
     * [section]
     * key3 = value3
     *
     * # empty lines are omitted
     * key4 = value4 # whitespace is trimmed
     * \endcode
     * The keys are @a key1, @a key2, @a section.key3, @a section.key4.
     *
     * @param in stream to parse from
     * @param delim delimiter for key, value; defaults to '='
     * @param allowRedefinition if `true`, redefinition of field will overwrite the old value; if `false`, it will
     * throw ConfigParseException
     * @throws ConfigParseException on parse error (no delimiter of a duplicate field)
     * @throws std::invalid_argument when delim = '#' (comment)
     * @return Config object to be deleted manualy after use
     */
    static Config parse(std::istream &in, char delim = '=', bool allowRedefinition = false);

    bool hasField(const std::string &field) const;
    std::size_t size() const { return this->keys.size(); }
    bool empty() const { return this->keys.empty(); }

    std::string getString(const std::string &field) const;
    int getInt(const std::string &field) const;
    unsigned long getUnsignedLong(const std::string &field) const;
    double getDouble(const std::string &field) const;
    float getFloat(const std::string &field) const;
    bool getBoolean(const std::string &field) const;

    /**
     * Returns keys in a config, preserving order from an input
     * @return keys in a config, preserving order from an input
     */
    std::vector<std::string> getKeys() const { return this->keys; }

    std::vector<std::string> getRootSections() const { return this->rootSections; }
    bool hasRootSection(const std::string &section) const;

    Config fetchSubconfig(const std::string &rootSection) const;
};


#endif // _CONFIG_H