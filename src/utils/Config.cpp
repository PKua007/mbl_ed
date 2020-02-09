//--------------------------------------------------------------------------------------------
// Simple configuration file parser
//--------------------------------------------------------------------------------------------
// (C)PKua 2017
//--------------------------------------------------------------------------------------------

#include "Config.h"
#include "Utils.h"
#include "Assertions.h"

#include <istream>


Config Config::parse(std::istream &in, char delim, bool allowRedefinition) {
    Expects(delim != '#' && delim != ';');

    Config result;
    std::size_t lineNum = 0;
    std::string line;
    std::string currentSection;
    while (std::getline(in, line)) {
        lineNum++;
        stripComment(line);
        trim(line);
        if (line.empty())
            continue;

        if (isItSectionEntry(line, lineNum)) {
            currentSection = line.substr(1, line.size() - 2) + ".";
            continue;
        }

        auto field = splitField(line, delim, lineNum, currentSection);
        if (result.hasParam(field.key) && !allowRedefinition) {
            throw ConfigParseException("Redefinition of field \"" + field.key + "\" in line "
                                       + std::to_string(lineNum));
        }

        result.fieldMap[field.key] = field.value;
        if (!result.hasParam(field.key))
            result.keys.push_back(field.key);
    }

    return result;
}

bool Config::isItSectionEntry(std::string &line, std::size_t lineNum) {
    if (line.front() != '[')
        return false;
    if (line.back() != ']')
        throw ConfigParseException("Malformed section entry in line " + std::to_string(lineNum));
    return line.front() == '[';
}

Config::Field Config::splitField(const std::string &line, char delim, std::size_t line_num,
                                 const std::string &currentSection)
{
    Field keyValue;
    std::size_t pos = line.find(delim);
    if (pos == std::string::npos)
        throw ConfigParseException("No '" + std::string(1, delim) + "' sign in line " + std::to_string(line_num));

    keyValue.key = currentSection + line.substr(0, pos);
    trim(keyValue.key);
    keyValue.value = (pos == line.length() - 1) ? "" : line.substr(pos + 1);
    trim(keyValue.value);
    return keyValue;
}

void Config::stripComment(std::string &line) {
    std::size_t pos = line.find_first_of(";#");
    if (pos != std::string::npos)
        line.erase(pos);
}

std::string Config::getString(const std::string & field) const {
    auto iter = this->fieldMap.find(field);
    if (iter == this->fieldMap.end())
        throw ConfigNoFieldException("No \"" + field + "\" field in config");
    return (*iter).second;
}

int Config::getInt(const std::string & field) const {
    return std::stoi(this->getString(field));
}

unsigned long Config::getUnsignedLong(const std::string & field) const {
    auto str = this->getString(field);
    if (std::stoi(str) < 0)
        throw std::invalid_argument("unsigned long field negative");
    return std::stoul(str);
}

double Config::getDouble(const std::string & field) const {
    return std::stod(this->getString(field));
}

float Config::getFloat(const std::string & field) const
{
    return std::stof(this->getString(field));
}

bool Config::getBoolean(const std::string &field) const {
    if (this->getString(field) == "true")
        return true;
    else if (this->getString(field) == "false")
        return false;
    else
        throw ConfigParseException("Cannot read boolean from " + field + ". It is neither 'true' nor 'false'");
}
