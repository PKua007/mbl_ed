//
// Created by Piotr Kubala on 16/01/2020.
//

#ifndef MBL_ED_FILEUTILS_H
#define MBL_ED_FILEUTILS_H

#include <fstream>
#include <memory>

/**
 * @brief Exception thrown when FileOstreamProvider or FileIstreamProvider fail to open the file.
 */
class FileException : public std::runtime_error {
public:
    explicit FileException(const std::string &what) : std::runtime_error(what) { }
};

/**
 * @brief A class, which provides std::ostream from a file with given name.
 *
 * Virtual openFile function enables to mock this class.
 */
class FileOstreamProvider {
private:
    std::string fileDescription = "";

public:
    virtual ~FileOstreamProvider() { }

    /**
     * @brief Returns current file description.
     * @return current file description
     */
    const std::string& getFileDescription() const { return fileDescription; }

    /**
     * @brief Sets new file description.
     * @param fileDescription new file description
     */
    void setFileDescription(const std::string &fileDescription) { this->fileDescription = fileDescription; }

    /**
     * @brief Opens file with name @a filename to read.
     *
     * If it fails it throws FileException with appropriate message constructed with getFileDescription.
     *
     * @return std::ostream from the file
     */
    virtual std::unique_ptr<std::ostream> openFile(const std::string &filename) const {
        auto file = std::unique_ptr<std::ofstream>(new std::ofstream(filename));
        if (!(*file))
            throw FileException("Cannot open file " + filename + " to read: " + this->fileDescription);
        return file;
    }
};

/**
 * @brief A class, which provides std::istream from a file with given name.
 *
 * Virtual openFile function enables to mock this class.
 */
class FileIstreamProvider {
private:
    std::string fileDescription = "";

public:
    virtual ~FileIstreamProvider() { }

    /**
     * @brief Returns current file description.
     * @return current file description
     */
    const std::string& getFileDescription() const { return fileDescription; }

    /**
     * @brief Sets new file description.
     * @param fileDescription new file description
     */
    void setFileDescription(const std::string &fileDescription) { this->fileDescription = fileDescription; }

    /**
     * @brief Opens file with name @a filename to write.
     *
     * If it fails it throws FileException with appropriate message constructed with getFileDescription.
     *
     * @return std::istream from the file
     */
    virtual std::unique_ptr<std::istream> openFile(const std::string &filename) const {
        auto file = std::unique_ptr<std::ifstream>(new std::ifstream(filename));
        if (!(*file))
            throw FileException("Cannot open file " + filename + " to write: " + this->fileDescription);
        return file;
    }
};

#endif //MBL_ED_FILEUTILS_H
