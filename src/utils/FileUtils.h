//
// Created by Piotr Kubala on 16/01/2020.
//

#ifndef MBL_ED_FILEUTILS_H
#define MBL_ED_FILEUTILS_H

#include <fstream>
#include <memory>
#include <filesystem>

/**
 * @brief Exception thrown when FileOstreamProvider or FileIstreamProvider fail to open the file.
 */
class FileException : public std::runtime_error {
public:
    explicit FileException(const std::string &what) : std::runtime_error(what) { }
};

/**
 * @brief Base class for FileOstreamProvider and FileIstreamProvider.
 */
class FileStreamProvider {
private:
    std::string fileDescription = "";

public:
    /**
     * @brief Returns current file description.
     * @return current file description
     */
    [[nodiscard]] const std::string& getFileDescription() const { return fileDescription; }

    /**
     * @brief Sets new file description.
     * @param fileDescription_ new file description
     */
    void setFileDescription(const std::string &fileDescription_) { this->fileDescription = fileDescription_; }
};

/**
 * @brief A class, which provides std::ostream from a file with given name.
 * @details Virtual openOutputFile function enables to mock this class.
 */
class FileOstreamProvider : virtual public FileStreamProvider {
public:
    virtual ~FileOstreamProvider() = default;

    /**
     * @brief Opens file with name @a filename to read (in binary, if `binary == true`).
     * @details If it fails it throws FileException with appropriate message constructed with getFileDescription.
     * @return std::ostream from the file
     */
    [[nodiscard]] virtual std::unique_ptr<std::ostream> openOutputFile(const std::string &filename, bool binary) const {
        std::unique_ptr<std::ofstream> file;
        if (binary)
            file = std::make_unique<std::ofstream>(filename, std::ios::out | std::ios::binary);
        else
            file = std::make_unique<std::ofstream>(filename);
        if (!(*file))
            throw FileException("Cannot open file " + filename + " to read: " + this->getFileDescription());
        return file;
    }

    /**
     * @brief The same as openOutputFile(const std::string&, bool), with with `binary == false`.
     */
    [[nodiscard]] std::unique_ptr<std::ostream> openOutputFile(const std::string &filename) const {
        return this->openOutputFile(filename, false);
    }
};

/**
 * @brief A class, which provides std::istream from a file with given name.
 * @details Virtual openOutputFile function enables to mock this class.
 */
class FileIstreamProvider : virtual public FileStreamProvider {
public:
    virtual ~FileIstreamProvider() = default;

    /**
     * @brief Opens file with name @a filename to write (in binary, if `binary == true`).
     * @details If it fails it throws FileException with appropriate message constructed with getFileDescription.
     * @return std::istream from the file
     */
    [[nodiscard]] virtual std::unique_ptr<std::istream> openInputFile(const std::string &filename, bool binary) const {
        std::unique_ptr<std::ifstream> file;
        if (binary)
            file = std::make_unique<std::ifstream>(filename, std::ios::in | std::ios::binary);
        else
            file = std::make_unique<std::ifstream>(filename);
        if (!(*file))
            throw FileException("Cannot open file " + filename + " to write: " + this->getFileDescription());
        return file;
    }

    /**
     * @brief The same as openInputFile(const std::string&, bool), with with `binary == false`.
     */
    [[nodiscard]] std::unique_ptr<std::istream> openInputFile(const std::string &filename) const {
        return this->openInputFile(filename, false);
    }
};

/**
 * @brief Combines FileIstreamProvider, FileOstreamProvider and adds listing files in directory and removing files.
 * @details Virtual methods enable mocking.
 */
class FilesystemManipulator : public FileIstreamProvider, public FileOstreamProvider {
public:
    /**
     * @brief Returns list of std::filesystem::paths of all non-directory fiels in a given directory.
     */
    [[nodiscard]] virtual std::vector<std::filesystem::path> listFilesInDirectory(const std::string &directory) const {
        std::vector<std::filesystem::path> files;
        for (const auto &entry : std::filesystem::directory_iterator(directory)) {
            if (!entry.is_directory())
                files.push_back(entry.path());
        }
        return files;
    }

    /**
     * @brief Deletes a file of a given name.
     */
    virtual void deleteFile(const std::string &filename) const {
        std::filesystem::remove(filename);
    }
};

#endif //MBL_ED_FILEUTILS_H
