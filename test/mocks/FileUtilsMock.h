//
// Created by Piotr Kubala on 17/01/2020.
//

#ifndef MBL_ED_FILEUTILSMOCK_H
#define MBL_ED_FILEUTILSMOCK_H

#include <catch2/trompeloeil.hpp>

#include "utils/FileUtils.h"

class FileOstreamProviderMock : public FileOstreamProvider {
public:
    MAKE_CONST_MOCK2(openOutputFile, std::unique_ptr<std::ostream>(const std::string &, bool), override);
};

class FilesystemManipulatorMock : public FilesystemManipulator {
public:
    MAKE_CONST_MOCK2(openOutputFile, std::unique_ptr<std::ostream>(const std::string &, bool), override);
    MAKE_CONST_MOCK2(openInputFile, std::unique_ptr<std::istream>(const std::string &, bool), override);
    MAKE_CONST_MOCK1(listFilesInDirectory, std::vector<std::filesystem::path>(const std::string &), override);
    MAKE_CONST_MOCK1(deleteFile, void(const std::string &), override);
};

#endif //MBL_ED_FILEUTILSMOCK_H
