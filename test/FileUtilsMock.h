//
// Created by Piotr Kubala on 17/01/2020.
//

#ifndef MBL_ED_FILEUTILSMOCK_H
#define MBL_ED_FILEUTILSMOCK_H

#include <catch2/trompeloeil.hpp>

#include "FileUtils.h"

class FileOstreamProviderMock : public FileOstreamProvider {
public:
    MAKE_CONST_MOCK1(openFile, std::unique_ptr<std::ostream>(const std::string &), override);
};

#endif //MBL_ED_FILEUTILSMOCK_H
