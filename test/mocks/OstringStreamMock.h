//
// Created by Piotr Kubala on 17/01/2020.
//

#ifndef MBL_ED_OSTRINGSTREAMMOCK_H
#define MBL_ED_OSTRINGSTREAMMOCK_H

#include "../../extern/Catch2/single_include/catch2/catch.hpp"
#include <sstream>

class OstringStreamMock : public std::ostringstream {
private:
    std::string expectedContentOnDestruction;
public:
    explicit OstringStreamMock(std::string expectedContentOnDestruction)
            : expectedContentOnDestruction{std::move(expectedContentOnDestruction)}
    { }

    ~OstringStreamMock() override {
        REQUIRE(this->str() == this->expectedContentOnDestruction);
    }
};

#endif //MBL_ED_OSTRINGSTREAMMOCK_H
