//
// Created by pkua on 15.05.19.
//

#ifndef RSA3D_ASSERTIONS_H
#define RSA3D_ASSERTIONS_H

#include <stdexcept>

// Cpp Core Guidelines-style assertions for design by contract
// https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#i6-prefer-expects-for-expressing-preconditions

#define ASSERTIONS_S(x) #x
#define ASSERTIONS_S_(x) ASSERTIONS_S(x)
#define __WHERE__ __FILE__ ":" ASSERTIONS_S_(__LINE__)

// Preconditions check (argument validation)
#define Expects(cond) if (!(cond)) throw std::invalid_argument(__WHERE__ ": Precondition " #cond " failed")

// Postconditions check (results assertion)
#define Ensures(cond) if (!(cond)) throw std::runtime_error(__WHERE__ ": Postcondition " #cond " failed")

// Runtime assertion. Why duplicate assert from cassert? Because we don't want to disable is in release mode and
// be more C++ and throw exception
#define Assert(cond) if (!(cond)) throw std::runtime_error(__WHERE__ ": Assertion " #cond " failed")

// Additional macros for validating things like input from file - wrong input shouldn't be considered as assertion
// fail, because it is not programmer's fault ;)

#define Validate(cond) if (!(cond)) throw ValidationException(__WHERE__ ": Validation " #cond " failed")
#define ValidateMsg(cond, msg) if (!(cond)) throw ValidationException(msg)

/**
 * @brief An exceptions thrown by Validate and ValidateMsg macros.
 */
struct ValidationException : public std::domain_error {
    explicit ValidationException(const std::string &msg) : std::domain_error{msg} { }
};

#endif //RSA3D_ASSERTIONS_H
