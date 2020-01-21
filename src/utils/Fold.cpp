//
// Created by Piotr Kubala on 21/01/2020.
//

#include "Fold.h"
#include "Assertions.h"

Fold::operator std::string() const {
    Expects(width_ > 1);
    Expects( margin_ < width_);
    if (this->str.empty())
        return "";

    std::ostringstream formatted;
    std::size_t lineBeg{};

    std::string marginSpaces = std::string(this->margin_, ' ');
    while (lineBeg + this->width_ - this->margin_ - 2 < this->str.size()) {
        std::size_t newline = this->str.find('\n', lineBeg);
        if (newline != std::string::npos && newline - lineBeg < this->width_ - this->margin_) {
            formatted << marginSpaces << this->str.substr(lineBeg, newline - lineBeg + 1);
            lineBeg = newline + 1;
            continue;
        }

        if (lineBeg + this->width_ - this->margin_ < this->str.size() && this->str[lineBeg + this->width_ - this->margin_] == ' ') {
            formatted << marginSpaces << this->str.substr(lineBeg, this->width_ - this->margin_) << '\n';
            lineBeg = lineBeg + this->width_ - this->margin_ + 1;
            continue;
        }

        std::size_t spacePos = this->str.rfind(' ', lineBeg + this->width_ - this->margin_);
        if (spacePos == std::string::npos || spacePos < lineBeg) {
            formatted << marginSpaces << this->str.substr(lineBeg, this->width_ - this->margin_) << '\n';
            lineBeg = lineBeg + this->width_ - this->margin_;
        } else {
            formatted << marginSpaces << this->str.substr(lineBeg, spacePos - lineBeg + 1) << '\n';
            lineBeg = spacePos + 1;
        }
    }

    if (lineBeg < this->str.size())
        formatted << marginSpaces << this->str.substr(lineBeg);
    return formatted.str();
}

Fold &Fold::width(std::size_t width_) {
    this->width_ = width_;
    return *this;
}

Fold &Fold::margin(std::size_t margin_) {
    this->margin_ = margin_;
    return *this;
}

std::ostream &operator<<(std::ostream &out, const Fold &fold) {
    return out << static_cast<std::string>(fold);
}
