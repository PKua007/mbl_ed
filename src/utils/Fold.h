//
// Created by Piotr Kubala on 21/01/2020.
//

#ifndef MBL_ED_FOLD_H
#define MBL_ED_FOLD_H

#include <sstream>
#include <utility>

class Fold {
private:
    std::size_t width_{};
    std::size_t margin_{};
    std::string str;

    [[nodiscard]] std::size_t findEnterInCurrentLine(std::size_t lineBeg, std::size_t effectiveWidth) const;
    [[nodiscard]] std::size_t findLastSpaceInCurrentLine(std::size_t lineBeg, std::size_t effectiveWidth) const;
    [[nodiscard]] bool spaceOnBreak(std::size_t lineBeg, std::size_t effectiveWidth) const;
    [[nodiscard]] bool notLastLine(std::size_t lineBeg, std::size_t effectiveWidth) const;

public:
    explicit Fold(std::string str) : str{std::move(str)} { }
    Fold() : str{""} { }

    Fold &width(std::size_t width_);
    Fold &margin(std::size_t margin_);
    operator std::string() const;
};

std::ostream &operator<<(std::ostream &out, const Fold &fold);


#endif //MBL_ED_FOLD_H
