//
// Created by Piotr Kubala on 31/08/2020.
//

#ifndef MBL_ED_RESTORABLE_H
#define MBL_ED_RESTORABLE_H

#include <iosfwd>

class Restorable {
public:
    virtual ~Restorable() = default;

    virtual void storeState(std::ostream &binaryOut) const = 0;
    virtual void joinRestoredState(std::istream &binaryIn) = 0;
    virtual void clear() = 0;

    void restoreState(std::istream &binaryIn) {
        this->clear();
        this->joinRestoredState(binaryIn);
    };

};


#endif //MBL_ED_RESTORABLE_H
