//
// Created by Piotr Kubala on 31/08/2020.
//

#ifndef MBL_ED_RESTORABLE_H
#define MBL_ED_RESTORABLE_H

#include <iosfwd>

/**
 * @brief A class representing a set of data, which can be stored, restored and joined during restoration.
 */
class Restorable {
public:
    virtual ~Restorable() = default;

    /**
     * @brief Stores the state of the class in binary to @a binaryOut.
     */
    virtual void storeState(std::ostream &binaryOut) const = 0;

    /**
     * @brief Restores data stored in binary in @a binaryIn and joins this restored data with the data already
     * present in the class.
     * @details For example, if the class stores a set of energies, the restored energies should be appended to this
     * set.
     */
    virtual void joinRestoredState(std::istream &binaryIn) = 0;

    /**
     * @brief Resets the data.
     */
    virtual void clear() = 0;

    /**
     * @brief Overwrites existing data with one restored from @a binaryIn.
     */
    void restoreState(std::istream &binaryIn) {
        this->clear();
        this->joinRestoredState(binaryIn);
    };

};


#endif //MBL_ED_RESTORABLE_H
