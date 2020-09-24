//
// Created by Piotr Kubala on 05/09/2020.
//

#ifndef MBL_ED_RESTORABLEHELPER_H
#define MBL_ED_RESTORABLEHELPER_H

#include <vector>
#include <ostream>
#include <istream>
#include <type_traits>
#include <memory>

#include "utils/Assertions.h"
#include "Restorable.h"

/**
 * @brief A set of helper function for storing/restoring selected types of data for Restorable classes.
 */
class RestorableHelper {
private:
    template<typename T> struct is_smart_ptr : std::false_type {};
    template<typename T> struct is_smart_ptr<std::shared_ptr<T>> : std::true_type {};
    template<typename T> struct is_smart_ptr<std::unique_ptr<T>> : std::true_type {};
    template<typename T> inline static constexpr bool is_smart_ptr_v = is_smart_ptr<T>::value;

public:
    /**
     * @brief Restorable::storeState for vector data.
     */
    template<typename T>
    static void storeStateForVector(const std::vector<T> &vector, std::ostream &binaryOut) {
        std::size_t vectorSize = vector.size();
        binaryOut.write(reinterpret_cast<const char*>(&vectorSize), sizeof(vectorSize));
        binaryOut.write(reinterpret_cast<const char*>(vector.data()), sizeof(vector[0]) * vectorSize);
        Assert(binaryOut.good());
    }

    /**
     * @brief Restorable::joinRestoredState for vector data.
     * @details Restored vector is concatenated with @a vector.
     */
    template<typename T>
    static void joinRestoredStateForVector(std::vector<T> &vector, std::istream &binaryIn) {
        std::size_t vectorSizeRestored{};
        binaryIn.read(reinterpret_cast<char*>(&vectorSizeRestored), sizeof(vectorSizeRestored));
        Assert(binaryIn.good());

        std::vector<T> vectorRestored(vectorSizeRestored);
        binaryIn.read(reinterpret_cast<char*>(vectorRestored.data()), sizeof(vectorRestored[0]) * vectorSizeRestored);
        Assert(binaryIn.good());

        vector.reserve(vector.size() + vectorSizeRestored);
        vector.insert(vector.end(), vectorRestored.begin(), vectorRestored.end());
    }

    /**
     * @brief Restorable::storeState for vector data with static size containing other Restorables.
     */
    template<typename T>
    static void storeStateForStaticRestorableVector(const std::vector<T> &vector, std::ostream &binaryOut) {
        static_assert(std::is_base_of_v<Restorable, T>);

        std::size_t vectorSize = vector.size();
        binaryOut.write(reinterpret_cast<const char*>(&vectorSize), sizeof(vectorSize));
        Assert(binaryOut.good());
        for (const auto &entry : vector)
            entry.storeState(binaryOut);
    }

    /**
     * @brief Restorable::joinRestoredState for vector data with static size containing other Restorables.
     */
    template<typename T>
    static void joinRestoredStateForStaticRestorableVector(std::vector<T> &vector, std::istream &binaryIn) {
        static_assert(std::is_base_of_v<Restorable, T>);

        std::size_t vectorSizeRestored{};
        binaryIn.read(reinterpret_cast<char*>(&vectorSizeRestored), sizeof(vectorSizeRestored));
        Assert(binaryIn.good());
        Assert(vectorSizeRestored == vector.size());

        for (auto &entry : vector)
            entry.joinRestoredState(binaryIn);
    }

    /**
     * @brief Restorable::storeState for vector data with static size containing other Restorables as a smart pointer.
     */
    template<typename T, template<typename> typename SmartPointer>
    static void storeStateForStaticRestorableVector(const std::vector<SmartPointer<T>> &vector,
                                                    std::ostream &binaryOut)
    {
        static_assert(std::is_base_of_v<Restorable, T>);
        static_assert(is_smart_ptr_v<std::remove_cv_t<SmartPointer<T>>>);

        std::size_t vectorSize = vector.size();
        binaryOut.write(reinterpret_cast<const char*>(&vectorSize), sizeof(vectorSize));
        Assert(binaryOut.good());
        for (const auto &entry : vector)
            entry->storeState(binaryOut);
    }

    /**
     * @brief Restorable::joinRestoredState for vector data with static size containing other Restorables as a smart
     * pointer.
     */
    template<typename T, template<typename> typename SmartPointer>
    static void joinRestoredStateForStaticRestorableVector(const std::vector<SmartPointer<T>> &vector,
                                                           std::istream &binaryIn)
    {
        static_assert(std::is_base_of_v<Restorable, T>);
        static_assert(is_smart_ptr_v<std::remove_cv_t<SmartPointer<T>>>);

        std::size_t vectorSizeRestored{};
        binaryIn.read(reinterpret_cast<char*>(&vectorSizeRestored), sizeof(vectorSizeRestored));
        Assert(binaryIn.good());
        Assert(vectorSizeRestored == vector.size());

        for (auto &entry : vector)
            entry->joinRestoredState(binaryIn);
    }

    /**
     * @brief Restorable::storeState for a histogram, so a vector of static size containg bins (vectors) with dynamic
     * size.
     */
    template<typename T>
    static void storeStateForHistogram(const std::vector<std::vector<T>> &histogram, std::ostream &binaryOut) {
        std::size_t numBins = histogram.size();
        binaryOut.write(reinterpret_cast<const char*>(&numBins), sizeof(numBins));
        Assert(binaryOut.good());

        for (const auto &binEntries : histogram) {
            std::size_t numEntries = binEntries.size();
            binaryOut.write(reinterpret_cast<const char*>(&numEntries), sizeof(numEntries));
            binaryOut.write(reinterpret_cast<const char*>(binEntries.data()), sizeof(binEntries[0]) * numEntries);
            Assert(binaryOut.good());
        }
    }

    /**
     * @brief Restorable::joinRestoredState for a histogram, so a vector of static size containg bins (vectors) with
     * dynamic size.
     * @details Loaded histogram is asserted to have the same number of bins and @a histogram and both are concatenated.
     */
    template<typename T>
    static void joinRestoredStateForHistogram(std::vector<std::vector<T>> &histogram, std::istream &binaryIn) {
        std::size_t numBinsRestored{};
        binaryIn.read(reinterpret_cast<char*>(&numBinsRestored), sizeof(numBinsRestored));
        Assert(binaryIn.good());
        Assert(numBinsRestored == histogram.size());

        for (std::size_t i{}; i < numBinsRestored; i++) {
            std::size_t numEntriesRestored{};
            binaryIn.read(reinterpret_cast<char*>(&numEntriesRestored), sizeof(numEntriesRestored));
            Assert(binaryIn.good());

            std::vector<T> entriesRestored(numEntriesRestored);
            binaryIn.read(reinterpret_cast<char*>(entriesRestored.data()),
                          sizeof(entriesRestored[0]) * numEntriesRestored);
            Assert(binaryIn.good());

            std::vector<T> &binEntries = histogram[i];
            binEntries.reserve(binEntries.size() + numEntriesRestored);
            binEntries.insert(binEntries.end(), entriesRestored.begin(), entriesRestored.end());
        }
    }
};


#endif //MBL_ED_RESTORABLEHELPER_H
