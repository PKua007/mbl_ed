//
// Created by Piotr Kubala on 05/09/2020.
//

#ifndef MBL_ED_RESTORABLEHELPER_H
#define MBL_ED_RESTORABLEHELPER_H

#include <vector>
#include <ostream>
#include <istream>
#include <type_traits>

#include "utils/Assertions.h"
#include "Restorable.h"

class RestorableHelper {
public:
    template<typename T>
    static void storeStateForVector(const std::vector<T> &vector, std::ostream &binaryOut) {
        std::size_t vectorSize = vector.size();
        binaryOut.write(reinterpret_cast<const char*>(&vectorSize), sizeof(vectorSize));
        binaryOut.write(reinterpret_cast<const char*>(vector.data()), sizeof(vector[0]) * vectorSize);
        Assert(binaryOut.good());
    }

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

    template<typename T>
    static void storeStateForStaticRestorableVector(const std::vector<T> &vector, std::ostream &binaryOut) {
        static_assert(std::is_base_of_v<Restorable, T>);

        std::size_t vectorSize = vector.size();
        binaryOut.write(reinterpret_cast<const char*>(&vectorSize), sizeof(vectorSize));
        Assert(binaryOut.good());
        for (const auto &entry : vector)
            entry.storeState(binaryOut);
    }

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
