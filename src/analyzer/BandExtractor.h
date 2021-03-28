//
// Created by Piotr Kubala on 28/03/2021.
//

#ifndef MBL_ED_BANDEXTRACTOR_H
#define MBL_ED_BANDEXTRACTOR_H

#include <utility>
#include <variant>
#include <vector>

#include "core/FockBasis.h"
#include "core/Eigensystem.h"
#include "utils/Logger.h"

class BandExtractor {
public:
    using WidthMargin = double;
    using NumberOfEnergiesMargin = std::size_t;

    using Margin = std::variant<WidthMargin, NumberOfEnergiesMargin>;

    /**
     * @brief Range in the epsilon version - band of energies normalized to [0, 1]
     */
    struct EpsilonRange {
        EpsilonRange() = default;
        EpsilonRange(double epsilonMiddle, Margin epsilonMargin);

        double epsilonMiddle{};
        Margin epsilonMargin{};
    };

    /**
     * @brief Range in the vector version - band of energies normalized to [0, 1], middle given by @a middleVector
     */
    struct VectorRange {
        VectorRange() = default;
        VectorRange(FockBasis::Vector middleVector, Margin epsilonMargin);

        FockBasis::Vector middleVector{};
        Margin epsilonMargin{};
    };

    /**
     * @brief Rande in the CDF version - energies given by quantiles determined by @a cdfMiddle and @a cdfMargin
     */
    struct CDFRange {
        CDFRange() = default;
        CDFRange(double cdfMiddle, double cdfMargin);

        double cdfMiddle{};
        double cdfMargin{};
    };

    using Range = std::variant<EpsilonRange, VectorRange, CDFRange>;

private:
    Range range;
    std::string taskName;

    static double calculateEnergyOfFockState(const FockBasis::Vector &state, const Eigensystem &eigensystem);

    std::vector<std::size_t> getIndicesForVectorRange(const Eigensystem &eigensystem,
                                                      const arma::vec &normalizedEnergies,
                                                      const VectorRange &vectorRange, Logger &logger) const;

    std::vector<std::size_t> getIndicesForEpsilonRange(const Eigensystem &eigensystem,
                                                       const arma::vec &normalizedEnergies,
                                                       const EpsilonRange &epsilonRange, Logger &logger) const;

    std::vector<std::size_t> getIndicesForVectorRange(const Eigensystem &eigensystem,
                                                      const arma::vec &normalizedEnergies,
                                                      const CDFRange &cdfRange, Logger &logger) const;

public:
    explicit BandExtractor(Range range, std::string taskName = "Task")
            : range{std::move(range)}, taskName{std::move(taskName)}
    { }

    std::vector<std::size_t> getBandIndices(const Eigensystem &eigensystem,  Logger &logger) const;
};


#endif //MBL_ED_BANDEXTRACTOR_H
