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

/**
 * @brief A helper class performing extraction of energies from a selected band.
 * @details The band can be specified in a variety of ways, see Range variant. The class is extracted to be used in a
 * uniform way in many places.
 */
class BandExtractor {
public:
    /**
     * @brief Margin specified as a width of epsilon.
     */
    using WidthMargin = double;

    /**
     * @brief Margin specified as a number of energies around central value.
     */
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

    std::vector<std::size_t> getIndicesForCDFRange(const Eigensystem &eigensystem,
                                                   const arma::vec &normalizedEnergies,
                                                   const CDFRange &cdfRange, Logger &logger) const;

public:
    /**
     * @brief Constructs the extractor using a given band specification
     * @param range band specification, see Range variant and many of its components.
     * @param taskName name, given with first capital letter to inform about calculated ranges in getBandIndices() via
     * @a logger
     */
    explicit BandExtractor(Range range, std::string taskName = "Task")
            : range{std::move(range)}, taskName{std::move(taskName)}
    { }

    /**
     * @brief Returns the (subsequent) indices of states in @a eigensystem in the band specified in a constructor.
     */
    std::vector<std::size_t> getBandIndices(const Eigensystem &eigensystem,  Logger &logger) const;
};


#endif //MBL_ED_BANDEXTRACTOR_H
