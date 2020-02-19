//
// Created by Piotr Kubala on 17/02/2020.
//

#ifndef MBL_ED_CORRELATIONSTIMEEVOLUTION_H
#define MBL_ED_CORRELATIONSTIMEEVOLUTION_H

#include "analyzer/BulkAnalyzerTask.h"
#include "utils/Assertions.h"

class CorrelationsTimeEvolution : public BulkAnalyzerTask {
public:
    enum TimeScaleType {
        Linear,
        Logarithmic
    };

private:
    class SymmetricMatrix {
    private:
        std::size_t size_{};
        std::vector<double> elements;

        [[nodiscard]] std::size_t idx(std::size_t i, std::size_t j) const {
            Expects(i < this->size_);
            Expects(j < this->size_);

            std::size_t smaller = std::min(i, j);
            std::size_t bigger = std::max(i, j);

            return bigger*(bigger + 1)/2 + smaller;
        }

    public:
        SymmetricMatrix() = default;
        explicit SymmetricMatrix(std::size_t size_) : size_{size_}, elements(size_*(size_ + 1)/2) { }

        [[nodiscard]] std::size_t size() const { return this->size_; }
        [[nodiscard]] double &operator()(std::size_t i, std::size_t j) { return this->elements[this->idx(i, j)]; }
        [[nodiscard]] double operator()(std::size_t i, std::size_t j) const { return this->elements[this->idx(i, j)]; }
    };

    struct Observables {
        double time{};
        std::vector<double> ns;
        SymmetricMatrix nns;

        Observables() = default;
        explicit Observables(std::size_t numberOfSites, double time)
                : time{time}, ns(numberOfSites), nns(numberOfSites)
        { }
    };

    struct Correlations {
        std::size_t distance{};
        double G{};

        void addObservables(const Observables &observables, std::size_t borderSize);
        [[nodiscard]] std::string getHeader() const;
        friend std::ostream &operator<<(std::ostream &out, const Correlations &corelations);
    };
    friend std::ostream &operator<<(std::ostream &out, const Correlations &corelations);

    struct OnsiteFluctuations {
        std::size_t i{};
        double rho{};

        void addObservables(const Observables &observables);
        [[nodiscard]] std::string getHeader() const;
        friend std::ostream &operator<<(std::ostream &out, const OnsiteFluctuations &onsiteFluctuations);
    };
    friend std::ostream &operator<<(std::ostream &out, const OnsiteFluctuations &onsiteFluctuations);

    struct TimeEntry {
        double t{};
        std::vector<Correlations> correlations{};
        std::vector<Correlations> borderlessCorrelations{};
        std::vector<OnsiteFluctuations> onsiteFluctuations{};
        double x{};

        [[nodiscard]] std::string getHeader() const;
        friend std::ostream &operator<<(std::ostream &out, const TimeEntry &timeEntry);
    };
    friend std::ostream &operator<<(std::ostream &out, const TimeEntry &timeEntry);

    struct VectorEvolution {
        FockBase::Vector initialVector;
        std::vector<TimeEntry> timeEntries{};

        [[nodiscard]] std::string getHeader() const;
    };

    std::size_t borderSize{};
    std::vector<VectorEvolution> evolutions{};
    std::vector<double> times{};

    [[nodiscard]] std::size_t getNumberOfSites() const;
    [[nodiscard]] SymmetricMatrix numberOfParticlesObservable(const FockBase &fockBase, const arma::mat &eigenvectors,
                                                              std::size_t site) const;
    [[nodiscard]] SymmetricMatrix numberOfParticlesSquaredObservable(const FockBase &fockBase,
                                                                     const arma::mat &eigenvectors,
                                                                     std::size_t site1, std::size_t site2) const;
    [[nodiscard]] SymmetricMatrix calculateEvolutionTerms(SymmetricMatrix matrixElements, const FockBase &fockBase,
                                                          const arma::mat &eigenvectors, std::size_t initialIdx) const;

    [[nodiscard]] double calculateObservableValue(const SymmetricMatrix &evolutionTerms, double time,
                                                  const arma::vec &eigenenergies) const;

public:
    CorrelationsTimeEvolution(double minTime, double maxTime, std::size_t numSteps, TimeScaleType timeScaleType,
                              std::size_t borderSize, const std::vector<FockBase::Vector> &vectorsToEvolve);

    void analyze(const Eigensystem &eigensystem) override;
    [[nodiscard]] std::string getName() const override;
    void storeResult(std::ostream &out) const override;

    std::vector<Observables>
    performObservblesEvolution(size_t initialIdx, const FockBase &fockBase, const Eigensystem &eigensystem) const;
};


#endif //MBL_ED_CORRELATIONSTIMEEVOLUTION_H
