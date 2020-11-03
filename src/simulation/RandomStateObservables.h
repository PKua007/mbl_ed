//
// Created by Piotr Kubala on 19/10/2020.
//

#ifndef MBL_ED_RANDOMSTATEOBSERVABLES_H
#define MBL_ED_RANDOMSTATEOBSERVABLES_H

#include <memory>
#include <vector>

#include "RestorableSimulation.h"

#include "core/FockBasis.h"
#include "core/RND.h"
#include "core/PrimaryObservable.h"
#include "core/SecondaryObservable.h"

class RandomStateObservables : public RestorableSimulation {
private:
    std::unique_ptr<RND> rnd;
    std::shared_ptr<FockBasis> basis;
    std::vector<std::shared_ptr<PrimaryObservable>> primaryObservables;
    std::vector<std::shared_ptr<SecondaryObservable>> secondaryObservables;
    std::vector<std::shared_ptr<Observable>> storedObservables;

    std::vector<std::string> header;
    std::vector<std::vector<double>> normalizedValues;
    std::vector<std::vector<double>> notNormalizedValues;

    std::size_t numValues{};

    double nextGaussian(double var);
    [[nodiscard]] std::size_t countStoredObservableValues() const;
    [[nodiscard]] std::vector<std::string> generateStoredObservablesHeader() const;
    void addStateToObservables(const arma::cx_vec &state,std::vector<std::vector<double>> &observables) const;
    void appendObservable(std::vector<std::string> &result, const std::vector<double> &observableValues) const;

public:
    RandomStateObservables(std::unique_ptr <RND> rnd, std::shared_ptr<FockBasis> basis,
                           std::vector <std::shared_ptr<PrimaryObservable>> primaryObservables,
                           std::vector <std::shared_ptr<SecondaryObservable>> secondaryObservables,
                           std::vector <std::shared_ptr<Observable>> storedObservables)
            : rnd{std::move(rnd)}, basis{std::move(basis)}, primaryObservables{std::move(primaryObservables)},
              secondaryObservables{std::move(secondaryObservables)}, storedObservables{std::move(storedObservables)}
    {
        this->numValues = this->countStoredObservableValues();
        this->normalizedValues.resize(this->numValues);
        this->notNormalizedValues.resize(this->numValues);
        this->header = this->generateStoredObservablesHeader();
    }

    void storeState(std::ostream &binaryOut) const override;
    void joinRestoredState(std::istream &binaryIn) override;
    void clear() override;

    void seedRandomGenerators(unsigned long seed) override { this->rnd->seed(seed); }
    void performSimulation(std::size_t simulationIndex, std::size_t totalSimulations, Logger &logger) override;
    [[nodiscard]] std::string getTagName() const override { return "rnd_vec"; }

    [[nodiscard]] std::vector<std::string> getHeader() const { return this->header; };
    [[nodiscard]] std::vector<std::string> getValues() const;
};


#endif //MBL_ED_RANDOMSTATEOBSERVABLES_H
