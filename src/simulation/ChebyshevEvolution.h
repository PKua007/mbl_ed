//
// Created by pkua on 09.05.2020.
//

#ifndef MBL_ED_CHEBYSHEVEVOLUTION_H
#define MBL_ED_CHEBYSHEVEVOLUTION_H

#include <utility>

#include "evolution/TimeEvolution.h"
#include "simulation/SimulationsSpan.h"
#include "evolution/ChebyshevEvolver.h"
#include "core/HamiltonianGenerator.h"
#include "core/AveragingModel.h"
#include "core/RND.h"
#include "core/QuenchCalculator.h"
#include "simulation/RestorableSimulation.h"

/**
 * @brief A class performing time evolutions using Chebyshev expansion technique.
 * @details See CorrelationsTimeEvolution and its "slave" classes to see what is calculated. The hamiltonian generator
 * is prepared for each simulation according to a specific averaging model and, if desired, the quench stated is
 * prepared for evolution (see constructor). The template parameters default to standard classes and exist solely for
 * mocking purposes. See the default classes description for details of what they do.
 */
template<typename HamiltonianGenerator_t = HamiltonianGenerator, typename AveragingModel_t = AveragingModel,
         typename CorrelationsTimeEvolution_t = TimeEvolution,
         typename QuenchCalculator_t = QuenchCalculator, typename ChebyshevEvolver_t = ChebyshevEvolver>
class ChebyshevEvolution : public RestorableSimulation {
private:
    std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator;
    std::unique_ptr<AveragingModel_t> averagingModel;
    std::unique_ptr<RND> rnd;

private:
    std::unique_ptr<CorrelationsTimeEvolution_t> correlationsTimeEvolution;

    std::unique_ptr<QuenchCalculator_t> quenchCalculator;
    std::unique_ptr<HamiltonianGenerator_t> quenchHamiltonianGenerator;
    std::unique_ptr<RND> quenchRnd;


    auto prepareHamiltonianAndPossiblyQuenchVector(std::size_t simulationIndex, std::size_t totalSimulations,
                                                   Logger &logger) const
    {
        std::vector<arma::cx_vec> additionalVectors;

        this->averagingModel->setupHamiltonianGenerator(*this->hamiltonianGenerator, *this->rnd, simulationIndex,
                                                        totalSimulations);
        arma::sp_mat hamiltonian = this->hamiltonianGenerator->generate();

        if (this->quenchCalculator != nullptr) {
            this->averagingModel->setupHamiltonianGenerator(*this->quenchHamiltonianGenerator, *this->quenchRnd,
                                                            simulationIndex, totalSimulations);
            arma::sp_mat initialHamiltonian = this->quenchHamiltonianGenerator->generate();

            this->quenchCalculator->addQuench(initialHamiltonian, hamiltonian);

            // Yup, Armadillo wouldn't let you initialize arma::cx_vec by arma::vec easily, we have to do this nonsense
            const arma::vec &quenchedStateDouble = this->quenchCalculator->getLastQuenchedState();
            arma::cx_vec quenchedState(quenchedStateDouble.size());
            std::copy(quenchedStateDouble.begin(), quenchedStateDouble.end(), quenchedState.begin());
            additionalVectors.push_back(quenchedState);

            logger.info() << "Quenched state epsilon: " << this->quenchCalculator->getLastQuenchEpsilon();
            logger << "; quantum error: " << this->quenchCalculator->getLastQuenchEpsilonQuantumUncertainty() << ". ";
            logger << std::endl;
        }
        return std::make_pair(hamiltonian, additionalVectors);
    }

public:
    /**
     * @brief Constructor, where HamiltonianGenerator and RND for quench should be passed, or set to nullptr if
     * quench should not be done.
     * @details If quench is to be done, @a parameters.initialVectors has to have exactly one external vector slot
     * provided, 0 otherwise.
     */
    ChebyshevEvolution(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator,
                       std::unique_ptr<AveragingModel_t> averagingModel, std::unique_ptr<RND> rnd,
                       std::unique_ptr<CorrelationsTimeEvolution_t> correlationsTimeEvolution,
                       std::unique_ptr<QuenchCalculator_t> quenchCalculator,
                       std::unique_ptr<HamiltonianGenerator_t> quenchHamiltonianGenerator,
                       std::unique_ptr<RND> quenchRnd)
            : hamiltonianGenerator{std::move(hamiltonianGenerator)}, averagingModel{std::move(averagingModel)},
              rnd{std::move(rnd)}, correlationsTimeEvolution{std::move(correlationsTimeEvolution)},
              quenchCalculator{std::move(quenchCalculator)},
              quenchHamiltonianGenerator{std::move(quenchHamiltonianGenerator)}, quenchRnd{std::move(quenchRnd)}
    {
        if (this->quenchCalculator == nullptr) {
            Expects(this->correlationsTimeEvolution->countExternalVectors() == 0);
        } else {
            Expects(this->quenchHamiltonianGenerator != nullptr);
            Expects(this->quenchRnd != nullptr);
            Expects(this->correlationsTimeEvolution->countExternalVectors() == 1);
        }
    }

    /**
     * @brief Constructor with no quench.
     */
    ChebyshevEvolution(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator,
                       std::unique_ptr<AveragingModel_t> averagingModel, std::unique_ptr<RND> rnd,
                       std::unique_ptr<CorrelationsTimeEvolution_t> correlationsTimeEvolution)
            : ChebyshevEvolution(std::move(hamiltonianGenerator), std::move(averagingModel), std::move(rnd),
                                 std::move(correlationsTimeEvolution), nullptr, nullptr, nullptr)
    { }

    void printQuenchInfo(Logger &logger) {
        if (this->quenchCalculator != nullptr) {
            logger.info() << "Mean quench data: epsilon: " << this->quenchCalculator->getMeanEpsilon();
            logger << "; avg error: " << this->quenchCalculator->getEpsilonAveragingSampleError();
            logger << "; quantum error: " << this->quenchCalculator->getLastQuenchEpsilonQuantumUncertainty();
            logger << std::endl;
        }
    }

    void storeResults(std::ostream &out) const {
        this->correlationsTimeEvolution->storeResult(out);
    }

    void storeState(std::ostream &binaryOut) const override {
        bool doesDoQuench = (this->quenchCalculator != nullptr);
        binaryOut.write(reinterpret_cast<const char*>(&doesDoQuench), sizeof(doesDoQuench));
        Assert(binaryOut.good());
        if (doesDoQuench)
            this->quenchCalculator->storeState(binaryOut);
        this->correlationsTimeEvolution->storeState(binaryOut);
    }

    void joinRestoredState(std::istream &binaryIn) override {
        bool doesDoQuench = (this->quenchCalculator != nullptr);
        bool doesDoQuenchRestored{};
        binaryIn.read(reinterpret_cast<char*>(&doesDoQuenchRestored), sizeof(doesDoQuenchRestored));
        Assert(binaryIn.good());
        Assert(doesDoQuench == doesDoQuenchRestored);
        if (doesDoQuench)
            this->quenchCalculator->joinRestoredState(binaryIn);
        this->correlationsTimeEvolution->joinRestoredState(binaryIn);
    }

    void clear() override {
        this->correlationsTimeEvolution->clear();
        if (this->quenchCalculator != nullptr)
            this->quenchCalculator->clear();
    }

    void seedRandomGenerators(unsigned long seed) override {
        this->rnd->seed(seed);
        if (this->quenchRnd != nullptr)
            this->quenchRnd->seed(seed);
    }

    void performSimulation(std::size_t simulationIndex, std::size_t totalSimulations, Logger &logger) override {
        static_cast<void>(totalSimulations);

        arma::wall_clock wholeTimer;
        arma::wall_clock timer;

        wholeTimer.tic();
        logger.info() << "Performing evolution " << simulationIndex << "... " << std::endl;
        logger.verbose() << "Preparing hamiltonian started... " << std::endl;
        timer.tic();
        auto [hamiltonian, additionalVectors]
            = this->prepareHamiltonianAndPossiblyQuenchVector(simulationIndex, totalSimulations, logger);
        logger.info() << "Preparing hamiltonian done (" << timer.toc() << " s)" << std::endl;

        logger.verbose() << "Preparing evolver started... " << std::endl;
        timer.tic();
        ChebyshevEvolver_t evolver(hamiltonian, logger);
        logger.info() << "Preparing evolver done (" << timer.toc() << " s)." << std::endl;
        this->correlationsTimeEvolution->addEvolution(evolver, logger, additionalVectors);
        logger.info() << "Whole evolution took " << wholeTimer.toc() << " s." << std::endl;
    }

    [[nodiscard]] std::string getTagName() const override {
        return "chebyshev";
    }
};


#endif //MBL_ED_CHEBYSHEVEVOLUTION_H
