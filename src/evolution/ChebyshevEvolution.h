//
// Created by pkua on 09.05.2020.
//

#ifndef MBL_ED_CHEBYSHEVEVOLUTION_H
#define MBL_ED_CHEBYSHEVEVOLUTION_H

#include <utility>

#include "CorrelationsTimeEvolution.h"
#include "CorrelationsTimeEvolutionParameters.h"
#include "ChebyshevEvolver.h"
#include "simulation/RND.h"
#include "simulation/QuenchCalculator.h"

/**
 * @brief A class performing a whole series of time evolutions using Chebyshev expansion technique.
 * @details See CorrelationsTimeEvolution and its "slave" classes to see what is calculated. The hamiltonian generator
 * is prepared for each simulation according to a specific averaging model. The template parameters default to standard
 * classes and exist solely for mocking purposes. See the default classes description for details of what they do.
 */
template<typename HamiltonianGenerator_t = HamiltonianGenerator, typename AveragingModel_t = AveragingModel>
class ChebyshevEvolution {
private:
    std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator;
    std::unique_ptr<AveragingModel_t> averagingModel;
    std::unique_ptr<RND> rnd;
    std::unique_ptr<FileOstreamProvider> ostreamProvider;
    CorrelationsTimeEvolution correlationsTimeEvolution;
    CorrelationsTimeEvolutionParameters correlationsTimeEvolutionParameters;

    std::unique_ptr<HamiltonianGenerator_t> quenchHamiltonianGenerator;
    std::unique_ptr<RND> quenchRnd;

    std::size_t from{};
    std::size_t to{};
    std::size_t totalSimulations{};

    std::string fileSignature{};


    auto prepareHamiltonianAndAdditionalVectors(QuenchCalculator &quenchCalculator, std::size_t simulationIndex,
                                                std::ostream &logger) const
    {
        std::vector<arma::cx_vec> additionalVectors;

        this->averagingModel->setupHamiltonianGenerator(*this->hamiltonianGenerator, *this->rnd, simulationIndex,
                                                        this->totalSimulations);
        arma::sp_mat hamiltonian = this->hamiltonianGenerator->generate();

        if (this->quenchHamiltonianGenerator != nullptr) {
            this->averagingModel->setupHamiltonianGenerator(*this->quenchHamiltonianGenerator, *this->quenchRnd,
                                                            simulationIndex, this->totalSimulations);
            arma::sp_mat initialHamiltonian = this->quenchHamiltonianGenerator->generate();

            quenchCalculator.addQuench(initialHamiltonian, hamiltonian);

            // Yup, Armadillo wouldn't let you initialize arma::cx_vec by arma::vec easily, we have to do this nonsense
            const arma::vec &quenchedStateDouble = quenchCalculator.getLastQuenchedState();
            arma::cx_vec quenchedState(quenchedStateDouble.size());
            std::copy(quenchedStateDouble.begin(), quenchedStateDouble.end(), quenchedState.begin());
            additionalVectors.push_back(quenchedState);

            logger << "quenched state epsilon: " << quenchCalculator.getLastQuenchEpsilon() << "; quantum error: ";
            logger << quenchCalculator.getLastQuenchEpsilonQuantumUncertainty() << ". ";
        }
        return std::make_pair(hamiltonian, additionalVectors);
    }

public:
    /**
     * @brief Constructor with FileOstreamProvider injection for mocking.
     */
    ChebyshevEvolution(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator,
                       std::unique_ptr<AveragingModel_t> averagingModel, std::unique_ptr<RND> rnd,
                       std::unique_ptr<FileOstreamProvider> ostreamProvider, std::size_t from, std::size_t to,
                       std::size_t totalSimulations, const CorrelationsTimeEvolutionParameters &parameters,
                       std::string fileSignature)
            : hamiltonianGenerator{std::move(hamiltonianGenerator)}, averagingModel{std::move(averagingModel)},
              rnd{std::move(rnd)}, ostreamProvider{std::move(ostreamProvider)}, correlationsTimeEvolution(parameters),
              correlationsTimeEvolutionParameters{parameters}, from{from}, to{to}, totalSimulations{totalSimulations},
              fileSignature{std::move(fileSignature)}

    {
        Expects(this->totalSimulations > 0);
        Expects(this->from < this->to);
        Expects(this->to <= this->totalSimulations);
    }

    /**
     * @brief Constructor with default FileOstreamProvider behaviour. The rest of parameters are self-explanatory.
     */
    ChebyshevEvolution(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator,
                       std::unique_ptr<AveragingModel_t> averagingModel, std::unique_ptr<RND> rnd, std::size_t from,
                       std::size_t to, std::size_t totalSimulations,
                       const CorrelationsTimeEvolutionParameters &parameters, const std::string &fileSignature)
            : ChebyshevEvolution(std::move(hamiltonianGenerator), std::move(averagingModel), std::move(rnd),
                                 std::make_unique<FileOstreamProvider>(), from, to, totalSimulations, parameters,
                                 fileSignature)
    { }

    /**
     * @brief Installs hamiltonian generator and random number generator which will be used to perform quenches.
     * @details If it is used, one slot for CorrelationTimeEvolutionParameters::ExternalVector should be present
     * in @a parameters passed in the constructor.
     */
    void addQuenchHamiltonianGenerator(std::unique_ptr<HamiltonianGenerator_t> quenchHamiltonianGenerator,
                                       std::unique_ptr<RND> quenchRnd)
    {
        this->quenchHamiltonianGenerator = std::move(quenchHamiltonianGenerator);
        this->quenchRnd = std::move(quenchRnd);
    }

    /**
     * @brief Performs the simuations. The results (see CorrelationsTimeEvolution) are stored afterwards using
     * FileOstreamProvider from the constructor.
     * @details The name of the file is in the code, go check ;).
     */
    void perform(std::ostream &logger) {
        if (this->quenchHamiltonianGenerator == nullptr) {
            Assert(this->correlationsTimeEvolutionParameters.countExternalVectors() == 0);
        } else {
            Assert(this->correlationsTimeEvolutionParameters.countExternalVectors() == 1);
        }

        QuenchCalculator quenchCalculator;
        for (std::size_t i = this->from; i < this->to; i++) {
            arma::wall_clock wholeTimer;
            arma::wall_clock timer;

            wholeTimer.tic();
            logger << "[ChebyshevEvolution::perform] Performing evolution " << i << "... " << std::endl;
            logger << "[ChebyshevEvolution::perform] Preparing hamiltonian... " << std::flush;
            timer.tic();
            auto [hamiltonian, additionalVectors]
                = this->prepareHamiltonianAndAdditionalVectors(quenchCalculator, i, logger);
            logger << "done (" << timer.toc() << " s)" << std::endl;

            logger << "[ChebyshevEvolution::perform] Preparing evolver... " << std::endl;
            timer.tic();
            ChebyshevEvolver evolver(hamiltonian, logger);
            logger << "[ChebyshevEvolution::perform] Preparing evolver done (" << timer.toc() << " s)." << std::endl;
            this->correlationsTimeEvolution.addEvolution(evolver, logger, additionalVectors);
            logger << "[ChebyshevEvolution::perform] Whole evolution took " << wholeTimer.toc() << " s." << std::endl;
        }

        if (this->quenchHamiltonianGenerator != nullptr) {
            logger << "[ChebyshevEvolution::perform] Mean quench data: epsilon: " << quenchCalculator.getMeanEpsilon();
            logger << "; avg error: " << quenchCalculator.getEpsilonAveragingSampleError() << "; quantum error: ";
            logger << quenchCalculator.getLastQuenchEpsilonQuantumUncertainty() << std::endl;
        }

        std::string filename = this->fileSignature + "_evolution.txt";
        auto file = this->ostreamProvider->openFile(filename);
        this->correlationsTimeEvolution.storeResult(*file);
        logger << "[ChebyshevEvolution::perform] Observables stored to " << filename << std::endl;
    }
};


#endif //MBL_ED_CHEBYSHEVEVOLUTION_H
