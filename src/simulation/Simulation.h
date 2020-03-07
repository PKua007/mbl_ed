//
// Created by pkua on 08.11.2019.
//

#ifndef MBL_ED_SIMULATION_H
#define MBL_ED_SIMULATION_H


#include <iosfwd>
#include <random>
#include <memory>
#include <fstream>
#include <iterator>

#include <armadillo>
#include <utility>

#include "utils/Quantity.h"
#include "analyzer/tasks/MeanGapRatio.h"
#include "analyzer/Analyzer.h"
#include "Eigensystem.h"
#include "SimulationParameters.h"
#include "RND.h"

/**
 * @brief A class performing a series of diagonalizations and optionaly some analyzer tasks.
 * @details Each simulation consists of: generating new hamiltonian according to @a AveragingModel_t,
 * diagonalisation, optionally saving eigenvalues and optionally doing some AnalyzerTask -s.
 * @tparam HamiltonianGenerator_t the concrete hamiltonian generator to use
 * @tparam AveragingModel_t the concrete averaging model to use
 * @tparam Analyzer_t the concrete analyzer to use
 */
template<typename HamiltonianGenerator_t, typename AveragingModel_t, typename Analyzer_t = Analyzer>
class Simulation {
private:
    std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator;
    std::unique_ptr<RND> rnd;
    std::unique_ptr<FileOstreamProvider> ostreamProvider;
    SimulationParameters params;

    Eigensystem performSingleSimulation(std::size_t simulationIndex) {
        AveragingModel_t::setupHamiltonianGenerator(*hamiltonianGenerator, *rnd, simulationIndex,
                                                    this->params.totalSimulations);
        arma::mat hamiltonian = this->hamiltonianGenerator->generate();

        arma::vec armaEnergies;
        arma::mat armaEigvec;

        if (this->params.calculateEigenvectors) {
            Assert(arma::eig_sym(armaEnergies, armaEigvec, hamiltonian));
            return Eigensystem(armaEnergies, armaEigvec, this->hamiltonianGenerator->getFockBase());
        } else {
            Assert(arma::eig_sym(armaEnergies, hamiltonian));
            return Eigensystem(armaEnergies, this->hamiltonianGenerator->getFockBase());
        }
    }

    void doSaveEigenenergies(const Eigensystem &eigensystem, std::size_t index) const {
        std::ostringstream filenameStream;
        filenameStream << this->params.fileSignature << "_" << index << "_nrg.bin";
        std::string filename = filenameStream.str();
        auto out = this->ostreamProvider->openFile(filename);
        eigensystem.store(*out);
    }

public:
    /**
     * @brief The constructor with mockable eigenenergy file creating using own FileOstreamProvider.
     */
    Simulation(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator, std::unique_ptr<RND> rnd,
               std::unique_ptr<FileOstreamProvider> ostreamProvider, SimulationParameters simulationParameters)
            : hamiltonianGenerator{std::move(hamiltonianGenerator)}, rnd{std::move(rnd)},
              ostreamProvider{std::move(ostreamProvider)}, params{std::move(simulationParameters)}
    {
        Expects(this->params.totalSimulations > 0);
        Expects(this->params.from < this->params.to);
        Expects(this->params.to <= this->params.totalSimulations);
    }

    /**
     * @brief The non-mockable constructor which will pass FileOstreamProvider which actually creates files to store
     * eigenenrgies.
     */
    Simulation(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator, std::unique_ptr<RND> rnd,
               const SimulationParameters &simulationParameters)
            : Simulation(std::move(hamiltonianGenerator), std::move(rnd), std::make_unique<FileOstreamProvider>(),
                         simulationParameters)
    { }

    /**
     * @brief Perform as many simulations as specified in SimulationParameters from teh constructor.
     * @details
     * <p> Before each simulation new hamiltonian is prepared according to @a AveragingModel_t. Then the
     * diagonalization is performed. After that, optionally, eigenenergies are stored and some on-the-fly
     * AnalyzerTask -s are performed.
     * <p> Eigenenergy files will be named `[SimulationParameters::fileSignature]_[simulation index]_ngr.txt`.
     * @param logger the output stream to log some info on the progress
     * @param analyzer @a Analyzer_t object to user for on-the-fly analyzer tasks
     */
    void perform(std::ostream &logger, Analyzer_t &analyzer) {
        for (std::size_t i = this->params.from; i < this->params.to; i++) {
            logger << "[Simulation::perform] Performing diagonalization " << i << "... " << std::flush;
            arma::wall_clock timer;
            timer.tic();
            Eigensystem eigensystem = this->performSingleSimulation(i);
            logger << "done (" << timer.toc() << " s). Analyzing... " << std::flush;

            timer.tic();
            analyzer.analyze(eigensystem, logger);
            if (this->params.saveEigenenergies)
                this->doSaveEigenenergies(eigensystem, i);
            logger << "done (" << timer.toc() << " s)." << std::endl;
        }
    }
};


#endif //MBL_ED_SIMULATION_H
