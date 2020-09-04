//
// Created by pkua on 08.11.2019.
//

#ifndef MBL_ED_EXACTDIAGONALIZATION_H
#define MBL_ED_EXACTDIAGONALIZATION_H


#include <iosfwd>
#include <random>
#include <memory>
#include <fstream>
#include <iterator>

#include <armadillo>
#include <utility>

#include "utils/Quantity.h"
#include "analyzer/Analyzer.h"
#include "core/Eigensystem.h"
#include "core/HamiltonianGenerator.h"
#include "ExactDiagonalizationParameters.h"
#include "core/RND.h"
#include "core/AveragingModel.h"
#include "simulation/RestorableSimulation.h"

/**
 * @brief A class performing a series of diagonalizations and optionaly some analyzer tasks.
 * @details Each simulation consists of: generating new hamiltonian according to @a AveragingModel_t,
 * diagonalisation, optionally saving eigenvalues and optionally doing some AnalyzerTask -s. Also, note, that the
 * template parameters default to standard classes and exist solely for testing purposes. See the default classes'
 * descriptions for the details of what they do.
 * @tparam HamiltonianGenerator_t the concrete hamiltonian generator to use
 * @tparam Analyzer_t the concrete analyzer to use
 * @tparam AveragingModel_t the concrete analyzer model
 */
template<typename HamiltonianGenerator_t = HamiltonianGenerator, typename AveragingModel_t = AveragingModel,
         typename Analyzer_t = Analyzer>
class ExactDiagonalization : public RestorableSimulation {
private:
    std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator;
    std::unique_ptr<AveragingModel_t> averagingModel;
    std::unique_ptr<RND> rnd;
    std::unique_ptr<FileOstreamProvider> ostreamProvider;
    std::unique_ptr<Analyzer_t> analyzer;
    ExactDiagonalizationParameters params;

    void doSaveEigenenergies(const Eigensystem &eigensystem, std::size_t index) const {
        std::ostringstream filenameStream;
        filenameStream << this->params.fileSignature << "_" << index << "_nrg.bin";
        std::string filename = filenameStream.str();
        auto out = this->ostreamProvider->openOutputFile(filename);
        eigensystem.store(*out);
    }

public:
    /**
     * @brief The constructor with mockable eigenenergy file creating using own FileOstreamProvider.
     */
    ExactDiagonalization(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator,
                         std::unique_ptr<AveragingModel_t> averagingModel, std::unique_ptr<RND> rnd,
                         std::unique_ptr<FileOstreamProvider> ostreamProvider,
                         ExactDiagonalizationParameters simulationParameters,
                         std::unique_ptr<Analyzer_t> analyzer)
            : hamiltonianGenerator{std::move(hamiltonianGenerator)}, averagingModel{std::move(averagingModel)},
              rnd{std::move(rnd)}, ostreamProvider{std::move(ostreamProvider)}, analyzer{std::move(analyzer)},
              params{std::move(simulationParameters)}
    { }

    /**
     * @brief The non-mockable constructor which will pass FileOstreamProvider which actually creates files to store
     * eigenenrgies.
     */
    ExactDiagonalization(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator,
                         std::unique_ptr<AveragingModel_t> averagingModel, std::unique_ptr<RND> rnd,
                         const ExactDiagonalizationParameters &simulationParameters,
                         std::unique_ptr<Analyzer_t> analyzer)
            : ExactDiagonalization(std::move(hamiltonianGenerator), std::move(averagingModel), std::move(rnd),
                                   std::make_unique<FileOstreamProvider>(), simulationParameters, std::move(analyzer))
    { }

    const Analyzer_t &getAnalyzer() const { return *this->analyzer; }

    void storeState(std::ostream &binaryOut) const override { this->analyzer->storeState(binaryOut); }
    void joinRestoredState(std::istream &binaryIn) override { this->analyzer->joinRestoredState(binaryIn); }
    void clear() override { this->analyzer->clear(); }

    void seedRandomGenerators(unsigned long seed) override {
        this->rnd->seed(seed);
    }

    /**
     * @brief Perform as many simulations as specified in SimulationParameters from teh constructor.
     * @details
     * <p> Before each simulation new hamiltonian is prepared according to @a AveragingModel_t. Then the
     * diagonalization is performed. After that, optionally, eigenenergies are stored and some on-the-fly
     * AnalyzerTask -s are performed.
     * <p> Eigenenergy files will be named `[SimulationParameters::fileSignature]_[simulation index]_ngr.txt`.
     * @param logger the output stream to log some info on the progress
     */
    void performSimulation(std::size_t simulationIndex, std::size_t totalSimulations, std::ostream &logger) override {
        logger << "[Simulation::perform] Performing diagonalization " << simulationIndex << "... " << std::flush;
        arma::wall_clock timer;
        timer.tic();
        this->averagingModel->setupHamiltonianGenerator(*this->hamiltonianGenerator, *this->rnd, simulationIndex,
                                                        totalSimulations);
        Eigensystem eigensystem = this->hamiltonianGenerator->calculateEigensystem(this->params.calculateEigenvectors);
        logger << "done (" << timer.toc() << " s). Analyzing... " << std::flush;

        timer.tic();
        this->analyzer->analyze(eigensystem, logger);
        if (this->params.saveEigenenergies)
            this->doSaveEigenenergies(eigensystem, simulationIndex);
        logger << "done (" << timer.toc() << " s)." << std::endl;
    }

    [[nodiscard]] std::string getTagName() const override {
        return "ed";
    }
};


#endif //MBL_ED_EXACTDIAGONALIZATION_H
