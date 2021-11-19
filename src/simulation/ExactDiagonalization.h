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
 * @brief A class performing diagonalizations and optionaly some analyzer tasks.
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

    void doSaveEigensystem(const Eigensystem &eigensystem, std::size_t index) const {
        std::ostringstream filenamePrefixStream;
        filenamePrefixStream << this->params.fileSignature << "_" << index;
        std::string filename = filenamePrefixStream.str();

        switch (this->params.storeLevel) {
            case ExactDiagonalizationParameters::StoreLevel::NONE:
                break;
            case ExactDiagonalizationParameters::StoreLevel::EIGENENERGIES: {
                auto out = this->ostreamProvider->openOutputFile(filename + "_nrg.bin");
                eigensystem.store(*out, this->params.fileType);
                break;
            }
            case ExactDiagonalizationParameters::StoreLevel::EIGENSYSTEM: {
                auto energyOut = this->ostreamProvider->openOutputFile(filename + "_nrg.bin");
                auto vectorOut = this->ostreamProvider->openOutputFile(filename + "_st.bin");
                eigensystem.store(*energyOut, *vectorOut, this->params.fileType);
                break;
            }
        }
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
     * @brief Perform simulation @a simulationIndex out of @a totalSimulations.
     * @details
     * <p> Before the simulation new hamiltonian is prepared according to @a AveragingModel_t. Then the
     * diagonalization is performed. After that, optionally, eigenenergies are stored and some on-the-fly
     * AnalyzerTask -s are performed.
     * <p> Eigenenergy files will be named `[Parameter::fileSignature]_[simulation index]_ngr.txt`.
     */
    void performSimulation(std::size_t simulationIndex, std::size_t totalSimulations, Logger &logger) override {
        logger.verbose() << "Performing diagonalization " << simulationIndex << " started..." << std::endl;
        arma::wall_clock timer;
        timer.tic();
        this->averagingModel->setupHamiltonianGenerator(*this->hamiltonianGenerator, *this->rnd, simulationIndex,
                                                        totalSimulations);
        Eigensystem eigensystem = this->hamiltonianGenerator->calculateEigensystem(this->params.calculateEigenvectors);
        double diagonalizationTime = timer.toc();

        logger.verbose() << "Performing analysis started..." << std::endl;
        timer.tic();
        this->analyzer->analyze(eigensystem, logger);
        this->doSaveEigensystem(eigensystem, simulationIndex);
        double analyzingTime = timer.toc();

        logger.info() << "Diagonalization " << simulationIndex << " done (diagonalization: " << diagonalizationTime;
        logger << "s, analysis: " << analyzingTime << " s)." << std::endl;
    }

    [[nodiscard]] std::string getTagName() const override {
        return "ed";
    }
};


#endif //MBL_ED_EXACTDIAGONALIZATION_H
