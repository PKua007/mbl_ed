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
class ExactDiagonalization {
private:
    std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator;
    std::unique_ptr<AveragingModel_t> averagingModel;
    std::unique_ptr<RND> rnd;
    std::unique_ptr<FileOstreamProvider> ostreamProvider;
    ExactDiagonalizationParameters params;

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
    ExactDiagonalization(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator,
                         std::unique_ptr<AveragingModel_t> averagingModel, std::unique_ptr<RND> rnd,
                         std::unique_ptr<FileOstreamProvider> ostreamProvider, ExactDiagonalizationParameters simulationParameters)
            : hamiltonianGenerator{std::move(hamiltonianGenerator)}, averagingModel{std::move(averagingModel)},
              rnd{std::move(rnd)}, ostreamProvider{std::move(ostreamProvider)}, params{std::move(simulationParameters)}
    {
        Expects(this->params.simulationsSpan.total > 0);
        Expects(this->params.simulationsSpan.from < this->params.simulationsSpan.to);
        Expects(this->params.simulationsSpan.to <= this->params.simulationsSpan.total);
    }

    /**
     * @brief The non-mockable constructor which will pass FileOstreamProvider which actually creates files to store
     * eigenenrgies.
     */
    ExactDiagonalization(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator,
                         std::unique_ptr<AveragingModel_t> averagingModel, std::unique_ptr<RND> rnd,
                         const ExactDiagonalizationParameters &simulationParameters)
            : ExactDiagonalization(std::move(hamiltonianGenerator), std::move(averagingModel), std::move(rnd),
                                   std::make_unique<FileOstreamProvider>(), simulationParameters)
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
        for (std::size_t i = this->params.simulationsSpan.from; i < this->params.simulationsSpan.to; i++) {
            logger << "[Simulation::perform] Performing diagonalization " << i << "... " << std::flush;
            arma::wall_clock timer;
            timer.tic();
            this->averagingModel->setupHamiltonianGenerator(*this->hamiltonianGenerator, *this->rnd, i,
                                                            this->params.simulationsSpan.total);
            Eigensystem eigensystem = this->hamiltonianGenerator->calculateEigensystem(params.calculateEigenvectors);
            logger << "done (" << timer.toc() << " s). Analyzing... " << std::flush;

            timer.tic();
            analyzer.analyze(eigensystem, logger);
            if (this->params.saveEigenenergies)
                this->doSaveEigenenergies(eigensystem, i);
            logger << "done (" << timer.toc() << " s)." << std::endl;
        }
    }
};


#endif //MBL_ED_EXACTDIAGONALIZATION_H
