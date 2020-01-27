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

template<typename HamiltonianGenerator_t, typename AveragingModel_t, typename Analyzer_t = Analyzer>
class Simulation {
private:
    std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator;
    std::unique_ptr<FileOstreamProvider> ostreamProvider;
    SimulationParameters params;

    Eigensystem performSingleSimulation(std::size_t simulationIndex) {
        AveragingModel_t::setupHamiltonianGenerator(*hamiltonianGenerator, simulationIndex,
                                                    this->params.totalSimulations);
        arma::mat hamiltonian = this->hamiltonianGenerator->generate();

        arma::vec armaEnergies;
        arma::mat armaEigvec;

        if (this->params.calculateEigenvectors) {
            Assert(arma::eig_sym(armaEnergies, armaEigvec, hamiltonian));
            return Eigensystem(armaEnergies, armaEigvec);
        } else {
            Assert(arma::eig_sym(armaEnergies, hamiltonian));
            return Eigensystem(armaEnergies);
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
    Simulation(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator,
               std::unique_ptr<FileOstreamProvider> ostreamProvider, SimulationParameters simulationParameters)
            : hamiltonianGenerator{std::move(hamiltonianGenerator)}, ostreamProvider{std::move(ostreamProvider)},
              params{std::move(simulationParameters)}
    {
        Expects(this->params.totalSimulations > 0);
        Expects(this->params.from < this->params.to);
        Expects(this->params.to <= this->params.totalSimulations);
    }

    Simulation(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator,
               const SimulationParameters &simulationParameters)
            : Simulation(std::move(hamiltonianGenerator), std::make_unique<FileOstreamProvider>(), simulationParameters)
    { }

    void perform(std::ostream &logger, Analyzer_t &analyzer) {
        for (std::size_t i = this->params.from; i < this->params.to; i++) {
            logger << "[Simulation::perform] Performing simulation " << i << "... " << std::flush;
            Eigensystem eigensystem = this->performSingleSimulation(i);
            analyzer.analyze(eigensystem);
            if (this->params.saveEigenenergies)
                this->doSaveEigenenergies(eigensystem, i);
            logger << "done." << std::endl;
        }
    }
};


#endif //MBL_ED_SIMULATION_H
