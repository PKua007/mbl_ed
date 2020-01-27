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

template<typename HamiltonianGenerator_t, typename AveragingModel_t, typename Analyzer_t = Analyzer>
class Simulation {
private:
    std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator;
    std::unique_ptr<FileOstreamProvider> ostreamProvider;
    std::size_t from{};
    std::size_t to{};
    std::size_t totalSimulations{};
    bool calculateEigenvectors{};
    bool saveEigenenergies{};
    std::string fileSignature{};

    Eigensystem performSingleSimulation(std::size_t simulationIndex) {
        AveragingModel_t::setupHamiltonianGenerator(*hamiltonianGenerator, simulationIndex, this->totalSimulations);
        arma::mat hamiltonian = this->hamiltonianGenerator->generate();

        arma::vec armaEnergies;
        arma::mat armaEigvec;

        if (this->calculateEigenvectors) {
            Assert(arma::eig_sym(armaEnergies, armaEigvec, hamiltonian));
            return Eigensystem(armaEnergies, armaEigvec);
        } else {
            Assert(arma::eig_sym(armaEnergies, hamiltonian));
            return Eigensystem(armaEnergies);
        }
    }

    void doSaveEigenenergies(const Eigensystem &eigensystem, std::size_t index) const {
        std::ostringstream filenameStream;
        filenameStream << this->fileSignature << "_" << index << "_nrg.bin";
        std::string filename = filenameStream.str();
        auto out = this->ostreamProvider->openFile(filename);
        eigensystem.store(*out);
    }

public:
    Simulation(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator,
               std::unique_ptr<FileOstreamProvider> ostreamProvider, std::size_t from, std::size_t to,
               std::size_t totalSimulations, std::string fileSignature, bool calculateEigenvectors,
               bool saveEigenenergies)
            : hamiltonianGenerator{std::move(hamiltonianGenerator)}, ostreamProvider{std::move(ostreamProvider)},
              from{from}, to{to}, totalSimulations{totalSimulations}, calculateEigenvectors{calculateEigenvectors},
              saveEigenenergies{saveEigenenergies}, fileSignature{std::move(fileSignature)}
    {
        Expects(totalSimulations > 0);
        Expects(from < to);
        Expects(to <= totalSimulations);
    }

    Simulation(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator, std::size_t from, std::size_t to,
               std::size_t totalSimulations, const std::string &fileSignature, bool calculateEigenvectors,
               bool saveEigenenergies)
            : Simulation(std::move(hamiltonianGenerator), std::make_unique<FileOstreamProvider>(), from, to,
                         totalSimulations, fileSignature, calculateEigenvectors, saveEigenenergies)
    { }

    void perform(std::ostream &logger, Analyzer_t &analyzer) {
        for (std::size_t i = this->from; i < this->to; i++) {
            logger << "[Simulation::perform] Performing simulation " << i << "... " << std::flush;
            Eigensystem eigensystem = this->performSingleSimulation(i);
            analyzer.analyze(eigensystem);
            if (this->saveEigenenergies)
                this->doSaveEigenenergies(eigensystem, i);
            logger << "done." << std::endl;
        }
    }
};


#endif //MBL_ED_SIMULATION_H
