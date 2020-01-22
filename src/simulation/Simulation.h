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
    std::size_t numberOfSimulations{};
    bool saveEigenenergies{};
    std::string fileSignature{};

    Eigensystem performSingleSimulation(std::size_t simulationIndex) {
        AveragingModel_t::setupHamiltonianGenerator(*hamiltonianGenerator, simulationIndex, this->numberOfSimulations);
        arma::mat hamiltonian = this->hamiltonianGenerator->generate();

        arma::vec armaEnergies;
        arma::mat armaEigvec;

        Assert(arma::eig_sym(armaEnergies, armaEigvec, hamiltonian));
        return Eigensystem(armaEnergies, armaEigvec);
    }

    void doSaveEigenenergies(const std::vector<double> &eigenenergies, std::size_t index) const {
        std::ostringstream filenameStream;
        filenameStream << this->fileSignature << "_" << index << "_nrg.dat";
        std::string filename = filenameStream.str();
        auto out = this->ostreamProvider->openFile(filename);

        std::copy(eigenenergies.begin(), eigenenergies.end(), std::ostream_iterator<double>(*out, "\n"));
    }

public:
    Simulation(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator,
               std::unique_ptr<FileOstreamProvider> ostreamProvider, size_t numberOfSimulations,
               std::string fileSignature, bool saveEigenenergies)
            : hamiltonianGenerator{std::move(hamiltonianGenerator)}, ostreamProvider{std::move(ostreamProvider)},
              numberOfSimulations{numberOfSimulations}, saveEigenenergies{saveEigenenergies},
              fileSignature{std::move(fileSignature)}
    { }

    Simulation(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator, size_t numberOfSimulations,
               const std::string &fileSignature, bool saveEigenenergies)
            : Simulation(std::move(hamiltonianGenerator), std::make_unique<FileOstreamProvider>(), numberOfSimulations,
                         fileSignature, saveEigenenergies)
    { }

    void perform(std::ostream &logger, Analyzer_t &analyzer) {
        for (std::size_t i = 0; i < this->numberOfSimulations; i++) {
            logger << "[Simulation::perform] Performing simulation " << i << "... " << std::flush;
            Eigensystem eigensystem = this->performSingleSimulation(i);
            analyzer.analyze(eigensystem);
            if (this->saveEigenenergies)
                this->doSaveEigenenergies(eigensystem.getEigenenergies(), i);
            logger << "done." << std::endl;
        }
    }
};


#endif //MBL_ED_SIMULATION_H
