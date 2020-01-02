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

#include "Quantity.h"
#include "GapRatioCalculator.h"

namespace {
    struct FileOstreamProvider {
    public:
        static std::ofstream openFile(const std::string &filename) {
            std::ofstream out(filename);
            if (!out)
                throw std::runtime_error("Cannot open " + filename + " to store eigenvalues");
            return out;
        }
    };
}

template<typename HamiltonianGenerator_t, typename AveragingModel_t,
         typename GapRatioCalculator_t = GapRatioCalculator,
         typename FstreamProvider_t = FileOstreamProvider>
class Simulation {
private:
    std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator;
    std::size_t numberOfSimulations{};
    double relativeMiddleEnergy{};
    double relativeMargin{};
    bool saveEigenenergies{};

    Quantity meanGapRatio{};

    std::vector<double> performSingleSimulation(std::size_t i) {
        AveragingModel_t::setupHamiltonianGenerator(*hamiltonianGenerator, i, this->numberOfSimulations);
        arma::mat hamiltonian = this->hamiltonianGenerator->generate();
        arma::vec armaEigenenergies = arma::eig_sym(hamiltonian);

        std::vector<double> eigenenergies;
        eigenenergies.reserve(armaEigenenergies.size());
        std::copy(armaEigenenergies.begin(), armaEigenenergies.end(), std::back_inserter(eigenenergies));

        return eigenenergies;
    }

    void doSaveEigenenergies(const std::vector<double> &eigenenergies, std::size_t index) const {
        std::ostringstream filenameStream;
        filenameStream << this->hamiltonianGenerator->fileSignature() << "_" << index << "_nrg.dat";
        std::string filename = filenameStream.str();
        auto out = FstreamProvider_t::openFile(filename);

        std::copy(eigenenergies.begin(), eigenenergies.end(), std::ostream_iterator<double>(out, "\n"));
    }

public:
    Simulation(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator, size_t numberOfSimulations,
               double relativeMiddleEnergy, double relativeMargin, bool saveEigenenergies)
            : hamiltonianGenerator(std::move(hamiltonianGenerator)), numberOfSimulations{numberOfSimulations},
              relativeMiddleEnergy{relativeMiddleEnergy}, relativeMargin{relativeMargin},
              saveEigenenergies{saveEigenenergies}
    { }

    [[nodiscard]] Quantity getMeanGapRatio() const {
        return meanGapRatio;
    }

    void perform(std::ostream &logger) {
        GapRatioCalculator_t calculator(this->relativeMiddleEnergy, this->relativeMargin);

        for (std::size_t i = 0; i < this->numberOfSimulations; i++) {
            logger << "[Simulation::perform] Performing simulation " << i << "... " << std::flush;
            auto eigenenergies = this->performSingleSimulation(i);
            calculator.addEigenenergies(eigenenergies);
            if (this->saveEigenenergies)
                this->doSaveEigenenergies(eigenenergies, i);
            logger << "done." << std::endl;
        }

        this->meanGapRatio = calculator.calculateMean();
    }
};


#endif //MBL_ED_SIMULATION_H
