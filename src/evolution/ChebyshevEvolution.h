//
// Created by pkua on 09.05.2020.
//

#ifndef MBL_ED_CHEBYSHEVEVOLUTION_H
#define MBL_ED_CHEBYSHEVEVOLUTION_H

#include "CorrelationsTimeEvolution.h"
#include "CorrelationsTimeEvolutionParameters.h"
#include "ChebyshevEvolver.h"
#include "simulation/RND.h"

template<typename HamiltonianGenerator_t, typename AveragingModel_t>
class ChebyshevEvolution {
private:
    std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator;
    std::unique_ptr<RND> rnd;
    std::unique_ptr<FileOstreamProvider> ostreamProvider;
    CorrelationsTimeEvolution correlationsTimeEvolution;
    std::size_t from{};
    std::size_t to{};
    std::size_t totalSimulations{};
    std::string fileSignature{};

public:
    ChebyshevEvolution(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator, std::unique_ptr<RND> rnd,
                       std::unique_ptr<FileOstreamProvider> ostreamProvider, std::size_t from,
                       std::size_t to, std::size_t totalSimulations,
                       const CorrelationsTimeEvolutionParameters &parameters, const std::string &fileSignature)
            : hamiltonianGenerator{std::move(hamiltonianGenerator)}, rnd{std::move(rnd)},
              ostreamProvider{std::move(ostreamProvider)}, correlationsTimeEvolution(parameters),
              from{from}, to{to}, totalSimulations{totalSimulations}, fileSignature{fileSignature}
    {
        Expects(this->totalSimulations > 0);
        Expects(this->from < this->to);
        Expects(this->to <= this->totalSimulations);
    }

    ChebyshevEvolution(std::unique_ptr<HamiltonianGenerator_t> hamiltonianGenerator, std::unique_ptr<RND> rnd,
                       std::size_t from, std::size_t to, std::size_t totalSimulations,
                       const CorrelationsTimeEvolutionParameters &parameters, const std::string &fileSignature)
            : ChebyshevEvolution(std::move(hamiltonianGenerator), std::move(rnd), std::make_unique<FileOstreamProvider>(),
                                 from, to, totalSimulations, parameters, fileSignature)
    { }

    void perform(std::ostream &logger) {
        for (std::size_t i = this->from; i < this->to; i++) {
            arma::wall_clock wholeTimer;
            arma::wall_clock timer;

            wholeTimer.tic();
            logger << "[ChebyshevEvolution::perform] Performing evolution " << i << "... " << std::endl;
            logger << "[ChebyshevEvolution::perform] Preparing hamiltonian... " << std::flush;
            timer.tic();
            AveragingModel_t::setupHamiltonianGenerator(*hamiltonianGenerator, *rnd, i, this->totalSimulations);
            arma::sp_mat hamiltonian = this->hamiltonianGenerator->generate();
            logger << "done (" << timer.toc() << " s)" << std::endl;

            logger << "[ChebyshevEvolution::perform] Preparing evolver... " << std::endl;
            timer.tic();
            ChebyshevEvolver evolver(hamiltonian, logger);
            logger << "[ChebyshevEvolution::perform] Preparing evolver done (" << timer.toc() << " s)." << std::endl;
            this->correlationsTimeEvolution.addEvolution(evolver, logger);
            logger << "[ChebyshevEvolution::perform] Whole evolution took " << wholeTimer.toc() << " s." << std::endl;
        }

        std::string filename = this->fileSignature + "_evolution.txt";
        auto file = this->ostreamProvider->openFile(filename);
        this->correlationsTimeEvolution.storeResult(*file);
        logger << "[ChebyshevEvolution::perform] Observables stores to " << filename << std::endl;
    }
};


#endif //MBL_ED_CHEBYSHEVEVOLUTION_H
