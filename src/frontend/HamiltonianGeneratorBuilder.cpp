//
// Created by pkua on 08.06.2020.
//

#include <ZipIterator.hpp>

#include "HamiltonianGeneratorBuilder.h"

#include "utils/Assertions.h"
#include "utils/Utils.h"
#include "CavityConstantsReader.h"

#include "core/terms/OnsiteDisorder.h"
#include "core/terms/HubbardHop.h"
#include "core/terms/CavityLongInteraction.h"
#include "core/terms/QuasiperiodicDisorder.h"
#include "core/terms/HubbardOnsite.h"
#include "core/terms/ListOnsite.h"
#include "core/terms/LookupCavityZ2.h"
#include "core/terms/LookupCavityYZ.h"
#include "core/terms/LookupCavityY2.h"

#include "core/disorder_generators/UniformGenerator.h"

namespace {
    std::unique_ptr<HubbardHop> parse_J_for_hubbard_hop(const Config &termParams) {
        std::vector<std::size_t> distances;
        std::vector<double> Js;

        for (const auto &key : termParams.getKeys()) {
            Validate(startsWith(key, "J"));
            std::size_t distance{};
            if (key == "J")
                distance = 1;
            else
                distance = std::stoul(key.substr(1));
            Validate(distance > 0);
            distances.push_back(distance);
            Js.push_back(termParams.getDouble(key));
        }

        Validate(std::any_of(Js.begin(), Js.end(), [](double J) { return J != 0; }));

        return std::make_unique<HubbardHop>(distances, Js);
    }
}

std::unique_ptr<HamiltonianGenerator>
HamiltonianGeneratorBuilder::build(const Parameters &params, std::shared_ptr<FockBasis> fockBasis, RND &rnd)
{
    std::size_t numberOfSites = fockBasis->getNumberOfSites();
    auto generator = std::make_unique<HamiltonianGenerator>(fockBasis, params.usePeriodicBC);

    for (auto &term : params.hamiltonianTerms) {
        std::string termName = term.first;
        const auto &termParams = term.second;
        if (termName == "hubbardHop") {
            generator->addHoppingTerm(parse_J_for_hubbard_hop(termParams));
        } else if (termName == "hubbardOnsite") {
            double U = termParams.getDouble("U");
            Validate(U >= 0);
            generator->addDiagonalTerm(std::make_unique<HubbardOnsite>(U));
        } else if (termName == "onsiteDisorder") {
            double W = termParams.getDouble("W");
            Validate(W >= 0);
            auto disorderGenerator = std::make_unique<UniformGenerator>(-W, W);
            generator->addDiagonalTerm(
                std::make_unique<OnsiteDisorder>(std::move(disorderGenerator), numberOfSites, rnd)
            );
        } else if (termName == "listOnsite") {
            auto stringValues = explode(termParams.getString("values"), ',');
            Validate(stringValues.size() == numberOfSites);
            std::vector<double> values(stringValues.size());
            std::transform(stringValues.begin(), stringValues.end(), values.begin(),
                           [](auto s) { return std::stod(s); });
            generator->addDiagonalTerm(std::make_unique<ListOnsite>(values));
        } else if (termName == "quasiperiodicDisorder") {
            double W = termParams.getDouble("W");
            double beta = termParams.getDouble("beta");
            double phi0 = termParams.getDouble("phi0");
            Validate(W >= 0);
            Validate(beta > 0);
            generator->addDiagonalTerm(std::make_unique<QuasiperiodicDisorder>(W, beta, phi0));
        } else if (termName == "cavityLongInteractions") {
            double U1 = termParams.getDouble("U1");
            double beta = termParams.getDouble("beta");
            double phi0 = termParams.getDouble("phi0");
            double phi0Bias{};
            try {
                phi0Bias = termParams.getDouble("phi0Bias");
            } catch (ConfigNoFieldException &) { }
            Validate(beta > 0);
            generator->addDiagonalTerm(std::make_unique<CavityLongInteraction>(U1, beta, phi0, phi0Bias));
        } else if (termName == "lookupCavityZ2") {
            double U1 = termParams.getDouble("U1");
            Validate(U1 >= 0);
            std::string cavityConstantsFilename = termParams.getString("ccfile");
            std::ifstream cavityConstantsFile(cavityConstantsFilename);
            if (!cavityConstantsFile)
                throw std::runtime_error("Cannot open " + cavityConstantsFilename + " to read cavity constants");
            CavityConstants cavityConstants = CavityConstantsReader::load(cavityConstantsFile);
            generator->addDiagonalTerm(std::make_unique<LookupCavityZ2>(U1, cavityConstants));
        } else if (termName == "lookupCavityYZ") {
            double U1 = termParams.getDouble("U1");
            Validate(U1 >= 0);
            std::string cavityConstantsFilename = termParams.getString("ccfile");
            std::ifstream cavityConstantsFile(cavityConstantsFilename);
            if (!cavityConstantsFile)
                throw std::runtime_error("Cannot open " + cavityConstantsFilename + " to read cavity constants");
            CavityConstants cavityConstants = CavityConstantsReader::load(cavityConstantsFile);
            generator->addHoppingTerm(std::make_unique<LookupCavityYZ>(U1, cavityConstants));
        } else if (termName == "lookupCavityY2") {
            double U1 = termParams.getDouble("U1");
            Validate(U1 >= 0);
            std::string cavityConstantsFilename = termParams.getString("ccfile");
            std::ifstream cavityConstantsFile(cavityConstantsFilename);
            if (!cavityConstantsFile)
                throw std::runtime_error("Cannot open " + cavityConstantsFilename + " to read cavity constants");
            CavityConstants cavityConstants = CavityConstantsReader::load(cavityConstantsFile);
            generator->addDoubleHoppingTerm(std::make_unique<LookupCavityY2>(U1, cavityConstants));
        } else if (termName == "lookupCavityZ2_YZ") {
            double U1 = termParams.getDouble("U1");
            Validate(U1 >= 0);
            std::string cavityConstantsFilename = termParams.getString("ccfile");
            std::ifstream cavityConstantsFile(cavityConstantsFilename);
            if (!cavityConstantsFile)
                throw std::runtime_error("Cannot open " + cavityConstantsFilename + " to read cavity constants");
            CavityConstants cavityConstants = CavityConstantsReader::load(cavityConstantsFile);
            generator->addDiagonalTerm(std::make_unique<LookupCavityZ2>(U1, cavityConstants));
            generator->addHoppingTerm(std::make_unique<LookupCavityYZ>(U1, cavityConstants));
        } else if (termName == "lookupCavityZ2_YZ_Y2") {
            double U1 = termParams.getDouble("U1");
            Validate(U1 >= 0);
            std::string cavityConstantsFilename = termParams.getString("ccfile");
            std::ifstream cavityConstantsFile(cavityConstantsFilename);
            if (!cavityConstantsFile)
                throw std::runtime_error("Cannot open " + cavityConstantsFilename + " to read cavity constants");
            CavityConstants cavityConstants = CavityConstantsReader::load(cavityConstantsFile);
            generator->addDiagonalTerm(std::make_unique<LookupCavityZ2>(U1, cavityConstants));
            generator->addHoppingTerm(std::make_unique<LookupCavityYZ>(U1, cavityConstants));
            generator->addDoubleHoppingTerm(std::make_unique<LookupCavityY2>(U1, cavityConstants));
        } else {
            throw ValidationException("Unknown hamiltonian term: " + termName);
        }
    }

    return generator;
}
