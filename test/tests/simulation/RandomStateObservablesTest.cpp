//
// Created by Piotr Kubala on 02/11/2020.
//

#include <catch2/catch.hpp>


#include "simulation/RandomStateObservables.h"

#include "core/observables/OnsiteOccupations.h"
#include "core/FockBasisGenerator.h"

TEST_CASE("RandomStateObservables: validation test", "[validation]") {
    auto basis = std::shared_ptr<FockBasis>(FockBasisGenerator{}.generate(6, 6));
    auto onsiteOccupations = std::make_shared<OnsiteOccupations>(basis);
    RandomStateObservables randomStateObservables(std::make_unique<RND>(), basis, {onsiteOccupations}, {},
                                                  {onsiteOccupations});
    std::ostringstream loggerStream;
    Logger logger(loggerStream);

    const std::size_t totalSimulations = 100;
    for (std::size_t i{}; i < totalSimulations; i++)
        randomStateObservables.performSimulation(i, totalSimulations, logger);

    auto header = randomStateObservables.getHeader();
    auto values = randomStateObservables.getValues();
    REQUIRE(header.size() == values.size());

    std::map<std::string, double> fields;
    for (std::size_t i{}; i < header.size(); i++)
        fields[header[i]] = std::stod(values[i]);

    // We are after 3 sigma tolerance. Errors are also checked to lie in resonable bounds
    INFO("n_2 = " << fields["n_2"] << " +/- " << fields["dn_2"]);
    REQUIRE((fields["dn_2"] > 0.005 && fields["dn_2"] < 0.02));
    CHECK(fields["n_2"] == Approx(1).margin(3 * fields["dn_2"]));

    INFO("n_4_norm = " << fields["n_4_norm"] << " +/- " << fields["dn_4_norm"]);
    REQUIRE((fields["dn_4_norm"] > 0.005 && fields["dn_4_norm"] < 0.02));
    CHECK(fields["n_4_norm"] == Approx(1).margin(3 * fields["dn_4_norm"]));
}

TEST_CASE("RandomStateObservables: header names") {
    auto basis = std::shared_ptr<FockBasis>(FockBasisGenerator{}.generate(2, 2));
    auto onsiteOccupations = std::make_shared<OnsiteOccupations>(basis);
    RandomStateObservables randomStateObservables(std::make_unique<RND>(), basis, {onsiteOccupations}, {},
                                                  {onsiteOccupations});

    std::vector<std::string> expected = {"n_1", "dn_1", "n_2", "dn_2", "n_1_norm", "dn_1_norm", "n_2_norm", "dn_2_norm"};
    REQUIRE(randomStateObservables.getHeader() == expected);
}