//
// Created by Piotr Kubala on 09/02/2020.
//

#include <catch2/catch.hpp>

#include "core/terms/OnsiteDisorder.h"
#include "core/HamiltonianGenerator.h"
#include "core/FockBaseGenerator.h"
#include "core/RND.h"
#include "core/DisorderGenerator.h"

namespace {
    class SequenceDisorderGenerator : public DisorderGenerator {
    private:
        std::size_t index{};
        std::vector<double> sequence;

    public:
        explicit SequenceDisorderGenerator(std::vector<double> sequence) : sequence(std::move(sequence)) { }

        double generate(RND &rnd) override {
            static_cast<void>(rnd);
            return this->sequence.at(index++);
        }
    };
}

TEST_CASE("OnsiteDisorder: first sample") {
    HamiltonianGenerator generator(FockBaseGenerator{}.generate(3, 3), false);
    const auto &fockBase = *generator.getFockBase();
    RND rnd;
    auto disorderGenerator = std::make_unique<SequenceDisorderGenerator>(std::vector<double>{1, 2, 3});
    OnsiteDisorder onsiteDisorder(std::move(disorderGenerator), 3, rnd);

    REQUIRE(onsiteDisorder.calculate(fockBase[0], generator) == 3);
    REQUIRE(onsiteDisorder.calculate(fockBase[1], generator) == 4);
    REQUIRE(onsiteDisorder.calculate(fockBase[2], generator) == 5);
    REQUIRE(onsiteDisorder.calculate(fockBase[3], generator) == 5);
    REQUIRE(onsiteDisorder.calculate(fockBase[4], generator) == 6);
    REQUIRE(onsiteDisorder.calculate(fockBase[5], generator) == 7);
    REQUIRE(onsiteDisorder.calculate(fockBase[6], generator) == 6);
    REQUIRE(onsiteDisorder.calculate(fockBase[7], generator) == 7);
    REQUIRE(onsiteDisorder.calculate(fockBase[8], generator) == 8);
    REQUIRE(onsiteDisorder.calculate(fockBase[9], generator) == 9);
}

TEST_CASE("OnsiteDisorder: resample") {
    HamiltonianGenerator generator(FockBaseGenerator{}.generate(3, 3), false);
    const auto &fockBase = *generator.getFockBase();
    RND rnd;
    auto disorderGenerator = std::make_unique<SequenceDisorderGenerator>(std::vector<double>{0, 0, 0, 1, 2, 3});
    OnsiteDisorder onsiteDisorder(std::move(disorderGenerator), 3, rnd);

    onsiteDisorder.resampleOnsiteEnergies(rnd);

    REQUIRE(onsiteDisorder.calculate(fockBase[0], generator) == 3);
    REQUIRE(onsiteDisorder.calculate(fockBase[1], generator) == 4);
    REQUIRE(onsiteDisorder.calculate(fockBase[2], generator) == 5);
    REQUIRE(onsiteDisorder.calculate(fockBase[3], generator) == 5);
    REQUIRE(onsiteDisorder.calculate(fockBase[4], generator) == 6);
    REQUIRE(onsiteDisorder.calculate(fockBase[5], generator) == 7);
    REQUIRE(onsiteDisorder.calculate(fockBase[6], generator) == 6);
    REQUIRE(onsiteDisorder.calculate(fockBase[7], generator) == 7);
    REQUIRE(onsiteDisorder.calculate(fockBase[8], generator) == 8);
    REQUIRE(onsiteDisorder.calculate(fockBase[9], generator) == 9);
}