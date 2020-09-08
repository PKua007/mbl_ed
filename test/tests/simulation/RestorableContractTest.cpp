//
// Created by Piotr Kubala on 05/09/2020.
//

#include <catch2/catch.hpp>

#include "simulation/Restorable.h"

#include "object_mothers/HamiltonianGeneratorMother.h"
#include "core/averaging_models/UniformPhi0AveragingModel.h"
#include "simulation/RestorableHelper.h"

#include "analyzer/Analyzer.h"
#include "analyzer/tasks/BulkMeanGapRatio.h"
#include "analyzer/tasks/CDF.h"
#include "analyzer/tasks/DressedStatesFinder.h"
#include "analyzer/tasks/InverseParticipationRatio.h"
#include "analyzer/tasks/MeanGapRatio.h"
#include "analyzer/tasks/MeanInverseParticipationRatio.h"

namespace {
    template<typename T>
    struct RestorableAccessor;


    class WithHubbardQuasiperiodicEigensystems {
    private:
        UniformPhi0AveragingModel averagingModel{};
        RND rnd{};
        std::unique_ptr<HamiltonianGenerator> hamiltonianGenerator;

    protected:
        Eigensystem eigensystem1;
        Eigensystem eigensystem2;

        WithHubbardQuasiperiodicEigensystems(std::size_t N, std::size_t K, double J, double U, double W, double beta,
                                             double phi0)
        {
            this->hamiltonianGenerator = HamiltonianGeneratorMother(N, K).hubbardQuasiperiodic(J, U, W, beta, phi0);
            this->averagingModel.setupHamiltonianGenerator(*this->hamiltonianGenerator, this->rnd, 1, 5);
            this->eigensystem1 = this->hamiltonianGenerator->calculateEigensystem(true);
            this->averagingModel.setupHamiltonianGenerator(*this->hamiltonianGenerator, this->rnd, 2, 5);
            this->eigensystem2 = this->hamiltonianGenerator->calculateEigensystem(true);
        }
    };


    struct WithGetResultForBulkAnalyzerTask {
        std::string getResult(const BulkAnalyzerTask &task) {
            std::ostringstream out;
            task.storeResult(out);
            return out.str();
        }
    };


    struct WithLogger {
        std::ostringstream loggerOut;
        Logger logger;

        WithLogger() : logger(this->loggerOut) { }
    };


    template<>
    struct RestorableAccessor<BulkMeanGapRatio> : public WithHubbardQuasiperiodicEigensystems,
                                                  public WithGetResultForBulkAnalyzerTask, public WithLogger
    {
        RestorableAccessor() : WithHubbardQuasiperiodicEigensystems(4, 4, 1, 1, 10, 0.3, 0) {}

        BulkMeanGapRatio generateRestorable() { return BulkMeanGapRatio(4); }
        void addFirstEntry(BulkMeanGapRatio &bmgr) { bmgr.analyze(this->eigensystem1, this->logger); }
        void addSecondEntry(BulkMeanGapRatio &bmgr) { bmgr.analyze(this->eigensystem2, this->logger); }
    };


    template<>
    struct RestorableAccessor<CDF> : public WithHubbardQuasiperiodicEigensystems,
                                     public WithGetResultForBulkAnalyzerTask, public WithLogger
    {
        RestorableAccessor() : WithHubbardQuasiperiodicEigensystems(4, 4, 1, 1, 10, 0.3, 0) {}

        CDF generateRestorable() { return CDF(4); }
        void addFirstEntry(CDF &cdf) { cdf.analyze(this->eigensystem1, this->logger); }
        void addSecondEntry(CDF &cdf) { cdf.analyze(this->eigensystem2, this->logger); }
    };


    template<>
    struct RestorableAccessor<DressedStatesFinder> : public WithHubbardQuasiperiodicEigensystems,
                                                     public WithGetResultForBulkAnalyzerTask, public WithLogger
    {
        RestorableAccessor() : WithHubbardQuasiperiodicEigensystems(4, 4, 1, 1, 10, 0.3, 0) {}

        DressedStatesFinder generateRestorable() { return DressedStatesFinder(0.5, 0.3, 0.9); }
        void addFirstEntry(DressedStatesFinder &dsf) { dsf.analyze(this->eigensystem1, this->logger); }
        void addSecondEntry(DressedStatesFinder &dsf) { dsf.analyze(this->eigensystem2, this->logger); }
    };


    template<>
    struct RestorableAccessor<InverseParticipationRatio> : public WithHubbardQuasiperiodicEigensystems,
                                                           public WithGetResultForBulkAnalyzerTask, public WithLogger
    {
        RestorableAccessor() : WithHubbardQuasiperiodicEigensystems(4, 4, 1, 1, 10, 0.3, 0) {}

        InverseParticipationRatio generateRestorable() { return InverseParticipationRatio(0.5, 0.3); }
        void addFirstEntry(InverseParticipationRatio &ipr) { ipr.analyze(this->eigensystem1, this->logger); }
        void addSecondEntry(InverseParticipationRatio &ipr) { ipr.analyze(this->eigensystem2, this->logger); }
    };


    template<>
    struct RestorableAccessor<MeanGapRatio> : public WithHubbardQuasiperiodicEigensystems, public WithLogger {
        RestorableAccessor() : WithHubbardQuasiperiodicEigensystems(4, 4, 1, 1, 10, 0.3, 0) {}

        MeanGapRatio generateRestorable() { return MeanGapRatio(0.5, 0.3); }
        void addFirstEntry(MeanGapRatio &mgr) { mgr.analyze(this->eigensystem1, this->logger); }
        void addSecondEntry(MeanGapRatio &mgr) { mgr.analyze(this->eigensystem2, this->logger); }
        auto getResult(const MeanGapRatio &mgr) { return mgr.getResultFields(); }
    };


    template<>
    struct RestorableAccessor<MeanInverseParticipationRatio> : public WithHubbardQuasiperiodicEigensystems,
                                                               public WithLogger
    {
        RestorableAccessor() : WithHubbardQuasiperiodicEigensystems(4, 4, 1, 1, 10, 0.3, 0) {}

        MeanInverseParticipationRatio generateRestorable() { return MeanInverseParticipationRatio(0.5, 0.3); }
        void addFirstEntry(MeanInverseParticipationRatio &mipr) { mipr.analyze(this->eigensystem1, this->logger); }
        void addSecondEntry(MeanInverseParticipationRatio &mipr) { mipr.analyze(this->eigensystem2, this->logger); }
        auto getResult(const MeanInverseParticipationRatio &mipr) { return mipr.getResultFields(); }
    };


    class MeanEnergyCalculator : public InlineAnalyzerTask, public WithLogger {
    private:
        bool normalized{};
        std::vector<double> energies;

    public:
        explicit MeanEnergyCalculator(bool normalized) : normalized{normalized}
        { }

        void analyze(const Eigensystem &eigensystem, Logger &) override {
            arma::vec newEnergies;
            if (this->normalized)
                newEnergies = eigensystem.getNormalizedEigenenergies();
            else
                newEnergies = eigensystem.getEigenenergies();
            this->energies.insert(this->energies.end(), newEnergies.begin(), newEnergies.end());
        }

        std::string getName() const override { return "mean_nrg"; }
        std::vector<std::string> getResultHeader() const override { return {"mean_energy"}; }

        std::vector<std::string> getResultFields() const override {
            double mean = std::accumulate(this->energies.begin(), this->energies.end(), 0.0) / this->energies.size();
            return {std::to_string(mean)};
        }

        void storeState(std::ostream &binaryOut) const override {
            RestorableHelper::storeStateForVector(this->energies, binaryOut);
        }

        void joinRestoredState(std::istream &binaryIn) override {
            RestorableHelper::joinRestoredStateForVector(this->energies, binaryIn);
        }

        void clear() override {
            this->energies.clear();
        }
    };


    template<>
    struct RestorableAccessor<Analyzer> : public WithHubbardQuasiperiodicEigensystems, public WithLogger {
        RestorableAccessor() : WithHubbardQuasiperiodicEigensystems(4, 4, 1, 1, 10, 0.3, 0) {}

        Analyzer generateRestorable() {
            Analyzer analyzer;
            analyzer.addTask(std::make_unique<MeanEnergyCalculator>(true));
            analyzer.addTask(std::make_unique<MeanEnergyCalculator>(false));
            return analyzer;
        }

        void addFirstEntry(Analyzer &analyzer) { analyzer.analyze(this->eigensystem1, this->logger); }
        void addSecondEntry(Analyzer &analyzer) { analyzer.analyze(this->eigensystem2, this->logger); }
        auto getResult(const Analyzer &analyzer) { return analyzer.getInlineResultsFields(); }
    };
}

TEMPLATE_TEST_CASE("Restorable: contract test", "", BulkMeanGapRatio, CDF, DressedStatesFinder,
                   InverseParticipationRatio, MeanGapRatio, MeanInverseParticipationRatio, Analyzer)
{
    RestorableAccessor<TestType> restorableAccessor;
    TestType restorable = restorableAccessor.generateRestorable();

    SECTION("clearing") {
        restorableAccessor.addFirstEntry(restorable);
        auto result1 = restorableAccessor.getResult(restorable);

        restorableAccessor.addSecondEntry(restorable);
        restorable.clear();
        restorableAccessor.addFirstEntry(restorable);
        auto result2 = restorableAccessor.getResult(restorable);

        REQUIRE(result1 == result2);
    }

    SECTION("storing and joining restored") {
        restorableAccessor.addFirstEntry(restorable);
        restorableAccessor.addSecondEntry(restorable);
        auto normalResult = restorableAccessor.getResult(restorable);

        restorable.clear();
        restorableAccessor.addSecondEntry(restorable);
        std::stringstream simulation2;
        restorable.storeState(simulation2);

        restorable.clear();
        restorableAccessor.addFirstEntry(restorable);
        restorable.joinRestoredState(simulation2);
        auto restoredResult = restorableAccessor.getResult(restorable);

        REQUIRE(normalResult == restoredResult);
    }
}