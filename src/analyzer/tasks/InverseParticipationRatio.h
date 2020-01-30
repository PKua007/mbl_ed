//
// Created by Piotr Kubala on 29/01/2020.
//

#ifndef MBL_ED_INVERSEPARTICIPATIONRATIO_H
#define MBL_ED_INVERSEPARTICIPATIONRATIO_H

#include "analyzer/BulkAnalyzerTask.h"

class InverseParticipationRatio : public BulkAnalyzerTask {
private:
    struct Entry {
        Entry() = default;
        Entry(double energy, double ipr) : energy(energy), ipr(ipr) { }

        double energy{};
        double ipr{};

        friend std::ostream &operator<<(std::ostream &out, const Entry &entry);
    };

    friend std::ostream &operator<<(std::ostream &out, const Entry &entry);

    double relativeMiddleEnergy{};
    double relativeMargin{};
    std::vector<double> ratios{};
    std::vector<Entry> entries{};

public:
    InverseParticipationRatio(double relativeMiddleEnergy, double relativeMargin);

    void analyze(const Eigensystem &eigensystem) override;
    [[nodiscard]] std::string getName() const override;
    void storeResult(std::ostream &out) const override;
};


#endif //MBL_ED_INVERSEPARTICIPATIONRATIO_H
