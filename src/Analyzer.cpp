//
// Created by Piotr Kubala on 16/01/2020.
//

#include <functional>

#include "Analyzer.h"

void Analyzer::addTask(std::unique_ptr<AnalyzerTask> task) {
    this->tasks.push_back(std::move(task));
}

void Analyzer::analyze(const std::vector<double> &eigenenergies) {
    using namespace std::placeholders;
    std::for_each(this->tasks.begin(), this->tasks.end(), std::bind(&AnalyzerTask::analyze, _1, eigenenergies));
}

void Analyzer::printResults(std::ostream &out) const {
    for (auto &task : this->tasks) {
        out << task->getName() << " results:" << std::endl;
        task->printResult(out);
    }
}

void Analyzer::storeResults(const std::string &fileSignature) const {
    for (auto &task : this->tasks) {
        this->ostreamProvider->setFileDescription(task->getName());
        std::string name = fileSignature + "_" + task->getName() + ".txt";
        auto fileOut = this->ostreamProvider->openFile(name);
        task->printResult(*fileOut);
    }
}
