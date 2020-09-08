//
// Created by Piotr Kubala on 16/01/2020.
//

#include <functional>

#include "Analyzer.h"

#include "InlineAnalyzerTask.h"
#include "BulkAnalyzerTask.h"
#include "utils/Assertions.h"

void Analyzer::addTask(std::unique_ptr<AnalyzerTask> task) {
    this->tasks.push_back(std::move(task));
}

void Analyzer::analyze(const Eigensystem &eigensystem, Logger &logger) {
    for (auto &task : this->tasks)
        task->analyze(eigensystem, logger);
}

void Analyzer::storeBulkResults(const std::string &fileSignature) const {
    for (const auto &task : this->tasks) {
        try {
            auto &bulkTask = dynamic_cast<const BulkAnalyzerTask&>(*task);
            this->ostreamProvider->setFileDescription(task->getName());
            std::string name = fileSignature + "_" + task->getName() + ".txt";
            auto fileOut = this->ostreamProvider->openOutputFile(name);
            bulkTask.storeResult(*fileOut);
        } catch (std::bad_cast &e) { }
    }
}

std::vector<std::string> Analyzer::getInlineResultsHeader() const {
    std::vector<std::string> results;
    for (const auto &task : this->tasks) {
        try {
            auto &inlineTask = dynamic_cast<const InlineAnalyzerTask &>(*task);
            auto taskHeader = inlineTask.getResultHeader();
            results.insert(results.end(), taskHeader.begin(), taskHeader.end());
        } catch (std::bad_cast &e) { }
    }
    return results;
}

std::vector<std::string> Analyzer::getInlineResultsFields() const {
    std::vector<std::string> results;
    for (const auto &task : this->tasks) {
        try {
            auto &inlineTask = dynamic_cast<const InlineAnalyzerTask &>(*task);
            auto taskFields = inlineTask.getResultFields();
            results.insert(results.end(), taskFields.begin(), taskFields.end());
        } catch (std::bad_cast &e) { }
    }
    return results;
}

void Analyzer::storeState(std::ostream &binaryOut) const {
    std::size_t tasksSize = this->tasks.size();
    binaryOut.write(reinterpret_cast<const char*>(&tasksSize), sizeof(tasksSize));
    Assert(binaryOut.good());
    for (const auto &task : this->tasks)
        task->storeState(binaryOut);
}

void Analyzer::joinRestoredState(std::istream &binaryIn) {
    std::size_t tasksSizeRestored{};
    binaryIn.read(reinterpret_cast<char*>(&tasksSizeRestored), sizeof(tasksSizeRestored));
    Assert(binaryIn.good());
    Assert(this->tasks.size() == tasksSizeRestored);
    for (auto &task : this->tasks)
        task->joinRestoredState(binaryIn);
}

void Analyzer::clear() {
    for (auto &task : this->tasks)
        task->clear();
}
