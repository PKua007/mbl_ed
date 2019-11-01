#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <cstdlib>

#include "FockBaseGenerator.h"

void printCurrent(std::vector<int> vec) {
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<int>(std::cout, ", "));
    std::cout << std::endl;
}

int main() {
    FockBaseGenerator generator;
    auto base = generator.generate(3, 3);
    
    for (const auto &baseVector : base)
        printCurrent(baseVector);
    
    return EXIT_SUCCESS;
}