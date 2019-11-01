#include <iostream>
#include <iterator>
#include <algorithm>
#include <functional>
#include <armadillo>

void printCurrent(std::vector<int> vec) {
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<int>(std::cout, ", "));
    std::cout << std::endl;
}

int main()
{
    int N = 3;
    int M = 3;

    std::vector<int> current(M, 0);
    current[0] = N;
    printCurrent(current);

    while (current.back() != N) {
        int lastNonzeroK = M - 2;
        while (current[lastNonzeroK] == 0)
            lastNonzeroK--;

        current[lastNonzeroK]--;
        current[lastNonzeroK + 1] = N - std::accumulate(current.begin(), current.begin() + lastNonzeroK + 1, 0);
        std::fill(current.begin() + lastNonzeroK + 2, current.end(), 0);
        printCurrent(current);
    }

    return 0;
}