//
// Created by pkua on 02.05.2020.
//

#include <complex>

#include "ChebyshevEvolver.h"
#include "utils/Assertions.h"

using namespace std::complex_literals;

ChebyshevEvolver::ChebyshevEvolver(const arma::sp_mat &hamiltonian, double Nfactor)
        : hamiltonian{hamiltonian}, Nfactor{Nfactor}
{
    Expects(Nfactor > 0);
}

arma::cx_vec ChebyshevEvolver::evolveVector(const arma::cx_vec &initialState) {
    arma::sp_mat hamiltonianRescaled = (this->hamiltonian - arma::speye(arma::size(this->hamiltonian)) * this->b) / this->a;

    //this->chebyshevVectors.resize(N + 1);
    arma::cx_mat prevprev = initialState;
    arma::cx_mat prev = hamiltonianRescaled * initialState;
    arma::cx_mat current;
    arma::cx_vec result(arma::size(initialState), arma::fill::zeros);

    result += initialState * std::cyl_bessel_j(0, this->a * this->dt);
    result += 2. * -1i * std::cyl_bessel_j(1, this->a * this->dt) * prev;

    for (std::size_t i = 2; i <= N; i++) {
        //std::cout << "Building order " << i << std::endl;
        current = 2 * hamiltonianRescaled * prev - prevprev;
        result += 2. * std::pow(-1i, i) * std::cyl_bessel_j(i, this->a * this->dt) * current;
        prevprev = prev;
        prev = current;
    }

    result *= std::exp(-1i * this->b * this->dt);

    return result;
}

void ChebyshevEvolver::prepareFor(const arma::cx_vec &initialState, double maxTime,
                                  std::size_t steps) {
    Expects(maxTime > 0);
    Expects(steps >= 2);
    Expects(initialState.size() == this->hamiltonian.n_cols);

    this->t = 0;
    this->dt = maxTime / static_cast<double>(steps - 1);
    this->currentState = initialState;
    this->step = 0;
    this->steps = steps;

    std::size_t numVals = std::min<std::size_t>(20, this->hamiltonian.n_rows - 1);

    arma::vec minEigval, maxEigval;
    std::cout << "Calculating Emin... " << std::flush;
    arma::wall_clock timer;
    timer.tic();
    Assert(arma::eigs_sym(minEigval, this->hamiltonian, numVals, "sa"));
    std::cout << "done (" << timer.toc() << " s). Emax... " << std::flush;
    timer.tic();
    Assert(arma::eigs_sym(maxEigval, this->hamiltonian, numVals, "la"));
    std::cout << "done (" << timer.toc() << " s). " << std::endl;
    double Emin = minEigval.front();
    double Emax = maxEigval.back();
//    double Emin = -31.80130421278989;
//    double Emax = 306.9052833537735;
    this->a = (Emax - Emin) / 2;
    this->b = (Emax + Emin) / 2;

//    this->N = static_cast<size_t>(this->Nfactor * 2 * this->a * this->dt);
//    Assert(this->N >= 1);

    this->N = 1;
    arma::cx_vec evolved;
    double normDiff{};

    std::cout << "Optimizing Chebyshev expansion order:" << std::endl;
    do {
        this->N *= 2;
        std::cout << "Trying " << this->N << "... " << std::flush;
        evolved = this->evolveVector(initialState);
        normDiff = std::abs(1 - arma::norm(evolved));
        std::cout << "Norm difference: " << normDiff << std::endl;
    } while (normDiff > 1e-12);

    std::size_t minN = this->N / 2;
    std::size_t maxN = this->N;

    do {
        std::size_t midN = (minN + maxN) / 2;
        this->N = midN;
        std::cout << "Trying " << this->N << "... " << std::flush;
        evolved = evolveVector(initialState);
        normDiff = std::abs(1 - arma::norm(evolved));
        std::cout << "Norm difference: " << normDiff << std::endl;

        if (normDiff <= 1e-12)
            maxN = midN;
        else
            minN = midN;
    } while (maxN - minN > 1);

    std::cout << "Optimal orders needed: " << this->N << std::endl;
}

void ChebyshevEvolver::evolve() {
    // Actually this->step == this->steps - 1 here will give 1 step too much, but do not throw for convenience of use
    Assert(this->step < this->steps);
    this->step++;
    this->t += this->dt;

    this->currentState = this->evolveVector(this->currentState);
}

const arma::cx_vec &ChebyshevEvolver::getCurrentState() const {
    return this->currentState;
}
