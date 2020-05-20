//
// Created by pkua on 02.05.2020.
//

#include <complex>

#include "ChebyshevEvolver.h"
#include "utils/Assertions.h"

using namespace std::complex_literals;

ChebyshevEvolver::ChebyshevEvolver(const arma::sp_mat &hamiltonian, std::size_t N, double Nfactor)
        : hamiltonian{hamiltonian}, N{N}, Nfactor{Nfactor}
{
    Expects(Nfactor > 0);
}

void ChebyshevEvolver::rebuildChebychevVectors(const arma::cx_vec &initialState) {
    arma::sp_mat hamiltonianRescaled = (this->hamiltonian - arma::speye(arma::size(this->hamiltonian)) * this->b) / this->a;

    this->chebyshevVectors.resize(N + 1);
    this->chebyshevVectors[0] = initialState;
    this->chebyshevVectors[1] = hamiltonianRescaled * initialState;
    for (std::size_t i = 2; i <= N; i++) {
        std::cout << "Building order " << i << std::endl;
        this->chebyshevVectors[i] =
                2 * hamiltonianRescaled * this->chebyshevVectors[i - 1] - this->chebyshevVectors[i - 2];
    }
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

//    arma::vec minEigval, maxEigval;
//    std::cout << "sa... " << std::flush;
//    arma::wall_clock timer;
//    timer.tic();
//    Assert(arma::eigs_sym(minEigval, this->hamiltonian, 20, "sa"));
//    std::cout << "done (" << timer.toc() << "). la... " << std::flush;
//    timer.tic();
//    Assert(arma::eigs_sym(maxEigval, this->hamiltonian, 20, "la"));
//    std::cout << "done (" << timer.toc() << "). " << std::flush;
    //double Emin = minEigval.front();
    //double Emax = maxEigval.back();
    double Emin = -31.80130421278989;
    double Emax = 306.9052833537735;
    //std::cout.precision(16);
    //std::cout << Emin << " " << Emax << ". " << std::flush;
    this->a = (Emax - Emin) / 2;
    this->b = (Emax + Emin) / 2;

    this->maxTimeForN = this->N / (this->Nfactor * this->a * 2);
    std::cout << "maxTime: " << this->maxTimeForN << std::endl;
    Expects(this->maxTimeForN > this->dt);

    this->rebuildChebychevVectors(initialState);
}

void ChebyshevEvolver::evolve() {
    // Actually this->step == this->steps - 1 here will give 1 step too much, but do not throw for convenience of use
    Assert(this->step < this->steps);
    this->step++;

    if (this->t + this->dt > this->maxTimeForN) {
        this->t = 0;
        this->rebuildChebychevVectors(this->currentState);
    }

    this->t += this->dt;

    this->currentState.fill(arma::fill::zeros);
    for (std::size_t k = 1; k <= this->N; k++)
        this->currentState += std::pow(-1i, k) * std::cyl_bessel_j(k, this->a * this->t) * this->chebyshevVectors[k];
    this->currentState *= 2;
    this->currentState += this->chebyshevVectors[0] * std::cyl_bessel_j(0, this->a * this->t);
    this->currentState *= std::exp(-1i * this->b * this->t);
}

const arma::cx_vec &ChebyshevEvolver::getCurrentState() const {
    return this->currentState;
}

double ChebyshevEvolver::getMaxTimeForN() const {
    return this->maxTimeForN;
}
