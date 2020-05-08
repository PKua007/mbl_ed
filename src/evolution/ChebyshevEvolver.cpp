//
// Created by pkua on 02.05.2020.
//

#include <complex>

#include "ChebyshevEvolver.h"
#include "utils/Assertions.h"

using namespace std::complex_literals;

ChebyshevEvolver::ChebyshevEvolver(arma::sp_mat hamiltonian, double Nfactor) : hamiltonian{std::move(hamiltonian)}, Nfactor{Nfactor} {
    Expects(Nfactor > 0);
}

void ChebyshevEvolver::prepareFor(const arma::cx_vec &initialState, double tMax,
                                  std::size_t steps) {
    Expects(tMax > 0);
    Expects(steps >= 2);
    Expects(initialState.size() == this->hamiltonian.n_cols);

    this->t = 0;
    this->dt = tMax / static_cast<double>(steps - 1);
    this->currentState = initialState;

    arma::vec minEigval, maxEigval;
    Assert(arma::eigs_sym(minEigval, this->hamiltonian, 5, "sa"));
    Assert(arma::eigs_sym(maxEigval, this->hamiltonian, 5, "la"));
    double Emin = minEigval.front();
    double Emax = maxEigval.back();
    this->a = (Emax - Emin) / 2;
    this->b = (Emax + Emin) / 2;

    arma::sp_mat hamiltonianRescaled = (this->hamiltonian - arma::speye(arma::size(this->hamiltonian)) * this->b) / this->a;
    this->N = static_cast<std::size_t>(this->Nfactor * 2 * this->a * tMax);

    this->chebyshevVectors.resize(N + 1);
    this->chebyshevVectors[0] = initialState;
    this->chebyshevVectors[1] = hamiltonianRescaled * initialState;
    for (std::size_t i = 2; i <= N; i++)
        this->chebyshevVectors[i] = 2 * hamiltonianRescaled * this->chebyshevVectors[i - 1] - this->chebyshevVectors[i - 2];
}

void ChebyshevEvolver::evolve() {
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

double ChebyshevEvolver::getCurrentTime() const {
    return this->t;
}
