//
// Created by pkua on 02.05.2020.
//

#include <complex>
#include <chrono>

#include "ChebyshevEvolver.h"
#include "utils/Assertions.h"

using namespace std::complex_literals;

template<typename V, typename M>
inline void csrMult(V &Ax, const V &x, const M &Adata, const std::vector<std::size_t> &Aindices, const std::vector<std::size_t> &Aindptr)
{
    for (std::size_t i = 0; i < Ax.size(); i++)
    {
        std::complex<double> Ax_i = 0.0;
        for (std::size_t dataIdx = Aindptr[i]; dataIdx < Aindptr[i + 1]; dataIdx++)
        {
            Ax_i += Adata[dataIdx] * x[Aindices[dataIdx]];
        }
        Ax[i] = Ax_i;
    }
}

/**
 * @brief Perform one Chebyshev step by this->dt using formulas from paper:
 * Many-body localization in presence of cavity mediated long-range interactions
 */
arma::cx_vec ChebyshevEvolver::evolveState(const arma::cx_vec &state) {
    arma::sp_mat hamiltonianRescaled = (this->hamiltonian - arma::speye(arma::size(this->hamiltonian)) * this->b) / this->a;
    std::vector<double> Adata(hamiltonianRescaled.values, hamiltonianRescaled.values + hamiltonianRescaled.n_nonzero);
    std::vector<std::size_t> Aindices(hamiltonianRescaled.row_indices, hamiltonianRescaled.row_indices + hamiltonianRescaled.n_nonzero);
    std::vector<std::size_t> Aindptr(hamiltonianRescaled.col_ptrs, hamiltonianRescaled.col_ptrs + hamiltonianRescaled.n_rows + 1);

    arma::cx_vec prevPrevOrder = state;
    arma::cx_vec prevOrder(arma::size(state));
    csrMult(prevOrder, state, Adata, Aindices, Aindptr);
    arma::cx_vec currentOrder(arma::size(state));
    arma::cx_vec result(arma::size(state));
    arma::cx_vec hamVec(arma::size(state));

    result = state * std::cyl_bessel_j(0, this->a * this->dt);
    result += 2. * -1i * std::cyl_bessel_j(1, this->a * this->dt) * prevOrder;

    for (std::size_t i = 2; i <= N; i++) {
        csrMult(hamVec, prevOrder, Adata, Aindices, Aindptr);
        currentOrder = 2* hamVec - prevPrevOrder;
        result += 2. * std::pow(-1i, i) * std::cyl_bessel_j(i, this->a * this->dt) * currentOrder;
        prevPrevOrder = prevOrder;
        prevOrder = currentOrder;
    }
    result *= std::exp(-1i * this->b * this->dt);

    return result;
}

void ChebyshevEvolver::prepareFor(const arma::cx_vec &initialState, double maxTime, std::size_t maxSteps_) {
    Expects(maxTime > 0);
    Expects(maxSteps_ > 0);
    Expects(initialState.size() == this->hamiltonian.n_cols);

    this->t = 0;
    this->dt = maxTime / static_cast<double>(maxSteps_);
    this->currentState = initialState;
    this->currentStep = 0;
    this->maxSteps = maxSteps_;

    this->optimizeOrder(initialState);
}

/**
 * @brief Finds the optimal order of Chebychev expansion, requiring that norm leakage per step should not be greater
 * than MAXIMAL_NORM_LEAKAGE
 */
void ChebyshevEvolver::optimizeOrder(const arma::cx_vec &initialState) {
    arma::cx_vec evolvedState;
    double normLeakage{};
    double initialNorm = arma::norm(initialState);

    this->logger << "[ChebyshevEvolver::prepareFor] Optimizing Chebyshev expansion order:" << std::endl;

    // Exponentially find the upper limit for order
    this->N = 1;
    do {
        this->N *= 2;
        this->logger << "[ChebyshevEvolver::prepareFor] Trying " << this->N << "... " << std::flush;
        evolvedState = this->evolveState(initialState);
        normLeakage = std::abs(initialNorm - arma::norm(evolvedState));
        this->logger << "Norm leakage: " << normLeakage << std::endl;
    } while (normLeakage > MAXIMAL_NORM_LEAKAGE);

    // Now, find optimal order using bisection
    std::size_t minN = this->N / 2;
    std::size_t maxN = this->N;
    do {
        std::size_t midN = (minN + maxN) / 2;
        this->N = midN;
        this->logger << "[ChebyshevEvolver::prepareFor] Trying " << this->N << "... " << std::flush;
        evolvedState = this->evolveState(initialState);
        normLeakage = std::abs(initialNorm - arma::norm(evolvedState));
        this->logger << "Norm leakage: " << normLeakage << std::endl;

        if (normLeakage <= MAXIMAL_NORM_LEAKAGE)
            maxN = midN;
        else
            minN = midN;
    } while (maxN - minN > 1);

    this->logger << "[ChebyshevEvolver::prepareFor] Optimal orders needed: " << N << std::endl;
}

/**
 * @brief Finds the highest and lowest eigenvalues of the hamiltonian which are needed in the expansion using sparse
 * matrix techniques
 */
void ChebyshevEvolver::findSpectrumRange() {
    std::size_t numEigvals = std::min<std::size_t>(MIN_EIGVAL, this->hamiltonian.n_rows - 1);
    
    arma::vec minEigval, maxEigval;
    this->logger << "[ChebyshevEvolver::findSpectrumRange] Calculating Emin... " << std::flush;
    arma::wall_clock timer;
    timer.tic();
    Assert(arma::eigs_sym(minEigval, this->hamiltonian, numEigvals, "sa"));
    this->logger << "done (" << timer.toc() << " s). Emax... " << std::flush;
    timer.tic();
    Assert(arma::eigs_sym(maxEigval, this->hamiltonian, numEigvals, "la"));
    this->logger << "done (" << timer.toc() << " s). " << std::endl;
    double Emin = minEigval.front();
    double Emax = maxEigval.back();
//    double Emin = -31.80130421278989;
//    double Emax = 306.9052833537735;
    this->a = (Emax - Emin) / 2;
    this->b = (Emax + Emin) / 2;
}

void ChebyshevEvolver::evolve() {
    // Actually this->currentStep == this->steps here will give 1 step too much, but do not throw for convenience of use
    Assert(this->currentStep <= this->maxSteps);
    this->currentStep++;
    this->t += this->dt;

    this->currentState = this->evolveState(this->currentState);
}

const arma::cx_vec &ChebyshevEvolver::getCurrentState() const {
    return this->currentState;
}

ChebyshevEvolver::ChebyshevEvolver(const arma::sp_mat &hamiltonian, std::ostream &logger)
        : hamiltonian{hamiltonian}, logger{logger}
{
    this->findSpectrumRange();
}

double ChebyshevEvolver::getDt() const {
    return this->dt;
}
