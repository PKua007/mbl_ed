//
// Created by pkua on 02.05.2020.
//

#include <complex>
#include <chrono>

#include "ChebyshevEvolver.h"
#include "utils/Assertions.h"
#include "utils/OMPMacros.h"

using namespace std::complex_literals;

namespace {
    /**
     * @brief returns (matrix * vec)[elementIdx] using CSR matrix * dense vector multiplication scheme
     */
    inline std::complex<double> csrMatTimesCxVecElem(const arma::cx_vec &vec, const std::vector<double> &matrixData,
                                                     const std::vector<std::size_t> &matrixColIdx,
                                                     const std::vector<std::size_t> &matrixRowPtr,
                                                     std::size_t elementIdx)
    {
        std::complex<double> Ax_i = 0.0;
        for (std::size_t dataIdx = matrixRowPtr[elementIdx]; dataIdx < matrixRowPtr[elementIdx + 1]; dataIdx++)
            Ax_i += matrixData[dataIdx] * vec[matrixColIdx[dataIdx]];
        return Ax_i;
    }
}

/**
 * @brief Perform one Chebyshev step by this->dt.
 * @details It uses custom sparse matrix - dense vector multiplication, since this from Armadillo is lame.
 */
arma::cx_vec ChebyshevEvolver::evolveState(const arma::cx_vec &state) {
    // Rescales so that hamiltonian's eigenvalues lie in [-1, 1] range
    arma::sp_mat hamiltonianRescaled = (this->hamiltonian - arma::speye(arma::size(this->hamiltonian)) * this->b)
                                       / this->a;

    // We just steal sparse matrix data from Armadillo's matrix. Note, that it uses CSC format and we interpret it as
    // CSR, but it's ok since the matrix is Hermitean
    // details: https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format)
    std::vector<double> hamData(hamiltonianRescaled.values, hamiltonianRescaled.values + hamiltonianRescaled.n_nonzero);
    std::vector<std::size_t> hamColIdx(hamiltonianRescaled.row_indices,
                                      hamiltonianRescaled.row_indices + hamiltonianRescaled.n_nonzero);
    std::vector<std::size_t> hamRowPtr(hamiltonianRescaled.col_ptrs,
                                     hamiltonianRescaled.col_ptrs + hamiltonianRescaled.n_rows + 1);

    // Now we perform Chebyshev expansion summation as stated in paper:
    // Many-body localization in presence of cavity mediated long-range interactions
    // using Clenshaw algorithm from:
    // https://en.wikipedia.org/wiki/Clenshaw_algorithm
    // All vague variable names follow from this Wikipedia link.
    // In this case a_0 = J_0(this->a t), a_k = 2(-i)^k J_k(this->a t)
    arma::cx_vec bNext(arma::size(state), arma::fill::zeros);
    arma::cx_vec bNextNext(arma::size(state), arma::fill::zeros);
    arma::cx_vec B(arma::size(state));  // Capital b not to collide with this->b

    // Iteratively reach bNext = b_1, bNext = b_2
    for (std::size_t i = this->N; i > 0; i--) {
        std::complex<double> coeff = 2. * std::pow(-1i, i) * std::cyl_bessel_j(i, this->a * this->dt);

        _OMP_PARALLEL_FOR
        for (std::size_t j = 0; j < state.size(); j++)
            B[j] = coeff * state[j] + 2. * csrMatTimesCxVecElem(bNext, hamData, hamColIdx, hamRowPtr, j) - bNextNext[j];
        bNextNext = bNext;
        bNext = B;
    }

    // Now, compute p_n and store it in B
    std::complex<double> coeff = std::cyl_bessel_j(0, this->a * this->dt);

    _OMP_PARALLEL_FOR
    for (std::size_t i = 0; i < state.size(); i++) {
        B[i] = coeff * state[i] + csrMatTimesCxVecElem(bNext, hamData, hamColIdx, hamRowPtr, i) - bNextNext[i];
        B[i] *= std::exp(-1i * this->b * this->dt);
    }

    return B;
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

    this->logger.info() << "Optimizing Chebyshev expansion order:" << std::endl;

    // Exponentially find the upper limit for order
    this->N = 1;
    do {
        this->N *= 2;
        Assert(this->N <= 2048);
        evolvedState = this->evolveState(initialState);
        normLeakage = std::abs(initialNorm - arma::norm(evolvedState));
        Assert(!std::isnan(normLeakage));
        this->logger.info() << "Trying " << this->N << "... Norm leakage: " << normLeakage << std::endl;
    } while (normLeakage > MAXIMAL_NORM_LEAKAGE);

    // Now, find optimal order using bisection
    std::size_t minN = this->N / 2;
    std::size_t maxN = this->N;
    do {
        std::size_t midN = (minN + maxN) / 2;
        this->N = midN;
        evolvedState = this->evolveState(initialState);
        normLeakage = std::abs(initialNorm - arma::norm(evolvedState));
        Assert(!std::isnan(normLeakage));
        this->logger.info() << "Trying " << this->N << "... Norm leakage: " << normLeakage << std::endl;

        if (normLeakage <= MAXIMAL_NORM_LEAKAGE)
            maxN = midN;
        else
            minN = midN;
    } while (maxN - minN > 1);

    this->logger.info() << "Optimal orders needed: " << N << std::endl;
}

/**
 * @brief Finds the highest and lowest eigenvalues of the hamiltonian which are needed in the expansion using sparse
 * matrix techniques
 */
void ChebyshevEvolver::findSpectrumRange() {
    std::size_t numEigvals = std::min<std::size_t>(MIN_EIGVAL, this->hamiltonian.n_rows - 1);
    arma::vec minEigval, maxEigval;
    this->logger.verbose() << "Calculating spectrum bounds started... " << std::endl;

    arma::wall_clock timer;
    timer.tic();
    Assert(arma::eigs_sym(minEigval, this->hamiltonian, numEigvals, "sa"));
    double EminTime = timer.toc();

    timer.tic();
    Assert(arma::eigs_sym(maxEigval, this->hamiltonian, numEigvals, "la"));
    double EmaxTime = timer.toc();

    double Emin = minEigval.front();
    double Emax = maxEigval.back();
    this->a = (Emax - Emin) / 2;
    this->b = (Emax + Emin) / 2;
    this->logger.info() << "Calculating spectrum range done (Emin: " << EminTime << " s, ";
    this->logger << "Emax: " << EmaxTime << " s)." << std::endl;
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

ChebyshevEvolver::ChebyshevEvolver(const arma::sp_mat &hamiltonian, Logger &logger)
        : hamiltonian{hamiltonian}, logger{logger}
{
    this->findSpectrumRange();
}

double ChebyshevEvolver::getDt() const {
    return this->dt;
}
