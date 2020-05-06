//
// Created by pkua on 02.05.2020.
//

#include <cmath>
#include <complex>

#include <catch2/catch.hpp>

#include "simulation/FockBaseGenerator.h"
#include "simulation/HamiltonianGenerator.h"
#include "simulation/terms/HubbardHop.h"
#include "simulation/terms/HubbardOnsite.h"

using namespace std::complex_literals;

TEST_CASE("ChebycshevEvolver test") {
    FockBaseGenerator generator;
    auto basis = std::shared_ptr(generator.generate(3, 3));
    HamiltonianGenerator hamiltonianGenerator(basis, false);
    hamiltonianGenerator.addHoppingTerm(std::make_unique<HubbardHop>(1));
    hamiltonianGenerator.addDiagonalTerm(std::make_unique<HubbardOnsite>(2));

    auto H = hamiltonianGenerator.generate();

    arma::mat eigvec;
    arma::vec eigval;

    REQUIRE(arma::eig_sym(eigval, eigvec, arma::mat(H)));

    const double t = 2;

    arma::cx_vec diagonalEvolution = arma::exp(-1i * t * eigval);
    arma::cx_mat fockBasisEvolution = eigvec * arma::diagmat(diagonalEvolution) * eigvec.t();

    arma::cx_vec psi0(basis->size(), arma::fill::zeros);
    psi0[0] = 1;
    arma::cx_vec psiDiag = fockBasisEvolution * psi0;

    std::cout << "psi(0): ";
    psi0.as_row().raw_print(std::cout);
    std::cout << std::endl;
    std::cout << "diag psi(" << t << "): ";
    psiDiag.as_row().raw_print(std::cout);
    std::cout << std::endl;

    arma::vec sp_eig_min, sp_eig_max;
    arma::eigs_sym(sp_eig_min, H, 1, "sa");
    arma::eigs_sym(sp_eig_max, H, 1, "la");
    double Emin = sp_eig_min[0];
    double Emax = sp_eig_max[0];

    double a = (Emax - Emin) / 2;
    double b = (Emax + Emin) / 2;
    H = (H - arma::speye(arma::size(H))*b)/a;

    auto N = static_cast<std::size_t>(1.5 * 2 * a * t);

    std::vector<arma::cx_vec> chebVec(N + 1);
    chebVec[0] = psi0;
    chebVec[1] = H * psi0;
    for (std::size_t i = 2; i <= N; i++)
        chebVec[i] = 2 * H * chebVec[i - 1] - chebVec[i - 2];

    arma::cx_vec psiCheb(arma::size(psi0), arma::fill::zeros);
    for (std::size_t k = 1; k <= N; k++)
        psiCheb += std::pow(-1i, k) * std::cyl_bessel_j(k, a * t) * chebVec[k];
    psiCheb *= 2;
    psiCheb += psi0 * std::cyl_bessel_j(0, a * t);
    psiCheb *= std::exp(-1i * b * t);

    std::cout << "cheb psi(" << t << "): ";
    psiCheb.as_row().raw_print(std::cout);
    std::cout << std::endl << std::endl;

    std::cout << "diff: " << arma::norm(psiCheb - psiDiag) << std::endl;
}