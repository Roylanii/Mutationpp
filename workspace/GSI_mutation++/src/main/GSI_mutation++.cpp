#include "input.h"
//#include <catch.hpp>
#include <Eigen/Dense>
#include <stdlib.h>
#include <stdio.h>

using namespace Mutation;
using namespace Eigen;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities;

int main(int argc, char *argv[])
{
    const double tol = std::numeric_limits<double>::epsilon();

    // Mixture
    //MixtureOptions opts("seb_oxidation_NASA9_ChemNonEq1T");

    //Read input
    Input input("input.xml");

    Mutation::GlobalOptions::workingDirectory(input.workingDirectory());
    Mutation::GlobalOptions::dataDirectory(input.dataDirectory());
    MixtureOptions opts(input.mixturename());
    Mixture mix(opts);
    const double P_init = input.pressureinit();
    const double T_init = input.temperatureinit();
    // Setting number of iterations for the solver
    mix.setIterationsSurfaceBalance(input.maxiterationsinit());
    mix.setSubIterationsSurfaceBalance(input.subiterationsinit());
    mix.setIterationsPert_m(input.pertminit());
    mix.setIterationsPert_T(input.pertTinit());
    mix.setIterationsEps(input.tolinit());
    bool iNewtonhistory;
    std::istringstream(input.iNewtonhistory()) >> std::boolalpha >> iNewtonhistory;
    mix.setIterationsHistory(iNewtonhistory);
    // Setting up
    const size_t set_state_with_rhoi_T = 1;
    const size_t pos_T_trans = 0;
    size_t ns = mix.nSpecies();
    size_t nT = mix.nEnergyEqns();
    size_t neq = ns + nT;
    // Conditions T = 3000K and p = 100Pa
    VectorXd Teq = VectorXd::Constant(nT, T_init);

    mix.equilibrate(Teq(pos_T_trans), P_init);

    // Mass gradient
    VectorXd xi_e(ns);
    VectorXd xi_ey(ns);
    const double dx = input.distanceinit();
    xi_ey = Map<const VectorXd>(mix.Y(), ns);
    xi_e = Map<const VectorXd>(mix.X(), ns);
    // std::cout <<xi_e<<xi_ey;
    //mix.getSpeciesComposition(std::string("Gas"), xi_e.data(),Composition::MOLE);
    mix.setDiffusionModel(xi_ey.data(), dx);

    // Temperature gradient
    VectorXd T_e = Teq;
    if (mix.getGSIMechanism() == "phenomenological_mass_energy")
        mix.setGasFourierHeatFluxModel(T_e.data(), dx);
        // mix.setGasFourierHeatFluxModel( 6000000.);

    // Initial conditions of the surface are the ones in the first
    // physical cell
    VectorXd rhoi_s(ns);
    mix.densities(rhoi_s.data());
    VectorXd T_s = VectorXd::Constant(nT, (T_init ));
    mix.setSurfaceState(rhoi_s.data(), T_s.data(), set_state_with_rhoi_T);

    // Solve balance and request solution
    mix.solveSurfaceBalance();
    mix.getSurfaceState(rhoi_s.data(), T_s.data(), set_state_with_rhoi_T);
    double rho = rhoi_s.sum();

    // Verifying the solution gives low residual in the balance equations
    mix.setState(rhoi_s.data(), T_s.data(), set_state_with_rhoi_T);
    VectorXd xi_s(ns);
    // xi_s = Map<const VectorXd>(mix.X(), ns);
    // Compute diffusion velocities
    
    VectorXd vdi(ns);
    // VectorXd dxidx(ns);
    // dxidx = (xi_s - xi_e) / dx;
    // double E = 0.;
    // mix.stefanMaxwell(dxidx.data(), vdi.data(), E);
    // std::cout <<vdi;
    // mix.setDiffusionModel(xi_ey.data(), dx);
    xi_s =rhoi_s/rhoi_s.sum();
    mix.comSurfaceDiffusionVelocity(xi_s.data(),vdi.data());

    // Conductive heat flux
    VectorXd gas_conduction(nT);
    if (mix.getGSIMechanism() == "phenomenological_mass_energy")
    {
        
        gas_conduction(pos_T_trans) = mix.computeGasFourierHeatFlux(T_s.data());
        // VectorXd dTdx(nT);
        // VectorXd lambda(nT);
        // dTdx = (T_s - T_e) / dx;
        // mix.frozenThermalConductivityVector(lambda.data());
        // std::cout<<mix.computeGasFourierHeatFlux(T_s.data());
    }

    // Get surface production rates
    VectorXd wdot(ns);
    mix.setSurfaceState(rhoi_s.data(), T_s.data(), set_state_with_rhoi_T);
    mix.surfaceReactionRates(wdot.data());

    // Blowing flux (should be zero for catalysis)
    double mblow;
    mix.getMassBlowingRate(mblow);

    // Species and mixture enthalpies
    VectorXd v_hi(ns * nT);
    mix.getEnthalpiesMass(v_hi.data());
    VectorXd v_h(nT);
    v_h(pos_T_trans) = mix.mixtureHMass();//(rhoi_s / rho).dot(v_hi.head(ns));

    // Surface radiation
    VectorXd q_srad = VectorXd::Zero(nT);
    if (mix.getGSIMechanism() == "phenomenological_mass_energy")
        q_srad(pos_T_trans) = mix.getSurfaceRadiativeHeatFlux();

    // Chemical Energy Contribution
    VectorXd v_hi_rhoi_vi = VectorXd::Zero(nT);
    v_hi_rhoi_vi(pos_T_trans) = v_hi.head(ns).dot(rhoi_s.cwiseProduct(vdi));

    //solid conduction
    VectorXd solid_conduction(nT);
    solid_conduction(pos_T_trans) = mix.computeSolidHeat();

    // Building balance functions
    VectorXd F(ns);
    if (mix.getGSIMechanism() == "phenomenological_mass_energy")
        F.resize(neq);
    F.head(ns) = (rhoi_s / rho) * mblow + rhoi_s.cwiseProduct(vdi) - wdot;
    if (mix.getGSIMechanism() == "phenomenological_mass_energy")
        F.tail(nT) = gas_conduction - q_srad + mblow * v_h +
                     v_hi_rhoi_vi + solid_conduction;

    // Compute error and write
    //double err = F.lpNorm<Infinity>();
    std::cout << std::endl;
    std::cout << "The residual of gas surface interaction:" << std::endl;
    std::cout << std::setw(18) << "Res.AVG(mass)";
    std::cout << std::setw(18) << "Res.MAX(mass)";
    std::cout << std::setw(18) << "Res.AVG(energy)";
    std::cout << std::setw(18) << "Res.MAX(energy)";
    std::cout << std::endl;
    std::cout << std::setw(18) << F.head(ns).lpNorm<1>();
    std::cout << std::setw(18) << F.head(ns).lpNorm<Infinity>();
    std::cout << std::setw(18) << F.tail(nT).lpNorm<1>();
    std::cout << std::setw(18) << F.tail(nT).lpNorm<Infinity>();
    std::cout << std::endl;
    // Write a header line for the table
    //Temperature
    std::cout << std::setw(28) << "Surface balance temperature[K]" << std::endl;
    for (int j = 0; j < nT; ++j)
        std::cout << std::setw(28) << T_s(j);
    std::cout << std::endl;

    std::cout << std::setw(28) << "Surface blowing flux[kg/m^2-s]" << std::endl;
    std::cout << std::setw(28) << mblow;
    std::cout << std::endl;

    std::cout << "Surface mass properties" << std::endl;
    std::cout << setw(30) << "";
    for (int j = 0; j < ns; ++j)
        std::cout << std::setw(16) << mix.speciesName(j);
    std::cout << std::endl;

    //species density
    std::cout << setw(30) << "Species density[kg/m^3]";
    for (int j = 0; j < ns; ++j)
        std::cout << std::setw(16) << rhoi_s(j);
    std::cout << std::endl;
    //chemical source
    std::cout << setw(30) << "Chemical production[kg/m^2-s]";
    for (int j = 0; j < ns; ++j)
        std::cout << std::setw(16) << wdot(j);
    std::cout << std::endl;
    std::cout << setw(30) << "Surface species mole fraction";
    for (int j = 0; j < mix.nSpecies(); ++j)
        std::cout << std::setw(16) << mix.X()[j];
    std::cout << std::endl;
    std::cout << setw(30) << "Edge species mole fraction";
    for (int j = 0; j < mix.nSpecies(); ++j)
        std::cout << std::setw(16) << xi_e[j];
    std::cout << std::endl;

    cout.precision(4);
    //设定后续以科学计数法的方式输出浮点数
    cout.setf(ios::scientific);

    std::cout << "Surface energy properties" << std::endl;
    std::cout << std::setw(22) << "Gas conductive heat[J]";
    std::cout << std::setw(18) << "Radiative heat[J]";
    std::cout << std::setw(18) << "Blowing heat[J]";
    std::cout << std::setw(18) << "Diffusion heat[J]";
    std::cout << std::setw(26) << "Solid conductive heat[J]";
    std::cout << std::endl;
    for (int j = 0; j < nT; ++j)
    {
        std::cout << std::setw(22) << (gas_conduction(j));
        std::cout << std::setw(18) << (-q_srad(j));
        std::cout << std::setw(18) << (mblow * v_h(j));
        std::cout << std::setw(18) << (-v_hi_rhoi_vi(j));
        std::cout << std::setw(26) << (mblow * solid_conduction(j));
        std::cout << std::endl;
    }


}