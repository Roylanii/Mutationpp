#include "input.h"
//#include <catch.hpp>
#include <Eigen/Dense>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

using namespace Mutation;
using namespace Eigen;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities;
void getOptLong(int argc, char *argv[], std::map<int, std::string> &input)
{
    int opt;              // getopt_long() 的返回值
    int digit_optind = 0; // 设置短参数类型及是否需要参数

    int option_index = 0;
    // 设置短参数类型及是否需要参数
    const char *optstring = "a";
    /*
    struct option {
             const char * name;  // 参数的名称
             int has_arg; // 是否带参数值，有三种：no_argument， required_argument，optional_argument
             int * flag; // 为空时，函数直接将 val 的数值从getopt_long的返回值返回出去，
                         当非空时，val的值会被赋到 flag 指向的整型数中，而函数返回值为0
             int val; // 用于指定函数找到该选项时的返回值，或者当flag非空时指定flag指向的数据的值
        };
     */
    static struct option long_options[] = {
        {"workdir", required_argument, NULL, 'a'},
        {"mixture", required_argument, NULL, 'b'},
        {"pressure", required_argument, NULL, 'c'},
        {"temperature", required_argument, NULL, 'd'},
        {0, 0, 0, 0} // 添加 {0, 0, 0, 0} 是为了防止输入空值
    };

    while ((opt = getopt_long(argc,
                              argv,
                              optstring,
                              long_options,
                              &option_index)) != -1)
    {

        // printf("opt = %c\n", opt);                           // 命令参数，亦即 -a -b -n -r
        // printf("optarg = %s\n", optarg);                     // 参数内容
        // printf("optind = %d\n", optind);                     // 下一个被处理的下标值
        // printf("argv[optind - 1] = %s\n", argv[optind - 1]); // 参数内容
        // printf("option_index = %d\n", option_index);         // 当前打印参数的下标值
        // printf("\n");
        if (opt == 'a')
            input.insert(std::pair<int, std::string>(0, optarg));
        else if (opt == 'b')
            input.insert(std::pair<int, std::string>(1, optarg));
        else if (opt == 'c')
            input.insert(std::pair<int, std::string>(2, optarg));
        else if (opt == 'd')
            input.insert(std::pair<int, std::string>(3, optarg));
    }
}
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
    mix.setIterationsSurfaceBalance(input.maxstepinit());
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
    // xi_e = Map<const VectorXd>(mix.X(), ns);
    mix.getSpeciesComposition(std::string("Gas"), xi_e.data(),Composition::MOLE);
    const double dx = input.distanceinit();
    mix.setDiffusionModel(xi_e.data(), dx);

    // Temperature gradient
    VectorXd T_e = Teq;
    if (mix.getGSIMechanism() == "phenomenological_mass_energy")
        mix.setGasFourierHeatFluxModel(T_e.data(), dx);

    // Initial conditions of the surface are the ones in the first
    // physical cell
    VectorXd rhoi_s(ns);
    mix.densities(rhoi_s.data());
    VectorXd T_s = VectorXd::Constant(nT, (T_init / 2.0));
    mix.setSurfaceState(rhoi_s.data(), T_s.data(), set_state_with_rhoi_T);

    // Solve balance and request solution
    mix.solveSurfaceBalance();
    mix.getSurfaceState(rhoi_s.data(), T_s.data(), set_state_with_rhoi_T);
    double rho = rhoi_s.sum();

    // Verifying the solution gives low residual in the balance equations
    mix.setState(rhoi_s.data(), T_s.data(), set_state_with_rhoi_T);
    VectorXd xi_s(ns);
    xi_s = Map<const VectorXd>(mix.X(), ns);

    // Compute diffusion velocities
    VectorXd dxidx(ns);
    dxidx = (xi_s - xi_e) / dx;
    VectorXd vdi(ns);
    double E = 0.;
    mix.stefanMaxwell(dxidx.data(), vdi.data(), E);

    // Conductive heat flux
    VectorXd dTdx(nT);
    VectorXd lambda(nT);
    if (mix.getGSIMechanism() == "phenomenological_mass_energy")
    {
        dTdx = (T_s - T_e) / dx;
        mix.frozenThermalConductivityVector(lambda.data());
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
    v_h(pos_T_trans) = (rhoi_s / rho).dot(v_hi.head(ns));

    // Surface radiation
    VectorXd q_srad = VectorXd::Zero(nT);
    if (mix.getGSIMechanism() == "phenomenological_mass_energy")
        q_srad(pos_T_trans) = mix.getSurfaceRadiativeHeatFlux();

    // Chemical Energy Contribution
    VectorXd v_hi_rhoi_vi = VectorXd::Zero(nT);
    v_hi_rhoi_vi(pos_T_trans) = -v_hi.head(ns).dot(rhoi_s.cwiseProduct(vdi));

    //solid conduction
    MixtureOptions graphiteopt("graphite.xml");
    Mixture graphite(graphiteopt);
    const int set_state_PT = 1;
    graphite.setState(&P_init, &T_s[0], 1);
    double hcp = graphite.mixtureHMass(T_s[0]);

    // Building balance functions
    VectorXd F(ns);
    if (mix.getGSIMechanism() == "phenomenological_mass_energy")
        F.resize(neq);
    F.head(ns) = (rhoi_s / rho) * mblow + rhoi_s.cwiseProduct(vdi) - wdot;
    if (mix.getGSIMechanism() == "phenomenological_mass_energy")
        F.tail(nT) = -lambda.cwiseProduct(dTdx) - q_srad + mblow * v_h -
                     v_hi_rhoi_vi;

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
    std::cout << std::setw(24) << "Conductive heat[J]";
    std::cout << std::setw(24) << "Radiative heat[J]";
    std::cout << std::setw(24) << "Blowing heat[J]";
    std::cout << std::setw(24) << "Diffusion heat[J]";
    std::cout << std::endl;
    for (int j = 0; j < nT; ++j)
    {
        std::cout << std::setw(24) << (-lambda(j)*(dTdx(j)));
        std::cout << std::setw(24) << (-q_srad(j));
        std::cout << std::setw(24) << (mblow * v_h(j));
        std::cout << std::setw(24) << (-v_hi_rhoi_vi(j));
        std::cout << std::endl;
    }


}