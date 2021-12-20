#include "input.h"
//#include <catch.hpp>
#include <Eigen/Dense>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

using namespace Mutation;
using namespace Eigen;
void getOptLong(int argc, char *argv[], std::map<int, std::string> &input)
{
    int opt;              // getopt_long() 的返回值
    int digit_optind = 0; // 设置短参数类型及是否需要参数

    int option_index = 0;
    // 设置短参数类型及是否需要参数
    const char *optstring = "a";

    // 设置长参数类型及其简写，比如 --reqarg <==>-r
    /*
    struct option {
             const char * name;  // 参数的名称
             int has_arg; // 是否带参数值，有三种：no_argument， required_argument，optional_argument
             int * flag; // 为空时，函数直接将 val 的数值从getopt_long的返回值返回出去，
                     // 当非空时，val的值会被赋到 flag 指向的整型数中，而函数返回值为0
             int val; // 用于指定函数找到该选项时的返回值，或者当flag非空时指定flag指向的数据的值
        };
    其中：
        no_argument(即0)，表明这个长参数不带参数（即不带数值，如：--name）
            required_argument(即1)，表明这个长参数必须带参数（即必须带数值，如：--name Bob）
            optional_argument(即2)，表明这个长参数后面带的参数是可选的，（即--name和--name Bob均可）
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
        // else
        // {
        //     throw LogicError()
        //     << "input argument wrong";
        // }
    }
}
int main(int argc, char *argv[])
{
    const double tol = std::numeric_limits<double>::epsilon();

    //read input from argv
    // std::map<int, std::string> input;
    // getOptLong(argc, argv, input); 
    // std::map<int, std::string>::iterator map_iter;
    // map_iter = input.find(0);
    // if (map_iter != input.end())
    //     Mutation::GlobalOptions::workingDirectory(map_iter->second);
    // map_iter = input.find(1);
    // MixtureOptions opts(map_iter->second.c_str());
    // map_iter = input.find(2);
    // double P_init = stod(map_iter->second);
    // map_iter = input.find(3);
    // double T_init = stod(map_iter->second);


    // Mixture
    //MixtureOptions opts("seb_oxidation_NASA9_ChemNonEq1T");

    //Read input 
    Input input("input.xml");
    
    Mutation::GlobalOptions::workingDirectory(input.workingDirectory());
    Mutation::GlobalOptions::dataDirectory(input.dataDirectory());
    MixtureOptions opts(input.mixturename());
    Mixture mix(opts);
    double P_init = input.pressureinit();
    double T_init = input.temperatureinit();
    // Setting up
    const size_t set_state_with_rhoi_T = 1;
    const size_t pos_T_trans = 0;
    size_t ns = mix.nSpecies();
    size_t nT = mix.nEnergyEqns();
    size_t neq = ns + nT;
    // Conditions T = 3000K and p = 100Pa
    VectorXd Teq = VectorXd::Constant(nT, T_init);

    mix.equilibrate(Teq(pos_T_trans), P_init);
    // Setting number of iterations for the solver
    const int iter = input.maxstepinit();
    mix.setIterationsSurfaceBalance(iter);

    // Mass gradient
    VectorXd xi_e(ns);
    xi_e = Map<const VectorXd>(mix.X(), ns);
    double dx = input.distanceinit();
    mix.setDiffusionModel(xi_e.data(), dx);

    // Temperature gradient
    VectorXd T_e = Teq;
    mix.setGasFourierHeatFluxModel(T_e.data(), dx);

    // Initial conditions of the surface are the ones in the first
    // physical cell
    VectorXd rhoi_s(ns);
    mix.densities(rhoi_s.data());
    VectorXd T_s = VectorXd::Constant(nT, 1800.);
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
    dTdx = (T_s - T_e) / dx;
    VectorXd lambda(nT);
    mix.frozenThermalConductivityVector(lambda.data());

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
    const double sigma = 2. * pow(PI, 5) * pow(KB, 4) / (15 * pow(C0, 2) * pow(HP, 3));
    const double eps = .86;
    VectorXd q_srad = VectorXd::Zero(nT);
    q_srad(pos_T_trans) = sigma * eps * pow(T_s(pos_T_trans), 4);

    // Chemical Energy Contribution
    VectorXd v_hi_rhoi_vi = VectorXd::Zero(nT);
    v_hi_rhoi_vi(pos_T_trans) = -v_hi.head(ns).dot(rhoi_s.cwiseProduct(vdi));

    // Building balance functions
    VectorXd F(neq);
    F.head(ns) = (rhoi_s / rho) * mblow + rhoi_s.cwiseProduct(vdi) - wdot;
    F.tail(nT) = -lambda.cwiseProduct(dTdx) - q_srad + mblow * v_h -
                 v_hi_rhoi_vi;

    // Compute error
    double err = F.lpNorm<Infinity>();
    std::cout << "Solving gas surface interaction done!" << std::endl;
    std::cout << "The L_infty norm of the specie density error and energy error is :";
    std::cout << std::setw(13) << err; // Cp [J/kg-K]
    std::cout << std::endl;
    // Write a header line for the table
    //Temperature
    std::cout << "Surface temperature[K]" <<std::endl;
    for (int j = 0; j < nT; ++j)
        std::cout << std::setw(13) << T_s(j);
    std::cout << std::endl;

    std::cout << std::setw(22) << "Blowing flux[kg/m^2-s]" << std::endl;
    std::cout << std::setw(22) << mblow;
    std::cout << std::endl;

    std::cout << "Surface species property" <<std::endl;
    std::cout <<setw(28) << "";
    for (int j = 0; j < ns; ++j)
        std::cout << std::setw(13) << mix.speciesName(j);
    std::cout << std::endl;

    //species density 
    std::cout <<setw(28) << "Species density[kg/m^3]";
    for (int j = 0; j < ns; ++j)
        std::cout << std::setw(13) << rhoi_s(j);
    std::cout << std::endl;
    std::cout <<setw(28) << "Species mole fraction";
    for (int j = 0; j < mix.nSpecies(); ++j)
        std::cout << std::setw(13) << mix.X()[j];
    std::cout << std::endl;

    //chemical source
    std::cout <<setw(28) << "Chemical mass[kg/m^2-s]";
    for (int j = 0; j < ns; ++j)
        std::cout << std::setw(13) << wdot(j);
    std::cout << std::endl;    
    
}