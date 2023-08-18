/**
 * @file GSIRateLawSublimationKL.h
 *
 * @brief Class which computes the reaction rate constant for a sublimation
 *        surface reaction using Knudsen-Langmuir.
 */

/*
 * Copyright 2018-2020 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */


#include "AutoRegistration.h"
#include "Thermodynamics.h"
#include "Transport.h"
#include "Utilities.h"

#include "GSIRateLaw.h"

using namespace Eigen;
using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawSublimationKL : public GSIRateLaw
{
public:
    GSIRateLawSublimationKL(ARGS args)
        : GSIRateLaw(args),
          mv_prod(args.s_products),
          set_state_with_rhoi_T(1),
          pos_T_trans(0),
          idx_gas_prod(0)
    {
        assert(args.s_node_rate_law.tag() == "sublimationkl");

        args.s_node_rate_law.getAttribute("alpha", m_alpha,
            "The reaction probability for the sublimation reaction "
            "should be provided.");
        args.s_node_rate_law.getAttribute( "p_vap", m_ps,
            "The vapour coefficient Ps for the sublimation reaction "
            "should be provided.");
        args.s_node_rate_law.getAttribute("q_vap", m_qs,
            "The vapour coefficient Qs for the sublimation reaction "
            "should be provided.");
    }

//==============================================================================

    ~GSIRateLawSublimationKL(){}

//==============================================================================

    double forwardReactionRateCoefficient(
        const VectorXd& v_rhoi, const VectorXd& v_Twall) const
    {
    	double Twall = v_Twall(pos_T_trans);

        // double sat_vap_p = std::exp(-m_ps/Twall+m_qs);
        // m_thermo.setState(
        //     v_rhoi.data(), v_Twall.data(), set_state_with_rhoi_T);
        // double ps = v_rhoi(mv_prod[idx_gas_prod]) * RU * Twall/m_thermo.speciesMw(mv_prod[idx_gas_prod])/101325.0;

        // double sp_thermal_speed = m_transport.speciesThermalSpeed2(
        //                               mv_prod[idx_gas_prod]);
        // return (sat_vap_p - ps)*m_alpha*sp_thermal_speed / m_thermo.speciesMw(mv_prod[idx_gas_prod]);
        double sat_vap_p = std::exp(-m_ps/Twall+m_qs)*101325.0;
        double sat_vap_rho = sat_vap_p * m_thermo.speciesMw(
                                 mv_prod[idx_gas_prod])/(RU*Twall);
        m_thermo.setState(
            v_rhoi.data(), v_Twall.data(), set_state_with_rhoi_T);
        double sp_thermal_speed = m_transport.speciesThermalSpeed(
                                      mv_prod[idx_gas_prod]);

        return (sat_vap_rho - v_rhoi(mv_prod[idx_gas_prod]))*m_alpha*
                   sp_thermal_speed/4.
                   / m_thermo.speciesMw(mv_prod[idx_gas_prod]) / 101325.0;
    }

private:
    const size_t pos_T_trans;
    const size_t idx_gas_prod;
    const int set_state_with_rhoi_T;

    double m_alpha;
    double m_qs;
    double m_ps;

    const std::vector<int>& mv_prod;
};

ObjectProvider<
    GSIRateLawSublimationKL, GSIRateLaw>
    gsi_rate_law_sublimationkl("sublimationkl");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
