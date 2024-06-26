/**
 * @file GSIRateLawGammaT2.cpp
 *
 * @brief  Class which computes the reaction rate constant for a surface
 *        reaction constant according to a gamma type model with gamma
 *        as an exponential function of temperature.
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


#include "Thermodynamics.h"
#include "Transport.h"

#include "AutoRegistration.h"
#include "Utilities.h"

#include "GSIRateLaw.h"

using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawGammaT2 : public GSIRateLaw
{
public:
    GSIRateLawGammaT2(ARGS args)
        : GSIRateLaw(args),
          mv_react(args.s_reactants),
          pos_T_trans(0),
          idx_react(0)
    {
        assert(args.s_node_rate_law.tag() == "gamma_T2");

        args.s_node_rate_law.getAttribute( "aexp", m_aexp,
            "The pre-exponential coefficient for the reaction "
            "should be provided with gamma as a function of temperature.");
        args.s_node_rate_law.getAttribute( "aEa", m_aEa,
            "The activation energy for the reaction "
            "should be provided with gamma as a function of temperature.");
        args.s_node_rate_law.getAttribute( "acof", m_acof,
            "The const coeficient for the reaction "
            "should be provided with gamma as a function of temperature.");
        args.s_node_rate_law.getAttribute( "bexp", m_bexp,
            "The pre-exponential coefficient for the reaction "
            "should be provided with gamma as a function of temperature.");
        args.s_node_rate_law.getAttribute( "bEa", m_bEa,
            "The activation energy for the reaction "
            "should be provided with gamma as a function of temperature.");
        args.s_node_rate_law.getAttribute( "bcof", m_bcof,
            "The const coeficient for the reaction "
            "should be provided with gamma as a function of temperature.");
    }

//==============================================================================

    ~GSIRateLawGammaT2( ){ }

//==============================================================================

    double forwardReactionRateCoefficient(
        const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Twall) const
    {
    	double Twall = v_Twall(pos_T_trans);

    	const int set_state_with_rhoi_T = 1;
    	m_thermo.setState(
    	    v_rhoi.data(), v_Twall.data(), set_state_with_rhoi_T);
    	double m_sp_thermal_speed = m_transport.speciesThermalSpeed(
                                        mv_react[idx_react]);

        double alpha = (m_acof+m_aexp*std::exp(- m_aEa/Twall))/
                        (m_bcof+m_bexp*std::exp(- m_bEa/Twall));
        return  m_sp_thermal_speed/4.
                * alpha / m_thermo.speciesMw(
                mv_react[idx_react])*v_rhoi(mv_react[idx_react]);
    }

private:
    const size_t pos_T_trans;
    const size_t idx_react;

    double m_aexp;
    double m_bexp;
    double m_aEa;
    double m_bEa;
    double m_acof;
    double m_bcof;

    const std::vector<int>& mv_react;
};

ObjectProvider<
    GSIRateLawGammaT2, GSIRateLaw>
    gsi_rate_law_gamma_T2("gamma_T2");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
