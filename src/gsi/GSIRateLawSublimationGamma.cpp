/**
 * @file GSIRateLawSublimationGamma.h
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

class GSIRateLawSublimationGamma : public GSIRateLaw
{
public:
    GSIRateLawSublimationGamma(ARGS args)
        : GSIRateLaw(args),
          mv_prod(args.s_products),
          set_state_with_rhoi_T(1),
          pos_T_trans(0),
          idx_gas_prod(0)
    {
        assert(args.s_node_rate_law.tag() == "subgamma");

        args.s_node_rate_law.getAttribute("pre_exp", m_pre_exp,
            "The reaction probability for the sublimation reaction "
            "should be provided.");
        args.s_node_rate_law.getAttribute( "T", m_activ_en,
            "The vapour coefficient Ps for the sublimation reaction "
            "should be provided.");
    }

//==============================================================================

    ~GSIRateLawSublimationGamma(){}

//==============================================================================

    double forwardReactionRateCoefficient(
        const VectorXd& v_rhoi, const VectorXd& v_Twall) const
    {
      double Twall = v_Twall(pos_T_trans);
      m_thermo.setState(
    	    v_rhoi.data(), v_Twall.data(), set_state_with_rhoi_T);
      double m_sp_thermal_speed = m_transport.speciesThermalSpeed(
                                        mv_prod[idx_gas_prod]);
      return  m_sp_thermal_speed/(4.*RU*Twall)
          * m_pre_exp * std::exp(- m_activ_en/Twall)/ m_thermo.speciesMw(mv_prod[idx_gas_prod]);

    }

private:
    const size_t pos_T_trans;
    const size_t idx_gas_prod;
    const int set_state_with_rhoi_T;

    double m_pre_exp;
    double m_activ_en;

    const std::vector<int>& mv_prod;
};

ObjectProvider<
    GSIRateLawSublimationGamma, GSIRateLaw>
    gsi_rate_law_sublimationgamma("subgamma");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
