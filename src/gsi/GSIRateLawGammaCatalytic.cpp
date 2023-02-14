/**
 * @file GSIRateLawGammaCatalytic.cpp
 *
 * @brief Class which computes the reaction rate constant for a surface
 *        reaction constant according to a gamma type model with gamma
 *        constant.
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
#include "Errors.h"
#include "Thermodynamics.h"
#include "Transport.h"

#include "GSIRateLaw.h"

using namespace Eigen;

using namespace Mutation::Utilities;
using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawGammaCatalytic : public GSIRateLaw
{
public:
    GSIRateLawGammaCatalytic(ARGS args)
        : GSIRateLaw (args),
          mv_react(args.s_reactants)
    {
        assert(args.s_node_rate_law.tag() == "gamma_finite");

        args.s_node_rate_law.getAttribute( "model", m_model,
            "The catalytic model (scott:1, zoby:2, or stewart:3) for the reaction "
            "should be provided with gamma as a function of temperature.");

        int l_diff_reac = 1;
        for (int i_reac = 0; i_reac < mv_react.size() - 1; ++i_reac) {
            if (mv_react[i_reac] != mv_react[i_reac + 1]) l_diff_reac++;
        }
        mv_pre_coef.resize(2*l_diff_reac);
        mv_active_T.resize(2*l_diff_reac);
        mv_index.resize(l_diff_reac);
        for (int i_reac = 0; i_reac < l_diff_reac; ++i_reac) {
            if(mv_react[i_reac] == m_thermo.speciesIndex("N")){
                    mv_index[i_reac] = 1;
                }else{
                    mv_index[i_reac] = 0;
                }
        }
        for(int i_reac =0;i_reac < l_diff_reac;i_reac++){
            switch (m_model){
                case 1:
                if(mv_index[i_reac] == 1){
                    mv_pre_coef[2*i_reac] = 0.0714;
                    mv_pre_coef[2*i_reac+1] = 0.0714;
                    mv_active_T[2*i_reac]=-2219.0;
                    mv_active_T[2*i_reac+1]=-2219.0;
                }else{
                    mv_pre_coef[2*i_reac] = 16.0;
                    mv_pre_coef[2*i_reac+1] = 16.0;
                    mv_active_T[2*i_reac]=-2219.0;
                    mv_active_T[2*i_reac+1]=-2219.0;
                }
                break;
                case 2:
                if(mv_index[i_reac] == 1){
                    mv_pre_coef[2*i_reac] = 0.0714;
                    mv_pre_coef[2*i_reac+1] = 0.0714;
                    mv_active_T[2*i_reac]=-2219.0;
                    mv_active_T[2*i_reac+1]=-2219.0;
                }else{
                    mv_pre_coef[2*i_reac] = 0.00941;
                    mv_pre_coef[2*i_reac+1] = 0.00941;
                    mv_active_T[2*i_reac]=-658.9;
                    mv_active_T[2*i_reac+1]=-658.9;
                }
                break;
                case 3:
                if(mv_index[i_reac] == 1){
                    mv_pre_coef[2*i_reac] = 6.1e-2;
                    mv_pre_coef[2*i_reac+1] = 6.1e-4;
                    mv_active_T[2*i_reac]=-2480.0;
                    mv_active_T[2*i_reac+1]=5090.0;
                    
                }else{
                    mv_pre_coef[2*i_reac] = 40.0;
                    mv_pre_coef[2*i_reac+1] = 39e-9;
                    mv_active_T[2*i_reac]=-11440.0;
                    mv_active_T[2*i_reac+1]=21410.0;
                }
                break;

            }
        }

        m_idx_react = -1;
        m_idx_sp = -1;
        m_stoich_coef = -1;

        mv_imp_flux_out.resize(mv_index.size());
        mv_imp_flux_per_stoich_coef.resize(mv_index.size());

    }

//==============================================================================

    ~GSIRateLawGammaCatalytic(){}

//==============================================================================

    double forwardReactionRateCoefficient(
        const VectorXd& v_rhoi, const VectorXd& v_Twall) const
    {
        VectorXd mv_gamma;
        mv_gamma.resize(mv_index.size());
        for(int i_reac =0;i_reac < mv_index.size();i_reac++){
            if(v_Twall[0] < 1580){
                mv_gamma[i_reac] = mv_pre_coef[2*i_reac]*exp(mv_active_T[2*i_reac]/v_Twall[0]);
            }else {
                mv_gamma[i_reac] = mv_pre_coef[2*i_reac+1]*exp(mv_active_T[2*i_reac+1]/v_Twall[0]);
            }

        }
        
        m_idx_react = 0;
        for(int i_g = 0; i_g < mv_index.size(); i_g++) {
            getSpeciesIndexandStoichiometricCoefficient(
                m_idx_react, m_idx_sp, m_stoich_coef);
            mv_imp_flux_out[i_g] = computeSurfaceImpingingMassFlux(
                m_idx_sp, v_rhoi, v_Twall);

            mv_imp_flux_per_stoich_coef[i_g] =
                mv_imp_flux_out[i_g]/m_stoich_coef;
            mv_imp_flux_out[i_g] =
                mv_imp_flux_per_stoich_coef[i_g]*mv_gamma[i_g];

            m_idx_react += m_stoich_coef;
        }

        return getLimitingImpingingMassFlux();
    }

private:
    mutable int m_idx_react;
    mutable int m_idx_sp;
    mutable int m_stoich_coef;
    int m_model;
    std::vector<int> mv_index;

    std::vector<double> mv_pre_coef;
    std::vector<double> mv_active_T;

    mutable std::vector<double> mv_imp_flux_out;
    mutable std::vector<double> mv_imp_flux_per_stoich_coef;

    const std::vector<int>& mv_react;

//==============================================================================
private:
/**
 * Function which provides the species indexes and stoichiometric coefficients.
 */
    inline void getSpeciesIndexandStoichiometricCoefficient(
        int idx_react, int& idx_sp, int& stoich_coef) const
    {
        idx_sp = mv_react[idx_react];
        stoich_coef = 1;
        idx_react++;

        while(idx_react < mv_react.size()) {
            if (idx_sp != mv_react[idx_react]) {
                break;
            }
            stoich_coef++;
            idx_react++;
        }
    }

//==============================================================================
/**
 * Function that returns the impinging mass flux of a species on a surface.
 */
    inline double computeSurfaceImpingingMassFlux(
        const int& idx_sp,
        const VectorXd& v_rhoi,
        const VectorXd& v_Twall) const
    {
    	const int set_state_with_rhoi_T = 1;
    	m_thermo.setState(
    	    v_rhoi.data(), v_Twall.data(), set_state_with_rhoi_T);
    	double m_sp_thermal_speed = m_transport.speciesThermalSpeed(idx_sp);

        return m_sp_thermal_speed / 4.
        	* v_rhoi(idx_sp) / m_thermo.speciesMw(idx_sp);
    }

//==============================================================================
/**
 * Function that returns the limiting impinging mass flux of the different
 * species reacting on the surface. Essential for heteronuclear reactions.
 */
    inline double getLimitingImpingingMassFlux() const
    {
        return mv_imp_flux_out[min_element(
            mv_imp_flux_per_stoich_coef.begin(),
            mv_imp_flux_per_stoich_coef.end()) -
            mv_imp_flux_per_stoich_coef.begin()];
    }
};

ObjectProvider<
    GSIRateLawGammaCatalytic, GSIRateLaw>
    gsi_rate_law_gamma_finite("gamma_finite");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
