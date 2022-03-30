/**
 * @file Surface.h
 *
 * @brief Declaration of the abstract Surface class.
 */

/*
 * Copyright 2014-2020 von Karman Institute for Fluid Dynamics (VKI)
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


#ifndef SURFACE_H
#define SURFACE_H

#include <Eigen/Dense>

namespace Mutation { namespace Thermodynamics { class Thermodynamics; }}
namespace Mutation { namespace Transport { class Transport; }}
namespace Mutation { namespace Utilities { namespace IO { class XmlElement; }}}

namespace Mutation {
    namespace GasSurfaceInteraction {

class SurfaceProperties;
class SurfaceState;

//==============================================================================

/**
 * Structure which stores the necessary inputs for the Surface
 * class.
 */
struct DataSurface {
    Mutation::Thermodynamics::Thermodynamics& s_thermo;
    Mutation::Transport::Transport& s_transport;
    const std::string& s_gsi_mechanism;
    const Mutation::Utilities::IO::XmlElement& xml_feats;
    const Mutation::Utilities::IO::XmlElement& xml_surf_chem;
    const Mutation::Utilities::IO::XmlElement& xml_surf_rad;
    SurfaceState& s_surf_state;
};

//==============================================================================

/**
 * This is the abstract class for the solution of the mass and energy balances
 * at the wall.
 */
class Surface
{
public:
    typedef const DataSurface& ARGS;

//==============================================================================

	/// Returns name of this type.
	static std::string typeName() { return "Surface"; }

//==============================================================================

    /**
     * Destructor
     */
    virtual ~Surface(){ }

//==============================================================================

    /**
     * Purely virtual function which returns the surface reaction rates for
     * all species in the gas mixture in \f$ kg/m^2 \f$.
     */
    virtual void computeSurfaceReactionRates(
        Eigen::VectorXd& v_surf_chem_rates) = 0;

//==============================================================================

    /**
     * Purely virtual function which returns the surface reaction rates Jacobian
     * for all species in the gas mixture in the mole fraction and the temperature$.
     */
    virtual void computeSurfaceReactionRatesJacobian(
        Eigen::MatrixXd& m_surf_chem_rates_jacobian) = 0;

//==============================================================================

    /**
     * Returns the reaction rates for all the surface reactions occuring or
     * an error if no surface reactions are considered.
     */
    virtual Eigen::VectorXd computeSurfaceReactionRatesPerReaction(){
        throw LogicError()
        << "computeSurfaceReactionRatesPerReaction can be called only "
        << "when surface reactions can occur!";
    }

//==============================================================================

    /**
     * Returns the reaction rates for all the surface reactions occuring or
     * an error if no surface reactions are considered.
     */
    virtual int nSurfaceReactions(){
        throw LogicError()
        << "nSurfaceReactions can be called only "
        << "when surface reactions can occur!";
    }
//==============================================================================

    /**
     * Function which set ups the diffusion model in order to compute
     * the gradient of mole fractions. Requires as input a mole fraction
     * vector for the chemical state of the gas near the surface and a distance
     * these mole fractions are computed in meters.
     *
     * @param v_mole_frac_edge   mole fractions at a distance from the surface
     *                           used to compute the gradient for diffusion
     * @param dx                 distance from the surface in m
     */
    virtual void setDiffusionModel(
        const Eigen::VectorXd& v_mole_frac_edge, const double& dx)
        {
            throw LogicError()
            << "setDiffusionModel can be called only when computing "
            << "the the diffusion flux to the surface are needed!";
        }

//==============================================================================
     /**
     * Function which set ups the conductive heat flux model necessary for the
     * energy balance. Requires as input the temperature vector at a distance
     * from the surface and the distance in meters.
     *
     * @param v_T_edge         temperatures at a distance from the surface
     *                         used to compute the gradient for diffusion
     * @param dx               distance from the surface in m
     */
    virtual void setGasFourierHeatFluxModel(
        const Eigen::VectorXd& v_T, const double& dx)
    {
        throw LogicError()
        << "setConductiveHeatFluxModel can be called only when solving "
        << "the surface energy balance!";
    }
    virtual void setGasFourierHeatFluxModel(const double& constHeat)
    {
        throw LogicError()
        << "setGasFourierHeatFluxModel can be called only when solving "
        << "the surface energy balance!";
    }

//==============================================================================

    virtual void setGasRadHeatFlux(const double& gas_rad_heat_flux)
    {
        throw LogicError()
        << "setGasRadHeatFlux can be called only when solving "
        << "the surface energy balance!";
    }

//==============================================================================

    /**
     * Virtual function to be called in order to solve the
     * surface balance.
     */
    virtual void solveSurfaceBalance()
    {
        throw LogicError()
        << "solveSurfaceBalance can be called only when solving "
        << "the surface balance!";
    }

//==============================================================================
// added by zhangjingchao

    /**
     * @brief solve SEB and SMB to compute balance mass fraction gradient
     *  and tenmperature gradient
     * 
     * @param p_dy mass fraction gradient
     * @param p_dT temperature gradient
     */
    virtual void solveSurfaceGradient(double*p_dy,double* p_dT)
    {
        throw LogicError()
        << "solveSurfaceGradien can be called only when solving "
        << "the surface balance!";
    }
// These virtual function are added to prevent calling from surface object
    
    /**
     * @brief Get the Surface Residual with the given surface state
     * 
     * @param p_mass input surface state(species partial density)
     * @param p_energy surface temperature
     * @param p_res in return, surface residual
     */
    virtual void getSurfaceRes(double* const p_res)
    {
        throw LogicError()
        << "getSurfaceRes can be called only when solving "
        << "the surface balance!";
    }
    
    /**
     * Purely virtual function to be called to change the number of
     * iterations performed in the surface balance.
     */
    virtual void setIterationsSurfaceBalance(const int& iter)
    {
        throw LogicError()
        << "setIterationsSurfaceBalance can be called only when solving "
        << "the surface  balance!";
    }

    virtual void setSubIterationsSurfaceBalance(const int& iter)
    {
        throw LogicError()
        << "setSubIterationsSurfaceBalance can be called only when solving "
        << "the surface balance!";
    }

    virtual double getSurfaceRadiativeHeatFlux()
    {
        throw LogicError()
        << "getSufaceRadiativeHeatFlux can be called only when solving "
        << "the surface energy balance!";
    }

    /**
     * Purely virtual function to be called to write iterational history
     * iterations performed in the surface balance.
     */
    virtual void setIterationsHistory(const bool& iof)
    {
        throw LogicError()
        << "setIterationsHistory can be called only when solving "
        << "the surface energy balance!";
    }

    /**
     * Purely virtual function to be called to change the number of
     * iterations performed in the surface balance.
     */
    virtual void setIterationsEps(const double& eps)
    {
        throw LogicError()
        << "setIterationsEps can be called only when solving "
        << "the surface energy balance!";
    }
    /**
     * Purely virtual function to be called to change the increment of
     * species density performed in the surface balance.
     */
    virtual void setIterationsPert_m(const double& pert_m)
    {
        throw LogicError()
        << "setIterationsPert_m can be called only when solving "
        << "the surface energy balance!";
    }
    /**
     * Purely virtual function to be called to change the increment of
     * temperature performed in the surface balance.
     */
    virtual void setIterationsPert_T(const double& Pert_T)
    {
        // throw LogicError()
        // << "setIterationsPert_T can be called only when solving "
        // << "the surface energy balance!";
        std::cout << "WARNING! setIterationsPert_T can be called only when \
        solving the surface energy balance! This setting is inactive!";
    }

    virtual double computeGasFourierHeatFlux(
        const Eigen::VectorXd &v_T)
    {
        throw LogicError()
            << "computeGasFourierHeatFlux can be called only when solving "
            << "the surface energy balance!";
    }
    virtual void comSurfaceDiffusionVelocity(const Eigen::VectorXd& v_x, double* vdi)
    {
        {
        throw LogicError()
            << "comSurfaceDiffusionVelocity can be called only when solving "
            << "the surface energy balance!";
    }
    }

    virtual double computeGasFourierHeatFlux()
    {
        throw LogicError()
            << "computeGasFourierHeatFlux can be called only when solving "
            << "the surface energy balance!";
    }

    virtual double computeSolidHeat()
    {
        throw LogicError()
            << "computeSolidHeat can be called only when solving "
            << "the surface energy balance!";
    }

//==============================================================================

    /**
     * Purely virtual function returning the total mass blowing flux
     * due to surface and bulk phase processes.
     */
    virtual double massBlowingRate() = 0;

//==============================================================================

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // SURFACE_H
