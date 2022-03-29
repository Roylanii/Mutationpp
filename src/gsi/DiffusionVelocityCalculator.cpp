/**
 * @file DiffusionVelocityCalculator.cpp
 *
 * @brief Class which computes the diffusion velocities needed by
 *        the surface balances.
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


#include "Errors.h"
#include "Thermodynamics.h"
#include "Transport.h"

#include "DiffusionVelocityCalculator.h"

using namespace Eigen;

namespace Mutation {
    namespace GasSurfaceInteraction {

DiffusionVelocityCalculator::DiffusionVelocityCalculator(
    const Mutation::Thermodynamics::Thermodynamics& thermo,
    Mutation::Transport::Transport& transport)
        : m_ns(thermo.nSpecies()),
          m_nT(thermo.nEnergyEqns()),
          mv_dxidx(thermo.nSpecies()),
          mv_mole_frac_edge(thermo.nSpecies()),
          m_transport(transport),
          m_dx(0.),
          m_is_diff_set(false)
{ }

//==============================================================================

DiffusionVelocityCalculator::~DiffusionVelocityCalculator(){}

//==============================================================================

void DiffusionVelocityCalculator::setDiffusionModel(
    const VectorXd& v_mole_frac_edge, const double& dx)
{
    mv_mole_frac_edge = v_mole_frac_edge;

    if (dx <= 0.) {
    	throw LogicError()
        << "Calling DiffusionVelocityCalculator::setDiffusionModel() with a "
        << "distance less or equal to zero. The distance dx should always be "
           "positive.";
    }
    m_dx = dx;

    m_is_diff_set = true;
}

//==============================================================================

void DiffusionVelocityCalculator::computeDiffusionVelocities(
    const VectorXd& v_mole_frac,
    VectorXd& v_diff_velocities)
{
    if (!m_is_diff_set) {
    	throw LogicError()
        << "Calling DiffusionVelocityCalculator::computeDrivingForces() before "
        << "calling DiffusionVelocityCalculator::setDiffusionCalculator().";
    }
    mv_dxidx = (v_mole_frac - mv_mole_frac_edge)/m_dx;

    double electric_field = 0.E0;
    m_transport.stefanMaxwell(
        mv_dxidx.data(),
        v_diff_velocities.data(),
        electric_field);
}

void DiffusionVelocityCalculator::computeDiffusionVelocitiesYi(
    const VectorXd& v_mass_frac,
    VectorXd& v_diff_velocities)
{
    if (!m_is_diff_set) {
    	throw LogicError()
        << "Calling DiffusionVelocityCalculator::computeDrivingForces() before "
        << "calling DiffusionVelocityCalculator::setDiffusionCalculator().";
    }
    // mv_mole_frac_edge is species mass fraction actually
    // for simplcity, use the same variable to avoid new varaable 
    mv_dxidx = (v_mass_frac - mv_mole_frac_edge)/m_dx;
    VectorXd di(m_ns);
    m_transport.averageDiffusionCoeffs(di.data());
    v_diff_velocities = -di.cwiseProduct(mv_dxidx);
}

Eigen::MatrixXd DiffusionVelocityCalculator::computeDiffusionVelocitiesJacobian(
    const VectorXd& v_mole_frac, const double& pert_m, const double& pert_T)
{
    if (!m_is_diff_set) {
    	throw LogicError()
        << "Calling DiffusionVelocityCalculator::computeDrivingForces() before "
        << "calling DiffusionVelocityCalculator::setDiffusionCalculator().";
    }
    MatrixXd mm_jacobian(m_ns,m_ns+m_nT);
    VectorXd v_x=v_mole_frac;
    VectorXd v_f_unpert(m_ns);
    VectorXd v_f(m_ns);
    mv_dxidx = (v_mole_frac - mv_mole_frac_edge)/m_dx;

    double electric_field = 0.E0;
    m_transport.stefanMaxwell(
        mv_dxidx.data(),
        v_f_unpert.data(),
        electric_field);
    
    for (int i_s=0; i_s<m_ns;++i_s)
    {
        v_x=v_mole_frac;
        v_x(i_s) += pert_m;
        mv_dxidx = (v_x - mv_mole_frac_edge)/m_dx;
        m_transport.stefanMaxwell(
            mv_dxidx.data(),
            v_f.data(),
            electric_field);
        mm_jacobian.col(i_s) = (v_f-v_f_unpert)/pert_m;
    }

    for (int i_T=0; i_T<m_nT; ++i_T)
    {
        mv_dxidx = (v_mole_frac - mv_mole_frac_edge)/m_dx;
        m_transport.stefanMaxwell(
            mv_dxidx.data(),
            v_f.data(),
            electric_field,
            pert_T);
        mm_jacobian.col(m_ns+i_T) = (v_f-v_f_unpert)/pert_T;
    }
    return mm_jacobian;
    
}

    } // namespace GasSurfaceInteraction
} // namespace Mutation
