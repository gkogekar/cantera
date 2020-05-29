//! @file FlowBase1D.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/FlowBase1D.h"
#include "cantera/base/ctml.h"
#include "cantera/transport/TransportBase.h"
#include "cantera/numerics/funcs.h"

using namespace std;

namespace Cantera
{

FlowBase1D::FlowBase1D(ThermoPhase* ph, size_t nsp, size_t points) :
    m_press(-1.0),
    m_nsp(nsp),
    m_thermo(0),
    m_kin(0),
    m_trans(0),
    m_epsilon_left(0.0),
    m_epsilon_right(0.0),
    m_do_soret(false),
    m_do_multicomponent(false),
    m_do_radiation(false),
    m_kExcessLeft(0),
    m_kExcessRight(0),
    m_zfixed(Undef),
    m_tfixed(-1.)
{
	if (ph->type() == "IdealGas") {
        m_thermo = static_cast<IdealGasPhase*>(ph);
    } else {
        throw CanteraError("FlowBase1D::FlowBase1D",
                           "Unsupported phase type: need 'IdealGasPhase'");
    }
    m_type = cFlowType;
    m_points = points;

    if (ph == 0) {
        return; // used to create a dummy object
    }

    // make a local copy of the species molecular weight vector
    m_wt = m_thermo->molecularWeights();

    // enable all species equations by default
    m_do_species.resize(m_nsp, true);

    // but turn off the energy equation at all points
    m_do_energy.resize(m_points,false);

    m_diff.resize(m_nsp*m_points);
    m_multidiff.resize(m_nsp*m_nsp*m_points);
    m_flux.resize(m_nsp,m_points);
    m_wdot.resize(m_nsp,m_points, 0.0);
    m_ybar.resize(m_nsp);
    m_qdotRadiation.resize(m_points, 0.0);

    vector_fp gr;
    for (size_t ng = 0; ng < m_points; ng++) {
        gr.push_back(1.0*ng/m_points);
    }
    setupGrid(m_points, gr.data());

    // Find indices for radiating species
    m_kRadiating.resize(2, npos);
    m_kRadiating[0] = m_thermo->speciesIndex("CO2");
    m_kRadiating[1] = m_thermo->speciesIndex("H2O");
}

void FlowBase1D::resize(size_t ncomponents, size_t points)
{
    Domain1D::resize(ncomponents, points);
	m_rho.resize(m_points, 0.0);
    m_wtm.resize(m_points, 0.0);
    m_cp.resize(m_points, 0.0);
    m_visc.resize(m_points, 0.0);
    m_tcon.resize(m_points, 0.0);
	m_diff.resize(m_nsp*m_points);
    if (m_do_multicomponent) {
        m_multidiff.resize(m_nsp*m_nsp*m_points);
        m_dthermal.resize(m_nsp, m_points, 0.0);
    }
    m_flux.resize(m_nsp,m_points);
    m_wdot.resize(m_nsp,m_points, 0.0);
    m_do_energy.resize(m_points,false);
    m_qdotRadiation.resize(m_points, 0.0);
    m_fixedtemp.resize(m_points);
    m_dz.resize(m_points-1);
    m_z.resize(m_points);
}

void FlowBase1D::setupGrid(size_t n, const doublereal* z)
{
	resize(m_nv, n);
	m_z[0] = z[0];
    for (size_t j = 1; j < m_points; j++) {
        if (z[j] <= z[j-1]) {
            throw CanteraError("FlowBase1D::setupGrid",
                               "grid points must be monotonically increasing");
        }
        m_z[j] = z[j];
        m_dz[j-1] = m_z[j] - m_z[j-1];
    }
}

void FlowBase1D::setTransport(Transport& trans)
{
    m_trans = &trans;
    m_do_multicomponent = (m_trans->transportType() == "Multi");

    m_diff.resize(m_nsp*m_points);
    if (m_do_multicomponent) {
        m_multidiff.resize(m_nsp*m_nsp*m_points);
        m_dthermal.resize(m_nsp, m_points, 0.0);
    }
}

void FlowBase1D::resetBadValues(double* xg)
{
	double* x = xg + loc();
	for (size_t j = 0; j < m_points; j++) {
		double* Y = x + m_nv * j + c_offset_Y;
		m_thermo->setMassFractions(Y);
		m_thermo->getMassFractions(Y);
	}
}

void FlowBase1D::_getInitialSoln(double* x)
{
	for (size_t j = 0; j < m_points; j++) {
		T(x, j) = m_thermo->temperature();
		m_thermo->getMassFractions(&Y(x, 0, j));
	}
}

void FlowBase1D::_finalize(const doublereal* x)
{
	if (!m_do_multicomponent && m_do_soret) {
		throw CanteraError("FlowBase1D::_finalize",
			"Thermal diffusion (the Soret effect) is enabled, and requires "
			"using a multicomponent transport model.");
	}

	size_t nz = m_zfix.size();
	bool e = m_do_energy[0];
	for (size_t j = 0; j < m_points; j++) {
		if (e || nz == 0) {
			m_fixedtemp[j] = T(x, j);
		}
		else {
			double zz = (z(j) - z(0)) / (z(m_points - 1) - z(0));
			double tt = linearInterp(zz, m_zfix, m_tfix);
			m_fixedtemp[j] = tt;
		}
	}
	if (e) {
		solveEnergyEqn();
	}

	if (domainType() == cFreeFlow) {
		// If the domain contains the temperature fixed point, make sure that it
		// is correctly set. This may be necessary when the grid has been modified
		// externally.
		if (m_tfixed != Undef) {
			for (size_t j = 0; j < m_points; j++) {
				if (z(j) == m_zfixed) {
					return; // fixed point is already set correctly
				}
			}

			for (size_t j = 0; j < m_points - 1; j++) {
				// Find where the temperature profile crosses the current
				// fixed temperature.
				if ((T(x, j) - m_tfixed) * (T(x, j + 1) - m_tfixed) <= 0.0) {
					m_tfixed = T(x, j + 1);
					m_zfixed = z(j + 1);
					return;
				}
			}
		}
	}
}

void FlowBase1D::updateProperties(size_t jg, double* x, size_t jmin, size_t jmax)
{
	// properties are computed for grid points from j0 to j1
	size_t j0 = std::max<size_t>(jmin, 1) - 1;
	size_t j1 = std::min(jmax + 1, m_points - 1);

	updateThermo(x, j0, j1);
	if (jg == npos || m_force_full_update) {
		// update transport properties only if a Jacobian is not being

		// evaluated, or if specifically requested
		updateTransport(x, j0, j1);
	}
	if (jg == npos) {
		double* Yleft = x + index(c_offset_Y, jmin);
		m_kExcessLeft = distance(Yleft, max_element(Yleft, Yleft + m_nsp));
		double* Yright = x + index(c_offset_Y, jmax);
		m_kExcessRight = distance(Yright, max_element(Yright, Yright + m_nsp));
	}

	// update the species diffusive mass fluxes whether or not a
	// Jacobian is being evaluated
	updateDiffFluxes(x, j0, j1);
}

void FlowBase1D::updateDiffFluxes(const doublereal* x, size_t j0, size_t j1)
{
	if (m_do_multicomponent) {
		for (size_t j = j0; j < j1; j++) {
			double dz = z(j + 1) - z(j);
			for (size_t k = 0; k < m_nsp; k++) {
				doublereal sum = 0.0;
				for (size_t m = 0; m < m_nsp; m++) {
					sum += m_wt[m] * m_multidiff[mindex(k, m, j)] * (X(x, m, j + 1) - X(x, m, j));
				}
				m_flux(k, j) = sum * m_diff[k + j * m_nsp] / dz;
			}
		}
	}
	else {
		for (size_t j = j0; j < j1; j++) {
			double sum = 0.0;
			double wtm = m_wtm[j];
			double rho = density(j);
			double dz = z(j + 1) - z(j);
			for (size_t k = 0; k < m_nsp; k++) {
				m_flux(k, j) = m_wt[k] * (rho*m_diff[k + m_nsp * j] / wtm);
				m_flux(k, j) *= (X(x, k, j) - X(x, k, j + 1)) / dz;
				sum -= m_flux(k, j);
			}
			// correction flux to insure that \sum_k Y_k V_k = 0.
			for (size_t k = 0; k < m_nsp; k++) {
				m_flux(k, j) += sum * Y(x, k, j);
			}
		}
	}

	if (m_do_soret) {
		for (size_t m = j0; m < j1; m++) {
			double gradlogT = 2.0 * (T(x, m + 1) - T(x, m)) /
				((T(x, m + 1) + T(x, m)) * (z(m + 1) - z(m)));
			for (size_t k = 0; k < m_nsp; k++) {
				m_flux(k, m) -= m_dthermal(k, m)*gradlogT;
			}
		}
	}
}

void FlowBase1D::setGas(const doublereal* x, size_t j)
{
    m_thermo->setTemperature(T(x,j));
    const doublereal* yy = x + m_nv*j + c_offset_Y;
    m_thermo->setMassFractions_NoNorm(yy);
    m_thermo->setPressure(m_press);
}

void FlowBase1D::setGasAtMidpoint(const doublereal* x, size_t j)
{
    m_thermo->setTemperature(0.5*(T(x,j)+T(x,j+1)));
    const doublereal* yyj = x + m_nv*j + c_offset_Y;
    const doublereal* yyjp = x + m_nv*(j+1) + c_offset_Y;
    for (size_t k = 0; k < m_nsp; k++) {
        m_ybar[k] = 0.5*(yyj[k] + yyjp[k]);
    }
    m_thermo->setMassFractions_NoNorm(m_ybar.data());
    m_thermo->setPressure(m_press);
}

void FlowBase1D::updateTransport(doublereal* x, size_t j0, size_t j1)
{
     if (m_do_multicomponent) {
        for (size_t j = j0; j < j1; j++) {
            setGasAtMidpoint(x,j);
            doublereal wtm = m_thermo->meanMolecularWeight();
            doublereal rho = m_thermo->density();
            m_visc[j] = (m_dovisc ? m_trans->viscosity() : 0.0);
            m_trans->getMultiDiffCoeffs(m_nsp, &m_multidiff[mindex(0,0,j)]);

            // Use m_diff as storage for the factor outside the summation
            for (size_t k = 0; k < m_nsp; k++) {
                m_diff[k+j*m_nsp] = m_wt[k] * rho / (wtm*wtm);
            }

            m_tcon[j] = m_trans->thermalConductivity();
            if (m_do_soret) {
                m_trans->getThermalDiffCoeffs(m_dthermal.ptrColumn(0) + j*m_nsp);
            }
        }
    } else { // mixture averaged transport
        for (size_t j = j0; j < j1; j++) {
            setGasAtMidpoint(x,j);
            m_visc[j] = (m_dovisc ? m_trans->viscosity() : 0.0);
            m_trans->getMixDiffCoeffs(&m_diff[j*m_nsp]);
            m_tcon[j] = m_trans->thermalConductivity();
        }
    }
}

void FlowBase1D::showSolution(const doublereal* x)
{
    writelog("    Pressure:  {:10.4g} Pa\n", m_press);

    Domain1D::showSolution(x);

    if (m_do_radiation) {
        writeline('-', 79, false, true);
        writelog("\n          z      radiative heat loss");
        writeline('-', 79, false, true);
        for (size_t j = 0; j < m_points; j++) {
            writelog("\n {:10.4g}        {:10.4g}", m_z[j], m_qdotRadiation[j]);
        }
        writelog("\n");
    }
}

void FlowBase1D::solveEnergyEqn(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (!m_do_energy[i]) {
                changed = true;
            }
            m_do_energy[i] = true;
        }
    } else {
        if (!m_do_energy[j]) {
            changed = true;
        }
        m_do_energy[j] = true;
    }
	if (changed) {
		needJacUpdate();
	}
}

void FlowBase1D::setBoundaryEmissivities(doublereal e_left, doublereal e_right)
{
    if (e_left < 0 || e_left > 1) {
        throw CanteraError("FlowBase1D::setBoundaryEmissivities",
            "The left boundary emissivity must be between 0.0 and 1.0!");
    } else if (e_right < 0 || e_right > 1) {
        throw CanteraError("FlowBase1D::setBoundaryEmissivities",
            "The right boundary emissivity must be between 0.0 and 1.0!");
    } else {
        m_epsilon_left = e_left;
        m_epsilon_right = e_right;
    }
}

void FlowBase1D::fixTemperature(size_t j)
{
    bool changed = false;
    if (j == npos) {
        for (size_t i = 0; i < m_points; i++) {
            if (m_do_energy[i]) {
                changed = true;
            }
            m_do_energy[i] = false;
        }
    } else {
        if (m_do_energy[j]) {
            changed = true;
        }
        m_do_energy[j] = false;
    }
	if (changed) {
		needJacUpdate();
	}
}

} // namespace
