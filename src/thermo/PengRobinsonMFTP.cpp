//! @file PengRobinsonMFTP.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/PengRobinsonMFTP.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"

#include <boost/math/tools/roots.hpp>

using namespace std;
namespace bmt = boost::math::tools;

namespace Cantera
{

const doublereal PengRobinsonMFTP::omega_a = 4.5723552892138218E-01;
const doublereal PengRobinsonMFTP::omega_b = 7.77960739038885E-02;
const doublereal PengRobinsonMFTP::omega_vc = 3.07401308698703833E-01;

PengRobinsonMFTP::PengRobinsonMFTP() :
	m_formTempParam(0),
	m_b_current(0.0),
	m_a_current(0.0),
	m_aAlpha_current(0.0),
    NSolns_(0),
    dpdV_(0.0),
    dpdT_(0.0)
{
	Vroot_[0] = 0.0;
	Vroot_[1] = 0.0;
	Vroot_[2] = 0.0;
}

PengRobinsonMFTP::PengRobinsonMFTP(const std::string& infile, const std::string& id_) :
    m_formTempParam(0),
    m_b_current(0.0),
    m_a_current(0.0),
	m_aAlpha_current(0.0),
    NSolns_(0),
    dpdV_(0.0),
    dpdT_(0.0)
{
    initThermoFile(infile, id_);
	Vroot_[0] = 0.0;
	Vroot_[1] = 0.0;
	Vroot_[2] = 0.0;
}

PengRobinsonMFTP::PengRobinsonMFTP(XML_Node& phaseRefRoot, const std::string& id_) :
    m_formTempParam(0),
    m_b_current(0.0),
	m_a_current(0.0), 
	m_aAlpha_current(0.0),
    NSolns_(0),
	dpdV_(0.0),
    dpdT_(0.0)
{
    importPhase(phaseRefRoot, this);
	Vroot_[0] = 0.0;
	Vroot_[1] = 0.0;
	Vroot_[2] = 0.0;
}

void PengRobinsonMFTP::calculateAlpha(const std::string& species, double a, double b, double w)
{
	size_t k = speciesIndex(species);
	if (k == npos) {
		throw CanteraError("PengRobinsonMFTP::setSpeciesCoeffs",
			"Unknown species '{}'.", species);
	}

	// Calculate value of kappa (independent of temperature)
	// w is an accentric factor of species and must be specified in the CTI file
 
	if (w <= 0.491) {
		kappa_vec_[k] = 0.37464 + 1.54226*w - 0.26992*w*w;
	}
	else {
		kappa_vec_[k] = 0.374642 + 1.487503*w - 0.164423*w*w + 0.016666*w*w*w;
	}

	//Calculate alpha (temperature dependent interaction parameter)
	double criTemp = speciesCritTemperature(a, b); // critical temperature of individual species
	double sqt_T_reduced = sqrt(temperature() / criTemp);
	double sqt_alpha = 1 + kappa_vec_[k] * (1 - sqt_T_reduced);
	alpha_vec_Curr_[k] = sqt_alpha*sqt_alpha;
}
	
void PengRobinsonMFTP::setSpeciesCoeffs(const std::string& species,
                                        double a, double b, double w)
{
	size_t k = speciesIndex(species);
	if (k == npos) {
		throw CanteraError("PengRobinsonMFTP::setSpeciesCoeffs",
			"Unknown species '{}'.", species);
	}	
	size_t counter = k + m_kk * k;
	a_coeff_vec(0, counter) = a;
	double aAlpha_k = a*alpha_vec_Curr_[k];
	aAlpha_coeff_vec(0, counter) = aAlpha_k;

	// standard mixing rule for cross-species interaction term
	for (size_t j = 0; j < m_kk; j++) {
		if (k == j) {
			continue;
		}
		double a0kj = sqrt(a_coeff_vec(0, j + m_kk * j) * a);
		double aAlpha_j = a*alpha_vec_Curr_[j];
		double a_Alpha = sqrt(aAlpha_j*aAlpha_k);
		if (a_coeff_vec(0, j + m_kk * k) == 0) {
			a_coeff_vec(0, j + m_kk * k) = a0kj;
			aAlpha_coeff_vec(0, j + m_kk * k) = a_Alpha;
			a_coeff_vec(0, k + m_kk * j) = a0kj;
			aAlpha_coeff_vec(0, k + m_kk * j) = a_Alpha;
		}
	}
	a_coeff_vec.getRow(0, a_vec_Curr_.data());
	aAlpha_coeff_vec.getRow(0, a_vec_Curr_.data());
	b_vec_Curr_[k] = b;
}

void PengRobinsonMFTP::setBinaryCoeffs(const std::string& species_i,
        const std::string& species_j, double a0, double alpha)
{
    size_t ki = speciesIndex(species_i);
    if (ki == npos) {
        throw CanteraError("PengRobinsonMFTP::setBinaryCoeffs",
            "Unknown species '{}'.", species_i);
    }
    size_t kj = speciesIndex(species_j);
    if (kj == npos) {
        throw CanteraError("PengRobinsonMFTP::setBinaryCoeffs",
            "Unknown species '{}'.", species_j);
    }

    size_t counter1 = ki + m_kk * kj;
    size_t counter2 = kj + m_kk * ki;
    a_coeff_vec(0, counter1) = a_coeff_vec(0, counter2) = a0;
	aAlpha_coeff_vec(0, counter1) = aAlpha_coeff_vec(0, counter2) = a0*alpha;
    a_vec_Curr_[counter1] = a_vec_Curr_[counter2] = a0;
	aAlpha_vec_Curr_[counter1] = aAlpha_vec_Curr_[counter2] = a0*alpha;
} 

// ------------Molar Thermodynamic Properties -------------------------

doublereal PengRobinsonMFTP::enthalpy_mole() const
{
    _updateReferenceStateThermo();
    doublereal h_ideal = RT() * mean_X(m_h0_RT);
    doublereal h_nonideal = hresid();
	return h_ideal + h_nonideal;
}

doublereal PengRobinsonMFTP::entropy_mole() const
{
    _updateReferenceStateThermo();
    doublereal sr_ideal = GasConstant * (mean_X(m_s0_R) - sum_xlogx() - std::log(pressure()/refPressure()));
    doublereal sr_nonideal = sresid();
	return sr_ideal + sr_nonideal;
}

doublereal PengRobinsonMFTP::cp_mole() const
{
    _updateReferenceStateThermo();
    doublereal TKelvin = temperature();
    doublereal sqt2 = sqrt(2);
    doublereal mv = molarVolume();
	doublereal vpb = mv + (1 + sqt2)*m_b_current;
	doublereal vmb = mv + (1 - sqt2)*m_b_current;
    pressureDerivatives();
    doublereal cpref = GasConstant * mean_X(m_cp0_R);
	doublereal dHdT_V = cpref + mv * dpdT_ - GasConstant + 1 / (2.0 * sqt2 *m_b_current) * log(vpb / vmb) * TKelvin *d2aAlpha_dT2();
	return dHdT_V - (mv + TKelvin * dpdT_ / dpdV_) * dpdT_;
}

doublereal PengRobinsonMFTP::cv_mole() const
{
    _updateReferenceStateThermo();
    doublereal TKelvin = temperature();
	doublereal mv = molarVolume();
	doublereal pp = pressure();
	pressureDerivatives(); 
	doublereal cp = cp_mole();
	return (cp_mole() + TKelvin* dpdT_* dpdT_ / dpdV_);
}

doublereal PengRobinsonMFTP::pressure() const
{
    _updateReferenceStateThermo();
    //  Get a copy of the private variables stored in the State object
    double T = temperature();
    double mv = meanMolecularWeight() / density();
	double den = mv * mv + 2 * mv * m_b_current - m_b_current * m_b_current;
	double pp = GasConstant * T / (mv - m_b_current) - m_aAlpha_current / den;
    return pp;
}

void PengRobinsonMFTP::calcDensity()
{
    // Calculate the molarVolume of the solution (m**3 kmol-1)
    const doublereal* const dtmp = moleFractdivMMW();
    getPartialMolarVolumes(m_tmpV.data());
    double invDens = dot(m_tmpV.begin(), m_tmpV.end(), dtmp);

    // Set the density in the parent State object directly, by calling the
    // Phase::setDensity() function.
    Phase::setDensity(1.0/invDens);
}

void PengRobinsonMFTP::setTemperature(const doublereal temp)
{
    Phase::setTemperature(temp);
    _updateReferenceStateThermo();
    updateAB();
}

void PengRobinsonMFTP::compositionChanged()
{
    MixtureFugacityTP::compositionChanged();
    updateAB();
}

void PengRobinsonMFTP::getActivityConcentrations(doublereal* c) const
{
    getActivityCoefficients(c);
    for (size_t k = 0; k < m_kk; k++) {
        c[k] *= moleFraction(k)*pressure()/RT();
    }
}

doublereal PengRobinsonMFTP::standardConcentration(size_t k) const
{
    getStandardVolumes(m_tmpV.data());
    return 1.0 / m_tmpV[k];
}

void PengRobinsonMFTP::getActivityCoefficients(doublereal* ac) const
{
    doublereal mv = molarVolume();
    doublereal T = temperature();
    doublereal sqt2 = sqrt(2);
    doublereal vpb_2 = mv + (1 + sqt2)*m_b_current;
    doublereal vmb_2 = mv + (1 - sqt2)*m_b_current;
	doublereal vmb = mv + (1 - sqt2)*m_b_current;
	doublereal pres = pressure();

    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            size_t counter = k + m_kk*i;
            m_pp[k] += moleFractions_[i] * aAlpha_vec_Curr_[counter];
        }
    }  
	doublereal num = 0;
	doublereal den = 2 * sqt2 * m_b_current * m_b_current;
	doublereal den2 = m_b_current*(mv * mv + 2 * mv * m_b_current - m_b_current * m_b_current);
    for (size_t k = 0; k < m_kk; k++) {
	    num = 2 * m_b_current * m_pp[k] - m_aAlpha_current* b_vec_Curr_[k];
        ac[k] = (-RT()*log(pres*mv/RT()) + RT() * log(mv / vmb)
                 + RT() * b_vec_Curr_[k] / vmb
                 - (num /den) * log(vpb_2/vmb_2)
                 - m_aAlpha_current* b_vec_Curr_[k] * mv/den2
                );
    }
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = exp(ac[k]/RT());
    }
}

// ---- Partial Molar Properties of the Solution -----------------

void PengRobinsonMFTP::getChemPotentials_RT(doublereal* muRT) const
{
    getChemPotentials(muRT);
    for (size_t k = 0; k < m_kk; k++) {
        muRT[k] *= 1.0 / RT();
    }
}

void PengRobinsonMFTP::getChemPotentials(doublereal* mu) const
{
    getGibbs_ref(mu);
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(SmallNumber, moleFraction(k));
        mu[k] += RT()*(log(xx));
    }

    doublereal mv = molarVolume();
    doublereal vmb = mv - m_b_current;
	doublereal sqt2 = sqrt(2);
	doublereal vpb_2 = mv + (1 + sqt2)*m_b_current;
	doublereal vmb_2 = mv + (1 - sqt2)*m_b_current;

    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            size_t counter = k + m_kk*i;
            m_pp[k] += moleFractions_[i] * aAlpha_vec_Curr_[counter];
        }
    }
    doublereal pres = pressure();
    doublereal refP = refPressure();
	doublereal num = 0;
	doublereal den = 2 * sqt2 * m_b_current * m_b_current;
	doublereal den2 = m_b_current*(mv * mv + 2 * mv * m_b_current - m_b_current * m_b_current);

    for (size_t k = 0; k < m_kk; k++) {
        num = 2 * m_b_current * m_pp[k] - m_aAlpha_current* b_vec_Curr_[k];
        
		mu[k] += (RT() * log(pres/refP) - RT() * log(pres * mv / RT())
                  + RT() * log(mv / vmb)
				  + RT() * b_vec_Curr_[k] / vmb
				  - (num /den) * log(vpb_2/vmb_2)
				  - m_aAlpha_current* b_vec_Curr_[k] * mv/den2
                 );
    }
}

void PengRobinsonMFTP::getPartialMolarEnthalpies(doublereal* hbar) const
{
    // First we get the reference state contributions
    getEnthalpy_RT_ref(hbar);
    scale(hbar, hbar+m_kk, hbar, RT());

    // We calculate dpdni_
    doublereal TKelvin = temperature();
    doublereal mv = molarVolume();
    doublereal sqt = sqrt(2);
    doublereal vmb = mv - m_b_current;
	doublereal sqt2 = sqrt(2);
	doublereal vpb_2 = mv + (1 + sqt2)*m_b_current;
	doublereal vmb_2 = mv + (1 - sqt2)*m_b_current;

    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            size_t counter = k + m_kk*i;
            m_pp[k] += moleFractions_[i] * aAlpha_vec_Curr_[counter];
        }
    }

	doublereal den = mv * mv + 2 * mv * m_b_current - m_b_current * m_b_current;
	doublereal den2 = den*den;
    for (size_t k = 0; k < m_kk; k++) {
	    dpdni_[k] = RT()/vmb + RT() * b_vec_Curr_[k] / (vmb * vmb) - 2.0 * m_pp[k] / den
			+ 2 * vmb * m_aAlpha_current * b_vec_Curr_[k] / den2;
    }

    doublereal daAlphadT = daAlpha_dT();
    doublereal fac = TKelvin * daAlphadT - m_aAlpha_current;

    pressureDerivatives();
    doublereal fac2 = mv + TKelvin * dpdT_ / dpdV_;
	doublereal fac3 = 2 * sqt * m_b_current *m_b_current;
    for (size_t k = 0; k < m_kk; k++) {
		double hE_v = mv * dpdni_[k] - RT() + (2 * m_b_current - b_vec_Curr_[k]) / fac3  * log(vpb_2 / vmb_2)*fac
					+ (mv * b_vec_Curr_[k]) /(m_b_current*den) * fac;
        hbar[k] = hbar[k] + hE_v;
        hbar[k] -= fac2 * dpdni_[k];
    }
}

void PengRobinsonMFTP::getPartialMolarEntropies(doublereal* sbar) const
{
    getEntropy_R_ref(sbar);
    scale(sbar, sbar+m_kk, sbar, GasConstant);
    doublereal TKelvin = temperature();
    doublereal mv = molarVolume();
    doublereal sqt2 = sqrt(2);
    doublereal vmb = mv - m_b_current;
	doublereal vpb_2 = mv + (1 + sqt2)*m_b_current;
	doublereal vmb_2 = mv + (1 - sqt2)*m_b_current;
	doublereal refP = refPressure();
	doublereal daAlphadT = daAlpha_dT();
	doublereal coeff1 = 0;
	doublereal den1 = 2 * sqt2 * m_b_current * m_b_current;
	doublereal den2 = mv * mv + 2 * mv * m_b_current - m_b_current * m_b_current;

    // Calculate sum(n_j (a alpha)_i,k * (1/alpha_k d/dT(alpha_k))) -> m_pp
	// Calculate sum(n_j (a alpha)_i,k * (1/alpha_i d/dT(alpha_i))) -> m_tmpV
	for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
		m_tmpV[k] = 0;
        for (size_t i = 0; i < m_kk; i++) {
            size_t counter = k + m_kk*i;
            m_pp[k] += moleFractions_[i] * aAlpha_vec_Curr_[counter];
			m_tmpV[k] += moleFractions_[i] * a_coeff_vec(1, counter) *(dalphadT_vec_Curr_[i] / alpha_vec_Curr_[i]);
        }
		m_pp[k] = m_pp[k] * dalphadT_vec_Curr_[k] / alpha_vec_Curr_[k];
    }
	

	for (size_t k = 0; k < m_kk; k++) {
		coeff1 = m_b_current * (m_pp[k] + m_tmpV[k]) - daAlphadT * b_vec_Curr_[k];
		sbar[k] += GasConstant * log(GasConstant * TKelvin / (refP * mv))
			+ GasConstant;
			+ GasConstant * log(mv / vmb)
			+ GasConstant * b_vec_Curr_[k] / vmb
			- coeff1* log(vpb_2 / vmb_2) / den1
			- b_vec_Curr_[k] * mv * daAlphadT / den2 / m_b_current;
	}
	pressureDerivatives();
    getPartialMolarVolumes(m_partialMolarVolumes.data());
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] -= m_partialMolarVolumes[k] * dpdT_;
    }
}

void PengRobinsonMFTP::getPartialMolarIntEnergies(doublereal* ubar) const
{
    getIntEnergy_RT(ubar);
    scale(ubar, ubar+m_kk, ubar, RT());
}

void PengRobinsonMFTP::getPartialMolarCp(doublereal* cpbar) const
{
    getCp_R(cpbar);
    scale(cpbar, cpbar+m_kk, cpbar, GasConstant);
}

void PengRobinsonMFTP::getPartialMolarVolumes(doublereal* vbar) const
{
    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = 0.0;
        for (size_t i = 0; i < m_kk; i++) {
            size_t counter = k + m_kk*i;
            m_pp[k] += moleFractions_[i] * aAlpha_vec_Curr_[counter];
        }
    }
    
	doublereal mv = molarVolume();
    doublereal vmb = mv - m_b_current;
    doublereal vpb = mv + m_b_current;
	doublereal fac = mv * mv + 2 * mv * m_b_current - m_b_current * m_b_current;
	doublereal fac2 = fac * fac;

    for (size_t k = 0; k < m_kk; k++) {
        doublereal num = (RT() + RT() * m_b_current/ vmb + RT() * b_vec_Curr_[k] / vmb
                          + RT() * m_b_current * b_vec_Curr_[k] /(vmb * vmb)
                          - 2 * mv * m_pp[k] / fac
                          + 2 * mv * vmb * m_aAlpha_current * b_vec_Curr_[k] / fac2
                         );
		doublereal denom = (pressure() + RT() * m_b_current / (vmb * vmb) 
							+ m_aAlpha_current/fac
							- 2 * mv* vpb *m_aAlpha_current / fac2
                           );
        vbar[k] = num / denom;
    }
}

doublereal PengRobinsonMFTP::speciesCritTemperature(double a, double b) const
{
	double pc, tc, vc;
	calcCriticalConditions(a, b, pc, tc, vc);
	return tc;
}

doublereal PengRobinsonMFTP::critTemperature() const
{
    double pc, tc, vc;
    calcCriticalConditions(m_a_current, m_b_current, pc, tc, vc);
    return tc;
}

doublereal PengRobinsonMFTP::critPressure() const
{
    double pc, tc, vc;
    calcCriticalConditions(m_a_current, m_b_current, pc, tc, vc);
    return pc;
}

doublereal PengRobinsonMFTP::critVolume() const
{
    double pc, tc, vc;
    calcCriticalConditions(m_a_current, m_b_current, pc, tc, vc);
    return vc;
}

doublereal PengRobinsonMFTP::critCompressibility() const
{
    double pc, tc, vc;
    calcCriticalConditions(m_a_current, m_b_current, pc, tc, vc);
    return pc*vc/tc/GasConstant;
}

doublereal PengRobinsonMFTP::critDensity() const
{
    double pc, tc, vc;
    calcCriticalConditions(m_a_current, m_b_current, pc, tc, vc);
    double mmw = meanMolecularWeight();
    return mmw / vc;
}

void PengRobinsonMFTP::setToEquilState(const doublereal* mu_RT)
{
    double tmp, tmp2;
    _updateReferenceStateThermo();
    getGibbs_RT_ref(m_tmpV.data());

    // Within the method, we protect against inf results if the exponent is too
    // high.
    //
    // If it is too low, we set the partial pressure to zero. This capability is
    // needed by the elemental potential method.
    doublereal pres = 0.0;
    double m_p0 = refPressure();
    for (size_t k = 0; k < m_kk; k++) {
        tmp = -m_tmpV[k] + mu_RT[k];
        if (tmp < -600.) {
            m_pp[k] = 0.0;
        } else if (tmp > 500.0) {
            tmp2 = tmp / 500.;
            tmp2 *= tmp2;
            m_pp[k] = m_p0 * exp(500.) * tmp2;
        } else {
            m_pp[k] = m_p0 * exp(tmp);
        }
        pres += m_pp[k];
    }
    // set state
    setState_PX(pres, &m_pp[0]);
}

bool PengRobinsonMFTP::addSpecies(shared_ptr<Species> spec)
{
    bool added = MixtureFugacityTP::addSpecies(spec);
    if (added) {
        a_vec_Curr_.resize(m_kk * m_kk, 0.0);
        b_vec_Curr_.push_back(0.0);
		a_vec_Curr_.push_back(0.0);
		aAlpha_vec_Curr_.resize(m_kk * m_kk, 0.0);
		aAlpha_vec_Curr_.push_back(0.0);
		kappa_vec_.push_back(0.0);

		alpha_vec_Curr_.push_back(0.0);
		a_coeff_vec.resize(1, m_kk * m_kk, 0.0);
		aAlpha_coeff_vec.resize(1, m_kk * m_kk, 0.0);
		dalphadT_vec_Curr_.push_back(0.0);
		d2alphadT2_.push_back(0.0);

        m_pp.push_back(0.0);
        m_tmpV.push_back(0.0);
        m_partialMolarVolumes.push_back(0.0);
        dpdni_.push_back(0.0);
    }
    return added;
}

vector<double> PengRobinsonMFTP::getCoeff(const std::string& iName)
{
	vector_fp spCoeff{ NAN, NAN };

	// Get number of species in the database
	// open xml file critProperties.xml
	XML_Node* doc = get_XML_File("critProperties.xml");
	size_t nDatabase = doc->nChildren();

	// Loop through all species in the database and attempt to match supplied
	// species to each. If present, calculate pureFluidParameters a_k and b_k
	// based on crit properties T_c and P_c:
	for (size_t isp = 0; isp < nDatabase; isp++) {
		XML_Node& acNodeDoc = doc->child(isp);
		std::string iNameLower = toLowerCopy(iName);
		std::string dbName = toLowerCopy(acNodeDoc.attrib("name"));

		// Attempt to match provided specie iName to current database species
		//  dbName:
		if (iNameLower == dbName) {
			// Read from database and calculate a and b coefficients
			double vParams;
			double T_crit, P_crit;

			if (acNodeDoc.hasChild("Tc")) {
				vParams = 0.0;
				XML_Node& xmlChildCoeff = acNodeDoc.child("Tc");
				if (xmlChildCoeff.hasAttrib("value"))
				{
					std::string critTemp = xmlChildCoeff.attrib("value");
					vParams = strSItoDbl(critTemp);
				}
				if (vParams <= 0.0) //Assuming that Pc and Tc are non zero.
				{
					throw CanteraError("PengRobinsonMFTP::GetCoeff",
						"Critical Temperature must be positive ");
				}
				T_crit = vParams;
			}
			if (acNodeDoc.hasChild("Pc")) {
				vParams = 0.0;
				XML_Node& xmlChildCoeff = acNodeDoc.child("Pc");
				if (xmlChildCoeff.hasAttrib("value"))
				{
					std::string critPressure = xmlChildCoeff.attrib("value");
					vParams = strSItoDbl(critPressure);
				}
				if (vParams <= 0.0) //Assuming that Pc and Tc are non zero.
				{
					throw CanteraError("PengRobinsonMFTP::GetCoeff",
						"Critical Pressure must be positive ");
				}
				P_crit = vParams;
			}

			//Assuming no temperature dependence
			spCoeff[0] = omega_a * pow(GasConstant, 2) * pow(T_crit, 2) / P_crit; //coeff a
			spCoeff[1] = omega_b * GasConstant * T_crit / P_crit; // coeff b
			break;
		}
	}
	return spCoeff;
}

void PengRobinsonMFTP::initThermoXML(XML_Node& phaseNode, const std::string& id)
{
	if (phaseNode.hasChild("thermo")) {
		XML_Node& thermoNode = phaseNode.child("thermo");
		std::string model = thermoNode["model"];
		if (model != "PengRobinson" && model != "PengRobinsonMFTP") {
			throw CanteraError("PengRobinsonMFTP::initThermoXML",
				"Unknown thermo model : " + model);
		}

		// Go get all of the coefficients and factors in the
		// activityCoefficients XML block
		if (thermoNode.hasChild("activityCoefficients")) {
			XML_Node& acNode = thermoNode.child("activityCoefficients");

			// Count the number of species with parameters provided in the
			//    input file:
			size_t nParams = 0;

			// Loop through the children and read out fluid parameters.  Process
			//   all the pureFluidParameters, first:
			for (size_t i = 0; i < acNode.nChildren(); i++) {
				XML_Node& xmlACChild = acNode.child(i);
				if (caseInsensitiveEquals(xmlACChild.name(), "purefluidparameters")) {
					readXMLPureFluid(xmlACChild);
					nParams += 1;
				}
			}

			// If any species exist which have undefined pureFluidParameters,
			// search the database in 'critProperties.xml' to find critical
			// temperature and pressure to calculate a and b.

			// Loop through all species in the CTI file
			size_t iSpecies = 0;

			for (size_t i = 0; i < m_kk; i++) {
				string iName = speciesName(i);

				// Get the index of the species
				iSpecies = speciesIndex(iName);

				// Check if a and b are already populated (only the diagonal elements of a).
				size_t counter = iSpecies + m_kk * iSpecies;

				// If not, then search the database:
				if (isnan(a_coeff_vec(0, counter))) {

					vector<double> coeffArray;

					// Search the database for the species name and calculate
					// coefficients a and b, from critical properties:
					// coeffArray[0] = a0, coeffArray[1] = b, coeffArray[2] = w;
					coeffArray = getCoeff(iName);

					// Check if species was found in the database of critical properties,
					// and assign the results
					if (!isnan(coeffArray[0])) {
						//Assuming no temperature dependence (i,e a1 = 0)
						setSpeciesCoeffs(iName, coeffArray[0], 0.0, coeffArray[1]);
					}
				}
			}

			// Loop back through the "activityCoefficients" children and process the
			// crossFluidParameters in the XML tree:
			for (size_t i = 0; i < acNode.nChildren(); i++) {
				XML_Node& xmlACChild = acNode.child(i);
				if (caseInsensitiveEquals(xmlACChild.name(), "crossfluidparameters")) {
					readXMLCrossFluid(xmlACChild);
				}
			}
		}
	}

	MixtureFugacityTP::initThermoXML(phaseNode, id);
}

void PengRobinsonMFTP::readXMLPureFluid(XML_Node& pureFluidParam)
{
	string xname = pureFluidParam.name();
	if (xname != "pureFluidParameters") {
		throw CanteraError("PengRobinsonMFTP::readXMLPureFluid",
			"Incorrect name for processing this routine: " + xname);
	}

	double a0 = 0.0;
	double a1 = 0.0;
	double b = 0.0;
	double w = 0.0;
	for (size_t iChild = 0; iChild < pureFluidParam.nChildren(); iChild++) {
		XML_Node& xmlChild = pureFluidParam.child(iChild);
		string nodeName = toLowerCopy(xmlChild.name());

		if (nodeName == "a_coeff") {
			vector_fp vParams;
			string iModel = toLowerCopy(xmlChild.attrib("model"));
			getFloatArray(xmlChild, vParams, true, "Pascal-m6/kmol2", "a_coeff");

			if (vParams.size() == 1) {
				a0 = vParams[0];
			}
			else if (vParams.size() == 2) {
				a0 = vParams[0];
				a1 = vParams[1];
			}
			else {
				throw CanteraError("PengRobinsonMFTP::readXMLPureFluid",
					"unknown model or incorrect number of parameters");
			}

		}
		else if (nodeName == "b_coeff") {
			b = getFloatCurrent(xmlChild, "toSI");
		}
		else if (nodeName == "accentricFactor") {
			w = getFloatCurrent(xmlChild);
		}
	}
	calculateAlpha(pureFluidParam.attrib("species"), a0, b, w);
	setSpeciesCoeffs(pureFluidParam.attrib("species"), a0, b, w);
}

void PengRobinsonMFTP::readXMLCrossFluid(XML_Node& CrossFluidParam)
{
	string xname = CrossFluidParam.name();
	if (xname != "crossFluidParameters") {
		throw CanteraError("PengRobinsonMFTP::readXMLCrossFluid",
			"Incorrect name for processing this routine: " + xname);
	}

	string iName = CrossFluidParam.attrib("species1");
	string jName = CrossFluidParam.attrib("species2");

	size_t num = CrossFluidParam.nChildren();
	for (size_t iChild = 0; iChild < num; iChild++) {
		XML_Node& xmlChild = CrossFluidParam.child(iChild);
		string nodeName = toLowerCopy(xmlChild.name());

		if (nodeName == "a_coeff") {
			vector_fp vParams;
			getFloatArray(xmlChild, vParams, true, "Pascal-m6/kmol2", "a_coeff");
			string iModel = toLowerCopy(xmlChild.attrib("model"));
			if (iModel == "constant" && vParams.size() == 1) {
				setBinaryCoeffs(iName, jName, vParams[0], 0.0);
			}
			else if (iModel == "linear_a") {
				setBinaryCoeffs(iName, jName, vParams[0], vParams[1]);
			}
			else {
				throw CanteraError("PengRobinsonMFTP::readXMLCrossFluid",
					"unknown model ({}) or wrong number of parameters ({})",
					iModel, vParams.size());
			}
		}
	}
}

void PengRobinsonMFTP::setParametersFromXML(const XML_Node& thermoNode)
{
    MixtureFugacityTP::setParametersFromXML(thermoNode);
    std::string model = thermoNode["model"];
}

doublereal PengRobinsonMFTP::sresid() const
{
    doublereal rho = density();
    doublereal mmw = meanMolecularWeight();
    doublereal molarV = mmw / rho;
    double hh = m_b_current / molarV;
    doublereal zz = z();
    doublereal alpha_1 = daAlpha_dT();
    doublereal T = temperature();
    doublereal sqt2 = sqrt(2.0);
	doublereal vpb = molarV + (1.0 + sqt2) *m_b_current;
	doublereal vmb = molarV + (1.0 - sqt2) *m_b_current;
	doublereal fac = alpha_1 / (2.0 * sqt2 * m_b_current);
	doublereal sresid_mol_R = log(zz*(1.0 - hh)) + fac * log(vpb / vmb) / GasConstant;
	return GasConstant * sresid_mol_R;
}

doublereal PengRobinsonMFTP::hresid() const
{
	doublereal rho = density();
    doublereal mmw = meanMolecularWeight();
    doublereal molarV = mmw / rho;
    doublereal zz = z();
	doublereal aAlpha_1 = daAlpha_dT();
	doublereal T = temperature();
	doublereal sqt2 = sqrt(2);
	doublereal vpb = molarV + (1 + sqt2) *m_b_current;
	doublereal vmb = molarV + (1 - sqt2) *m_b_current;
	doublereal fac = 1 / (2.0 * sqt2 * m_b_current);
	return GasConstant * T * (zz - 1.0) + fac * log(vpb / vmb) *(T * aAlpha_1 - m_aAlpha_current);
}

doublereal PengRobinsonMFTP::liquidVolEst(doublereal TKelvin, doublereal& presGuess) const
{
    double v = m_b_current * 1.1;
    double atmp;
    double btmp;
	double aAlphatmp;
	calculateAB(TKelvin, atmp, btmp, aAlphatmp);
    doublereal pres = std::max(psatEst(TKelvin), presGuess);
    double Vroot[3];
    bool foundLiq = false;
    int m = 0;
    while (m < 100 && !foundLiq) {
		int nsol = NicholsSolve(TKelvin, pres, atmp, btmp, aAlphatmp, Vroot);
        if (nsol == 1 || nsol == 2) {
            double pc = critPressure();
            if (pres > pc) {
                foundLiq = true;
            }
            pres *= 1.04;
        } else {
            foundLiq = true;
        }
    }

    if (foundLiq) {
        v = Vroot[0];
        presGuess = pres;
    } else {
        v = -1.0;
    }
    return v;
}

doublereal PengRobinsonMFTP::densityCalc(doublereal TKelvin, doublereal presPa, int phaseRequested, doublereal rhoguess)
{
    // It's necessary to set the temperature so that m_aAlpha_current is set correctly.
    setTemperature(TKelvin);
    double tcrit = critTemperature();
    doublereal mmw = meanMolecularWeight();
    if (rhoguess == -1.0) {
        if (phaseRequested != FLUID_GAS) {
            if (TKelvin > tcrit) {
                rhoguess = presPa * mmw / (GasConstant * TKelvin);
            } else {
                if (phaseRequested == FLUID_GAS || phaseRequested == FLUID_SUPERCRIT) {
                    rhoguess = presPa * mmw / (GasConstant * TKelvin);
                } else if (phaseRequested >= FLUID_LIQUID_0) {
                    double lqvol = liquidVolEst(TKelvin, presPa);
                    rhoguess = mmw / lqvol;
                }
            }
        } else {
            // Assume the Gas phase initial guess, if nothing is specified to
            // the routine
            rhoguess = presPa * mmw / (GasConstant * TKelvin);
        }
    }

    doublereal volguess = mmw / rhoguess;
	NSolns_ = NicholsSolve(TKelvin, presPa, m_a_current, m_b_current, m_aAlpha_current, Vroot_);

    doublereal molarVolLast = Vroot_[0];
    if (NSolns_ >= 2) {
        if (phaseRequested >= FLUID_LIQUID_0) {
            molarVolLast = Vroot_[0];
        } else if (phaseRequested == FLUID_GAS || phaseRequested == FLUID_SUPERCRIT) {
            molarVolLast = Vroot_[2];
        } else {
            if (volguess > Vroot_[1]) {
                molarVolLast = Vroot_[2];
            } else {
                molarVolLast = Vroot_[0];
            }
        }
    } else if (NSolns_ == 1) {
        if (phaseRequested == FLUID_GAS || phaseRequested == FLUID_SUPERCRIT || phaseRequested == FLUID_UNDEFINED) {
            molarVolLast = Vroot_[0];
        } else {
            return -2.0;
        }
    } else if (NSolns_ == -1) {
        if (phaseRequested >= FLUID_LIQUID_0 || phaseRequested == FLUID_UNDEFINED || phaseRequested == FLUID_SUPERCRIT) {
            molarVolLast = Vroot_[0];
        } else if (TKelvin > tcrit) {
            molarVolLast = Vroot_[0];
        } else {
            return -2.0;
        }
    } else {
        molarVolLast = Vroot_[0];
        return -1.0;
    }
    return mmw / molarVolLast;
}

doublereal PengRobinsonMFTP::densSpinodalLiquid() const
{
    double Vroot[3];
    double T = temperature();
	int nsol = NicholsSolve(T, pressure(), m_a_current, m_b_current, m_aAlpha_current, Vroot);
    if (nsol != 3) {
        return critDensity();
    }

    auto resid = [this, T](double v) {
        double pp;
        return dpdVCalc(T, v, pp);
    };

    boost::uintmax_t maxiter = 100;
    std::pair<double, double> vv = bmt::toms748_solve(
        resid, Vroot[0], Vroot[1], bmt::eps_tolerance<double>(48), maxiter);

    doublereal mmw = meanMolecularWeight();
    return mmw / (0.5 * (vv.first + vv.second));
}

doublereal PengRobinsonMFTP::densSpinodalGas() const
{
    double Vroot[3];
    double T = temperature();
	int nsol = NicholsSolve(T, pressure(), m_a_current, m_b_current, m_aAlpha_current, Vroot);
    if (nsol != 3) {
        return critDensity();
    }

    auto resid = [this, T](double v) {
        double pp;
        return dpdVCalc(T, v, pp);
    };

    boost::uintmax_t maxiter = 100;
    std::pair<double, double> vv = bmt::toms748_solve(
        resid, Vroot[1], Vroot[2], bmt::eps_tolerance<double>(48), maxiter);

    doublereal mmw = meanMolecularWeight();
    return mmw / (0.5 * (vv.first + vv.second));
}

doublereal PengRobinsonMFTP::pressureCalc(doublereal TKelvin, doublereal molarVol) const
{
	doublereal den = molarVol * molarVol + 2 * molarVol * m_b_current - m_b_current * m_b_current;
    double pres = GasConstant * TKelvin / (molarVol - m_b_current) - m_aAlpha_current / den;
    return pres;
}

doublereal PengRobinsonMFTP::dpdVCalc(doublereal TKelvin, doublereal molarVol, doublereal& presCalc) const
{
	doublereal den = molarVol * molarVol + 2 * molarVol * m_b_current - m_b_current * m_b_current; 
	presCalc = GasConstant * TKelvin / (molarVol - m_b_current) - m_aAlpha_current/ den;

    doublereal vpb = molarVol + m_b_current;
    doublereal vmb = molarVol - m_b_current;
	doublereal dpdv = -GasConstant * TKelvin / (vmb * vmb) + 2 *m_aAlpha_current * vpb / (den*den);
    return dpdv;
}

void PengRobinsonMFTP::pressureDerivatives() const
{
    doublereal TKelvin = temperature();
    doublereal mv = molarVolume();
    doublereal pres;

    dpdV_ = dpdVCalc(TKelvin, mv, pres);
    doublereal vmb = mv - m_b_current;
	doublereal den = mv * mv + 2 * mv * m_b_current - m_b_current * m_b_current;
	dpdT_ = (GasConstant / vmb - daAlpha_dT() / den);
}

void PengRobinsonMFTP::updateMixingExpressions()
{
    updateAB();
}

void PengRobinsonMFTP::updateAB()
{
    double temp = temperature();
	//Update aAlpha_i
	double sqt_alpha;
	double criTemp = critTemperature();
	double sqt_T_reduced = sqrt(temp / criTemp);

	// Update indiviual alpha
	for (size_t j = 0; j < m_kk; j++) {
		sqt_alpha = 1 + kappa_vec_[j] * (1 - sqt_T_reduced);
		alpha_vec_Curr_[j] = sqt_alpha*sqt_alpha;
	}
	
	//Update aAlpha_i,j
	for (size_t i = 0; i < m_kk; i++) {
		for (size_t j = 0; j < m_kk; j++) {
			size_t counter = i * m_kk + j;
			a_vec_Curr_[counter] = a_coeff_vec(0, counter);
			aAlpha_vec_Curr_[counter] = sqrt(alpha_vec_Curr_[i] * alpha_vec_Curr_[j]) * a_coeff_vec(0, counter);
		}
	}	

    m_b_current = 0.0;
    m_a_current = 0.0;
	m_aAlpha_current = 0.0;

    for (size_t i = 0; i < m_kk; i++) {
		m_b_current += moleFractions_[i] * b_vec_Curr_[i];
        for (size_t j = 0; j < m_kk; j++) {
            m_a_current += a_vec_Curr_[i * m_kk + j] * moleFractions_[i] * moleFractions_[j];
			m_aAlpha_current += aAlpha_vec_Curr_[i * m_kk + j] * moleFractions_[i] * moleFractions_[j];
        }
    }
}

void PengRobinsonMFTP::calculateAB(doublereal temp, doublereal& aCalc, doublereal& bCalc, doublereal& aAlphaCalc) const
{
	bCalc = 0.0;
	aCalc = 0.0;
	aAlphaCalc = 0.0;
	for (size_t i = 0; i < m_kk; i++) {
		bCalc += moleFractions_[i] * b_vec_Curr_[i];
		for (size_t j = 0; j < m_kk; j++) {
			size_t counter = i * m_kk + j;
			doublereal a_vec_Curr = a_coeff_vec(0, counter);
			aCalc += a_vec_Curr * moleFractions_[i] * moleFractions_[j];
			aAlphaCalc += aAlpha_vec_Curr_[counter] * moleFractions_[i] * moleFractions_[j];
		}
	}
}

doublereal PengRobinsonMFTP::daAlpha_dT() const
{
	doublereal daAlphadT = 0.0, temp, k, Tc = 0.0, sqtTr = 0.0;
	//we need species critical temperature
	double coeff1 = 1 / (critTemperature()*sqtTr);
	double coeff2 = sqtTr - 1;
	for (size_t i = 0; i < m_kk; i++) {
		size_t counter = i + m_kk * i;
		// Calculate first and double derivatives of alpha for individual species
		Tc = speciesCritTemperature(a_vec_Curr_[counter], b_vec_Curr_[i]);
		sqtTr = sqrt(temperature() / Tc); //we need species critical temperature
		coeff1 = 1 / (Tc*sqtTr);
		coeff2 = sqtTr - 1;
		k = kappa_vec_[i];
		dalphadT_vec_Curr_[i] = coeff1 *(k* k*coeff2 - k);
	}
	//Calculate mixture derivative
	for (size_t i = 0; i < m_kk; i++) {
		size_t counter1 = i + m_kk * i;
		for (size_t j = 0; j < m_kk; j++) {
			size_t counter2 = j * m_kk + j;
			temp = 0.5 * sqrt((a_vec_Curr_[counter1] * a_vec_Curr_[counter2]) / (alpha_vec_Curr_[i] * alpha_vec_Curr_[j]));
			daAlphadT += moleFractions_[i] * moleFractions_[j] * temp
				* (dalphadT_vec_Curr_[j] * alpha_vec_Curr_[i] + dalphadT_vec_Curr_[i] * alpha_vec_Curr_[j]);	
		}
	}
	return daAlphadT;
}

doublereal PengRobinsonMFTP::d2aAlpha_dT2() const
{
	doublereal daAlphadT = 0.0, temp, fac1, fac2, alphaij, alphai, alphaj, d2aAlphadT2 = 0.0, num;
	double k;
	double sqt_Tr = sqrt(temperature() / critTemperature()); //we need species critical temperature
	double coeff1 = 1 / (critTemperature()*critTemperature()*sqt_Tr);
	double coeff2 = sqt_Tr - 1;
	for (size_t i = 0; i < m_kk; i++) {
		// Calculate individual dAlpha_dTi
		size_t counter = i + m_kk * i;
		k = kappa_vec_[i];
		dalphadT_vec_Curr_[i] = coeff1 *(k* k*coeff2 - k);
		d2alphadT2_[i] = (k*k + k) * coeff1 / (2 * sqt_Tr*sqt_Tr);
	}

	//Calculate mixture derivative
	for (size_t i = 0; i < m_kk; i++) {
		size_t counter1 = i + m_kk * i;
		alphai = alpha_vec_Curr_[i];
		for (size_t j = 0; j < m_kk; j++) {
			size_t counter2 = j + m_kk * j;
			alphaj = alpha_vec_Curr_[j];
			alphaij = alphai * alphaj;
			temp = 0.5 * sqrt((a_vec_Curr_[counter1] * a_vec_Curr_[counter2]) / (alphaij));
			num = (dalphadT_vec_Curr_[j] * alphai + dalphadT_vec_Curr_[i] * alphaj);
			fac1 = -(0.5 / alphaij)*num*num;
			fac2 = alphaj * d2alphadT2_[counter1] + alphai *d2alphadT2_[counter2] + 2 * dalphadT_vec_Curr_[i] * dalphadT_vec_Curr_[j];
			d2aAlphadT2 += moleFractions_[i] * moleFractions_[j] * temp *(fac1 + fac2);
		}
	}
	return d2aAlphadT2;
}

void PengRobinsonMFTP::calcCriticalConditions(doublereal a, doublereal b,
        doublereal& pc, doublereal& tc, doublereal& vc) const
{
    if (b <= 0.0) {
        tc = 1000000.;
        pc = 1.0E13;
        vc = omega_vc * GasConstant * tc / pc;
        return;
    }
    if (a <= 0.0) {
        tc = 0.0;
        pc = 0.0;
        vc = 2.0 * b;
        return;
    }    
    tc = a * omega_b / (b * omega_a * GasConstant);
    pc = omega_b * GasConstant * tc / b;
    vc = omega_vc * GasConstant * tc / pc;
}

int PengRobinsonMFTP::NicholsSolve(double TKelvin, double pres, doublereal a, doublereal b, doublereal aAlpha,
                                   doublereal Vroot[3]) const
{
	double tmp;
	Vroot[0] = 0.0;
    Vroot[1] = 0.0;
    Vroot[2] = 0.0;
    if (TKelvin <= 0.0) {
        throw CanteraError("PengRobinsonMFTP::NicholsSolve()", "neg temperature");
    }

    // Derive the coefficients of the cubic polynomial (in terms of molar volume v) to solve.
	doublereal bsqr = b * b;
	doublereal RT_p = GasConstant * TKelvin / pres;
	doublereal aAlpha_p = aAlpha / pres;
    doublereal an = 1.0;
	doublereal bn = (b - RT_p);
	doublereal cn = -(2 * RT_p * b - aAlpha_p + 3 * bsqr);
	doublereal dn = (bsqr * RT_p + bsqr * b - aAlpha_p * b);

    double tc = a * omega_b / (b * omega_a * GasConstant);
    double pc = omega_b * GasConstant * tc / b;
    double vc = omega_vc * GasConstant * tc / pc;
    
	// Derive the center of the cubic, x_N
    doublereal xN = - bn /(3 * an);

    // Derive the value of delta**2. This is a key quantity that determines the
    // number of turning points
    doublereal delta2 = (bn * bn - 3 * an * cn) / (9 * an * an); //This is delta^2 term in Nichols method
    doublereal delta = 0.0;

    // Calculate a couple of ratios
	// Cubic equation in z : z^2 - (1-B) z^2 + (A-2B -3B^2)z - (AB-B^2-B^3) = 0
    doublereal ratio1 = 3.0 * an * cn / (bn * bn);
    doublereal ratio2 = pres * b / (GasConstant * TKelvin); // B
    if (fabs(ratio1) < 1.0E-7) {
		doublereal ratio3 = aAlpha / (GasConstant * TKelvin) * pres / (GasConstant * TKelvin); // A
        if (fabs(ratio2) < 1.0E-5 && fabs(ratio3) < 1.0E-5) {
            // A and B terms in cubic equation for z are almost zero, then z is mear to 1
			doublereal zz = 1.0;
            for (int i = 0; i < 10; i++) {
                doublereal znew = zz / (zz - ratio2) - ratio3 / (zz + ratio1);
                doublereal deltaz = znew - zz;
                zz = znew;
                if (fabs(deltaz) < 1.0E-14) {
                    break;
                }
            }
            doublereal v = zz * GasConstant * TKelvin / pres;
            Vroot[0] = v;
            return 1;
        }
    }

    int nSolnValues;
    double h2 = 4. * an * an * delta2 * delta2 * delta2; // h^2
    if (delta2 > 0.0) {
        delta = sqrt(delta2);
    }

    doublereal h = 2.0 * an * delta * delta2;
    doublereal yN = 2.0 * bn * bn * bn / (27.0 * an * an) - bn * cn / (3.0 * an) + dn; // y_N term
    doublereal desc = yN * yN - h2; // descriminant

	//check if y = h
    if (fabs(fabs(h) - fabs(yN)) < 1.0E-10) {
        if (desc != 0.0) {
            // this is for getting to other cases
            throw CanteraError("NicholsSolve()", "numerical issues");
        }
        desc = 0.0;
    }
		
    if (desc < 0.0) {
		// desc<0 then we have three distinct roots.
        nSolnValues = 3;
    } else if (desc == 0.0) {
		// desc=0 then we have two distinct roots (third one is repeated root)
        nSolnValues = 2;
        // We are here as p goes to zero.
    } else if (desc > 0.0) {
		// desc> 0 then we have one real root.
        nSolnValues = 1;
    }

    // One real root -> have to determine whether gas or liquid is the root
    if (desc > 0.0) {
        doublereal tmpD = sqrt(desc);
        doublereal tmp1 = (- yN + tmpD) / (2.0 * an);
        doublereal sgn1 = 1.0;
        if (tmp1 < 0.0) {
            sgn1 = -1.0;
            tmp1 = -tmp1;
        }
        doublereal tmp2 = (- yN - tmpD) / (2.0 * an);
        doublereal sgn2 = 1.0;
        if (tmp2 < 0.0) {
            sgn2 = -1.0;
            tmp2 = -tmp2;
        }
        doublereal p1 = pow(tmp1, 1./3.);
        doublereal p2 = pow(tmp2, 1./3.);
        doublereal alpha = xN + sgn1 * p1 + sgn2 * p2;
        Vroot[0] = alpha;
        Vroot[1] = 0.0;
        Vroot[2] = 0.0;
        tmp = an * Vroot[0] * Vroot[0] * Vroot[0] + bn * Vroot[0] * Vroot[0] + cn * Vroot[0] + dn;
    } else if (desc < 0.0) {
        doublereal tmp = - yN/h;
        doublereal val = acos(tmp);
        doublereal theta = val / 3.0;
        doublereal oo = 2. * Pi / 3.;
        doublereal alpha = xN + 2. * delta * cos(theta);
        doublereal beta = xN + 2. * delta * cos(theta + oo);
        doublereal gamma = xN + 2. * delta * cos(theta + 2.0 * oo);
        Vroot[0] = beta;
        Vroot[1] = gamma;
        Vroot[2] = alpha;

        for (int i = 0; i < 3; i++) {
            tmp = an * Vroot[i] * Vroot[i] * Vroot[i] + bn * Vroot[i] * Vroot[i] + cn * Vroot[i] + dn;
            if (fabs(tmp) > 1.0E-4) {
                for (int j = 0; j < 3; j++) {
                    if (j != i && fabs(Vroot[i] - Vroot[j]) < 1.0E-4 * (fabs(Vroot[i]) + fabs(Vroot[j]))) {
                        writelog("PengRobinsonMFTP::NicholsSolve(T = {}, p = {}):"
                                 " WARNING roots have merged: {}, {}\n",
                                 TKelvin, pres, Vroot[i], Vroot[j]);
                    }
                }
            }
        }
    } else if (desc == 0.0) {
        if (yN == 0.0 && h == 0.0) {
            Vroot[0] = xN;
            Vroot[1] = xN;
            Vroot[2] = xN;
        } else {
            // need to figure out whether delta is pos or neg
            if (yN > 0.0) {
                tmp = pow(yN/(2*an), 1./3.);
                if (fabs(tmp - delta) > 1.0E-9) {
                    throw CanteraError("PengRobinsonMFTP::NicholsSolve()", "unexpected");
                }
                Vroot[1] = xN + delta;
                Vroot[0] = xN - 2.0*delta; // liquid phase root
            } else {
                tmp = pow(yN/(2*an), 1./3.);
                if (fabs(tmp - delta) > 1.0E-9) {
                    throw CanteraError("PengRobinsonMFTP::NicholsSolve()", "unexpected");
                }
                delta = -delta;
                Vroot[0] = xN + delta;
                Vroot[1] = xN - 2.0*delta; // gas phase root
            }
        }
        for (int i = 0; i < 2; i++) {
            tmp = an * Vroot[i] * Vroot[i] * Vroot[i] + bn * Vroot[i] * Vroot[i] + cn * Vroot[i] + dn;
        }
    }

    // Unfortunately, there is a heavy amount of roundoff error due to bad
    // conditioning in this
    double res, dresdV = 0.0;
    for (int i = 0; i < nSolnValues; i++) {
        for (int n = 0; n < 20; n++) {
            res = an * Vroot[i] * Vroot[i] * Vroot[i] + bn * Vroot[i] * Vroot[i] + cn * Vroot[i] + dn;
            if (fabs(res) < 1.0E-14) {
                break;
            }
            dresdV = 3.0 * an * Vroot[i] * Vroot[i] + 2.0 * bn * Vroot[i] + cn;
            double del = - res / dresdV;
            Vroot[i] += del;
            if (fabs(del) / (fabs(Vroot[i]) + fabs(del)) < 1.0E-14) {
                break;
            }
            double res2 = an * Vroot[i] * Vroot[i] * Vroot[i] + bn * Vroot[i] * Vroot[i] + cn * Vroot[i] + dn;
            if (fabs(res2) < fabs(res)) {
                continue;
            } else {
                Vroot[i] -= del;
                Vroot[i] += 0.1 * del;
            }
        }
        if ((fabs(res) > 1.0E-14) && (fabs(res) > 1.0E-14 * fabs(dresdV) * fabs(Vroot[i]))) {
            writelog("PengRobinsonMFTP::NicholsSolve(T = {}, p = {}): "
                "WARNING root didn't converge V = {}", TKelvin, pres, Vroot[i]);
            writelogendl();
        }
    }

    if (nSolnValues == 1) {
        if (TKelvin > tc) {
            if (Vroot[0] < vc) {
                nSolnValues = -1;
            }
        } else {
            if (Vroot[0] < xN) {
                nSolnValues = -1;
            }
        }
    } else {
        if (nSolnValues == 2 && delta > 0.0) {
            nSolnValues = -2;
        }
    }
    return nSolnValues;
}

}
