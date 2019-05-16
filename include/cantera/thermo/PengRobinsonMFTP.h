//! @file PengRobinsonMFTP.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_PENGROBINSONMFTP_H
#define CT_PENGROBINSONMFTP_H

#include "MixtureFugacityTP.h"
#include "cantera/base/Array.h"

namespace Cantera
{
/**
 * Implementation of a multi-species Peng-Robinson equation of state
 *
 * @ingroup thermoprops
 */
class PengRobinsonMFTP : public MixtureFugacityTP
{
public:
    //! @name Constructors and Duplicators
    //! @{

    //! Base constructor.
    PengRobinsonMFTP();

    //! Construct and initialize a PengRobinsonMFTP object directly from an
    //! ASCII input file
    /*!
     * @param infile    Name of the input file containing the phase XML data
     *                  to set up the object
     * @param id        ID of the phase in the input file. Defaults to the empty
     *     string.
     */
    PengRobinsonMFTP(const std::string& infile, const std::string& id="");

    //! Construct and initialize a PengRobinsonMFTP object directly from an
    //! XML database
    /*!
     *  @param phaseRef XML phase node containing the description of the phase
     *  @param id       id attribute containing the name of the phase.  (default
     *      is the empty string)
     */
    PengRobinsonMFTP(XML_Node& phaseRef, const std::string& id = "");

    virtual std::string type() const {
        return "PengRobinson";
    }

    //! @name Molar Thermodynamic properties
    //! @{

    virtual doublereal enthalpy_mole() const;
    virtual doublereal entropy_mole() const;
    virtual doublereal cp_mole() const;
    virtual doublereal cv_mole() const;

    //! @}
    //! @name Mechanical Properties
    //! @{

    //! Return the thermodynamic pressure (Pa).
    /*!
     *  Since the mass density, temperature, and mass fractions are stored,
     *  this method uses these values to implement the
     *  mechanical equation of state \f$ P(T, \rho, Y_1, \dots, Y_K) \f$.
     *
     * \f[
     *    P = \frac{RT}{v-b_{mix}} - \frac{\left(\alpha a\right)_{mix}}{v^2 + 2b_{mix}v - b_{mix}^2 }
     * \f]
     *
     *  where:
     *
     * \f[
     *    \alpha = \left[ 1 + \left(0.37464 + 1.54226\omega - 0.26992\omega^2\right)\left(1-T_r^{0.5}\right)\right]^2
     * \f]
     */
    virtual doublereal pressure() const;

    // @}

protected:
    /**
     * Calculate the density of the mixture using the partial molar volumes and
     * mole fractions as input
     *
     * The formula for this is
     *
     * \f[
     * \rho = \frac{\sum_k{X_k W_k}}{\sum_k{X_k V_k}}
     * \f]
     *
     * where \f$X_k\f$ are the mole fractions, \f$W_k\f$ are the molecular
     * weights, and \f$V_k\f$ are the pure species molar volumes.
     *
     * Note, the basis behind this formula is that in an ideal solution the
     * partial molar volumes are equal to the species standard state molar
     * volumes. The species molar volumes may be functions of temperature and
     * pressure.
     */
    virtual void calcDensity();

    virtual void setTemperature(const doublereal temp);
    virtual void compositionChanged();

public:
    virtual void getActivityConcentrations(doublereal* c) const;

    //! Returns the standard concentration \f$ C^0_k \f$, which is used to
    //! normalize the generalized concentration.
    /*!
     * This is defined as the concentration by which the generalized
     * concentration is normalized to produce the activity. In many cases, this
     * quantity will be the same for all species in a phase. Since the activity
     * for an ideal gas mixture is simply the mole fraction, for an ideal gas
     * \f$ C^0_k = P/\hat R T \f$.
     *
     * @param k Optional parameter indicating the species. The default is to
     *          assume this refers to species 0.
     * @return
     *   Returns the standard Concentration in units of m3 kmol-1.
     */
    virtual doublereal standardConcentration(size_t k=0) const;

    //! Get the array of non-dimensional activity coefficients at the current
    //! solution temperature, pressure, and solution concentration.
    /*!
     * For all objects with the Mixture Fugacity approximation, we define the
     * standard state as an ideal gas at the current temperature and pressure of
     * the solution. The activities are based on this standard state.
     *
     * @param ac Output vector of activity coefficients. Length: m_kk.
     */
    virtual void getActivityCoefficients(doublereal* ac) const;

    /// @name  Partial Molar Properties of the Solution
    //@{

    //! Get the array of non-dimensional species chemical potentials.
    //! These are partial molar Gibbs free energies.
    /*!
     * \f$ \mu_k / \hat R T \f$.
     * Units: unitless
     *
     * We close the loop on this function, here, calling getChemPotentials() and
     * then dividing by RT. No need for child classes to handle.
     *
     * @param mu    Output vector of non-dimensional species chemical potentials
     *              Length: m_kk.
     */
    virtual void getChemPotentials_RT(doublereal* mu) const;

    virtual void getChemPotentials(doublereal* mu) const;
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;
    virtual void getPartialMolarEntropies(doublereal* sbar) const;
    virtual void getPartialMolarIntEnergies(doublereal* ubar) const;
    virtual void getPartialMolarCp(doublereal* cpbar) const;
    virtual void getPartialMolarVolumes(doublereal* vbar) const;

    
	virtual void calculateAlpha(const std::string& species, double a, double b, double w);
	//@}
    /// @name Critical State Properties.
    //@{

    virtual doublereal critTemperature() const;
    virtual doublereal critPressure() const;
    virtual doublereal critVolume() const;
    virtual doublereal critCompressibility() const;
    virtual doublereal critDensity() const;
	virtual doublereal speciesCritTemperature(double a, double b) const;

public:
    //@}
    //! @name Initialization Methods - For Internal use
    /*!
     * The following methods are used in the process of constructing
     * the phase and setting its parameters from a specification in an
     * input file. They are not normally used in application programs.
     * To see how they are used, see importPhase().
     */
    //@{

    virtual bool addSpecies(shared_ptr<Species> spec);
    virtual void setParametersFromXML(const XML_Node& thermoNode);
    virtual void setToEquilState(const doublereal* lambda_RT);
    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);

	//! Retrieve a and b coefficients by looking up tabulated critical parameters
	/*!
	*  If pureFluidParameters are not provided for any species in the phase,
	*  consult the critical properties tabulated in /thermo/critProperties.xml.
	*  If the species is found there, calculate pure fluid parameters a_k and b_k as:
	*  \f[ a_k = 0.4278*R**2*T_c^2/P_c \f]
	*
	*  and:
	*  \f[ b_k = 0.08664*R*T_c/P_c \f]
	*
	*  @param iName    Name of the species
	*/
	virtual std::vector<double> getCoeff(const std::string& iName);

    //! Set the pure fluid interaction parameters for a species
    /*!
     *  The "a" parameter for species *i* in the Peng-Robinson model is assumed
     *  to be a linear function of temperature:
     *  \f[ a = a_0 + a_1 T \f]
     *
     *  @param species   Name of the species
     *  @param a0        constant term in the expression for the "a" parameter
     *      of the specified species [Pa-m^6/kmol^2]
     *  @param a1        temperature-proportional term in the expression for the
     *      "a" parameter of the specified species [Pa-m^6/kmol^2/K]
     *  @param b         "b" parameter in the Peng-Robinson model [m^3/kmol]
	 *  @param alpha	 dimensionless function of T_r and \omega
	 *  @param omega	 acentric factor
     */
    void setSpeciesCoeffs(const std::string& species, double a, double b,
                              double w);

    //! Set values for the interaction parameter between two species
    /*!
     *  The "a" parameter for interactions between species *i* and *j* is
     *  assumed by default to be computed as:
     *  \f[ a_{ij} = \sqrt(a_{i,0} a_{j,0}) + \sqrt(a_{i,1} a_{j,1}) T \f]
     *
     *  This function overrides the defaults with the specified parameters:
     *  \f[ a_{ij} = a_{ij,0} + a_{ij,1} T \f]
     *
     *  @param species_i   Name of one species
     *  @param species_j   Name of the other species
     *  @param a0          constant term in the "a" expression [Pa-m^6/kmol^2]
     *  @param a1          temperature-proportional term in the "a" expression
     *      [Pa-m^6/kmol^2/K]
     */
    void setBinaryCoeffs(const std::string& species_i,
                         const std::string& species_j, double a0, double a1);

private:
    //! Read the pure species PengRobinson input parameters
    /*!
     *  @param pureFluidParam   XML_Node for the pure fluid parameters
     */
    void readXMLPureFluid(XML_Node& pureFluidParam);

    //! Read the cross species PengRobinson input parameters
    /*!
     *  @param pureFluidParam   XML_Node for the cross fluid parameters
     */
    void readXMLCrossFluid(XML_Node& pureFluidParam);

    // @}

protected:
    // Special functions inherited from MixtureFugacityTP
    virtual doublereal sresid() const;
    virtual doublereal hresid() const;

public:
    virtual doublereal liquidVolEst(doublereal TKelvin, doublereal& pres) const;
    virtual doublereal densityCalc(doublereal TKelvin, doublereal pressure, int phase, doublereal rhoguess);

    virtual doublereal densSpinodalLiquid() const;
    virtual doublereal densSpinodalGas() const;
    virtual doublereal pressureCalc(doublereal TKelvin, doublereal molarVol) const;
    virtual doublereal dpdVCalc(doublereal TKelvin, doublereal molarVol, doublereal& presCalc) const;

    //! Calculate dpdV and dpdT at the current conditions
    /*!
     *  These are stored internally.
     */
    void pressureDerivatives() const;

    virtual void updateMixingExpressions();

    //! Update the a and b parameters
    /*!
     *  The a and the b parameters depend on the mole fraction and the
     *  temperature. This function updates the internal numbers based on the
     *  state of the object.
     */
    void updateAB();

    //! Calculate the a and the b parameters given the temperature
    /*!
     * This function doesn't change the internal state of the object, so it is a
     * const function.  It does use the stored mole fractions in the object.
     *
     * @param temp  Temperature (TKelvin)
     * @param aCalc (output)  Returns the a value
     * @param bCalc (output)  Returns the b value.
     */
	void calculateAB(doublereal temp, doublereal& aCalc, doublereal& bCalc, doublereal& aAlpha) const;

    // Special functions not inherited from MixtureFugacityTP

	doublereal daAlpha_dT() const;
	doublereal d2aAlpha_dT2() const;

    void calcCriticalConditions(doublereal a, doublereal b,doublereal& pc, doublereal& tc, doublereal& vc) const;

    //! Solve the cubic equation of state
    /*!
     * The P-R equation of state may be solved via the following formula:
     *
     *     V**3 - V**2(RT/P - b)  - V(2bRT/P - \alpha a/P + 3*b*b) - (a \alpha b/p - b*b RT/P - b*b*b) = 0
     *
     * Returns the number of solutions found. If it only finds the liquid
     * branch solution, it will return a -1 or a -2 instead of 1 or 2.  If it
     * returns 0, then there is an error.
     */
	int NicholsSolve(double TKelvin, double pres, doublereal a, doublereal b, doublereal aAlpha,
                     doublereal Vroot[3]) const;

protected:
    //! Form of the temperature parameterization
    /*!
     *  - 0 = There is no temperature parameterization of a or b
     *  - 1 = The a_ij parameter is a linear function of the temperature
     */
    int m_formTempParam;

    //! Value of b in the equation of state
    /*!
     *  m_b is a function of the temperature and the mole fraction.
     */
    doublereal m_b_current;

    //! Value of a in the equation of state
    /*!
     *  a_b is a function of the temperature and the mole fraction.
     */
    doublereal m_a_current;
	doublereal m_aAlpha_current;

    vector_fp a_vec_Curr_;
    vector_fp b_vec_Curr_;
	vector_fp aAlpha_vec_Curr_;
	vector_fp alpha_vec_Curr_;
	vector_fp kappa_vec_;
	mutable vector_fp dalphadT_vec_Curr_;
	mutable vector_fp d2alphadT2_;

    Array2D a_coeff_vec;
	Array2D aAlpha_coeff_vec;
	
	int NSolns_;

    doublereal Vroot_[3];

    //! Temporary storage - length = m_kk.
    mutable vector_fp m_pp;

    //! Temporary storage - length = m_kk.
    mutable vector_fp m_tmpV;

    // Partial molar volumes of the species
    mutable vector_fp m_partialMolarVolumes;

    //! The derivative of the pressure wrt the volume
    /*!
     * Calculated at the current conditions. temperature and mole number kept
     * constant
     */
    mutable doublereal dpdV_;

    //! The derivative of the pressure wrt the temperature
    /*!
     *  Calculated at the current conditions. Total volume and mole number kept
     *  constant
     */
    mutable doublereal dpdT_;

    //! Vector of derivatives of pressure wrt mole number
    /*!
     *  Calculated at the current conditions. Total volume, temperature and
     *  other mole number kept constant
     */
    mutable vector_fp dpdni_;

public:
    //! Omega constant for a -> value of a in terms of critical properties
    /*!
     *  this was calculated from a small nonlinear solve
     */
    static const doublereal omega_a;

    //! Omega constant for b
    static const doublereal omega_b;

    //! Omega constant for the critical molar volume
    static const doublereal omega_vc;
};
}

#endif
