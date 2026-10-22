/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpinFlipHydrogenGasMix.hpp"
#include "Configuration.hpp"
#include "Constants.hpp"
#include "DisjointWavelengthGrid.hpp"
#include "FatalError.hpp"
#include "MaterialState.hpp"

////////////////////////////////////////////////////////////////////

namespace
{
    // special wavelengths
    constexpr double lambdaUV = 1000e-10;  // 1000 Angstrom
    constexpr double lambdaSF = Constants::lambdaSpinFlip();

    // wavelength range outside of which we consider absorption to be zero (approximately 20.47 - 21.74 cm)
    constexpr Range absorptionRange(lambdaSF * (1. - 0.03), lambdaSF * (1. + 0.03));
}

////////////////////////////////////////////////////////////////////

void SpinFlipHydrogenGasMix::setupSelfBefore()
{
    EmittingGasMix::setupSelfBefore();

    // assemble the specific state variables
    _metallicity = _variables.add(StateVariable::metallicity());
    _temperature = _variables.add(StateVariable::temperature());
    _neutralSurfaceDensity = _variables.add(
        StateVariable::custom("neutral hydrogen mass surface density", "masssurfacedensity").named("Sigma_HI"));
    _atomicFraction = _variables.add(StateVariable::custom("atomic hydrogen fraction", "").named("f_atom"));

    auto config = find<Configuration>();
    if (config->hasSecondaryEmission())
    {
        _indexUV = config->radiationFieldWLG()->bin(lambdaUV);
        if (_indexUV < 0) throw FATALERROR("Radiation field wavelength grid does not include 1000 Angstrom");
    }
}

////////////////////////////////////////////////////////////////////

bool SpinFlipHydrogenGasMix::hasExtraSpecificState() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

MaterialMix::DynamicStateType SpinFlipHydrogenGasMix::hasDynamicMediumState() const
{
    return DynamicStateType::Secondary;
}

////////////////////////////////////////////////////////////////////

bool SpinFlipHydrogenGasMix::hasLineEmission() const
{
    return true;
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> SpinFlipHydrogenGasMix::parameterInfo() const
{
    return {SnapshotParameter::metallicity(), SnapshotParameter::temperature(),
            SnapshotParameter::custom("neutral hydrogen mass surface density", "masssurfacedensity", "Msun/pc2")};
}

////////////////////////////////////////////////////////////////////

const StateVariableList& SpinFlipHydrogenGasMix::specificStateVariableInfo() const
{
    return _variables;
}

////////////////////////////////////////////////////////////////////

void SpinFlipHydrogenGasMix::initializeSpecificState(MaterialState* state, const Array& params) const
{
    // leave the properties at zero if the cell does not contain any material for this component
    if (state->numberDensity() > 0.)
    {
        // use the imported values, or for a geometric medium, the default values
        // make sure the temperature is at least the local universe CMB temperature
        bool imported = params.size() > 0;
        _metallicity(state) = imported ? params[0] : defaultMetallicity();
        _temperature(state) = max(Constants::Tcmb(), imported ? params[1] : defaultTemperature());
        _neutralSurfaceDensity(state) = imported ? params[2] : defaultNeutralSurfaceDensity();
        _atomicFraction(state) = 0.;
    }
}

////////////////////////////////////////////////////////////////////

UpdateStatus SpinFlipHydrogenGasMix::updateSpecificState(MaterialState* state, const Array& Jv) const
{
    if (_indexUV < 0) throw FATALERROR("State update should not be called if there is no radiation field");
    UpdateStatus status;

    // if the cell has no neutral hydrogen, then leave the atomic fraction at zero
    if (state->numberDensity() > 0.)
    {
        // get the neutral hydrogen surface mass density
        double SigmaHIpH2 = _neutralSurfaceDensity(state);

        // get the dust-to-gas ratio (metallicity scaled to solar value)
        double DMW = _metallicity(state) / 0.0127;

        // get the radiation field
        // scaled to the reference Milky Way radiation field at 1000 Angstrom
        // converted from 1e6 photons/cm2/s/sr/eV to internal units W/m2/m/sr
        constexpr double h = Constants::h();
        constexpr double c = Constants::c();
        constexpr double Qel = Constants::Qelectron();
        constexpr double JMW = 1e6 * 1e4 * (h * c * h * c / (lambdaUV * lambdaUV * lambdaUV)) / Qel;
        double UMW = Jv[_indexUV] / JMW;

        // get the length scale
        double Lcell = cbrt(state->volume());

        // perform the partitioning scheme
        constexpr double pc100 = 100. * Constants::pc();
        double S = Lcell / pc100;
        double S5 = pow(S, 5);
        double Dstar = 0.17 * (2. + S5) / (1. + S5);
        double g = sqrt(DMW * DMW + Dstar * Dstar);
        constexpr double front = 5e7 * Constants::Msun() / (1e6 * Constants::pc() * Constants::pc());
        double root = sqrt(0.001 + 0.1 * UMW);
        double Sigmac = front * root / (g * (1. + 1.69 * root));
        double alpha = 0.5 + 1. / (1. + sqrt(UMW * DMW * DMW / 600.));
        double Rmol = pow(SigmaHIpH2 / Sigmac, alpha);
        double fmol = Rmol / (Rmol + 1.);

        // set the atomic fraction
        _atomicFraction(state) = max(0., 1. - fmol);
        status.updateConverged();
    }
    return status;
}

////////////////////////////////////////////////////////////////////

double SpinFlipHydrogenGasMix::mass() const
{
    return Constants::Mproton();
}

////////////////////////////////////////////////////////////////////

namespace
{
    // returns the absorption cross section per neutral hydrogen atom for the given wavelength and gas temperature
    double crossSection(double lambda, double T)
    {
        constexpr double front = 3. * M_SQRT2 * M_2_SQRTPI / M_PI / 128. * Constants::EinsteinASpinFlip()
                                 * Constants::h() * Constants::c() * lambdaSF * lambdaSF / Constants::k();
        double Tspin = 6000. * (1 - exp(-0.0002 * T));
        double sigma = sqrt(Constants::k() / Constants::Mproton() * T);
        double u = Constants::c() * (lambda - lambdaSF) / lambda;
        double x = u / sigma;
        double expon = exp(-0.5 * x * x);
        return front / Tspin / sigma * expon;
    }
}

////////////////////////////////////////////////////////////////////

double SpinFlipHydrogenGasMix::sectionAbs(double lambda) const
{
    return absorptionRange.contains(lambda) ? crossSection(lambda, defaultTemperature()) : 0.;
}

////////////////////////////////////////////////////////////////////

double SpinFlipHydrogenGasMix::sectionSca(double /*lambda*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double SpinFlipHydrogenGasMix::sectionExt(double lambda) const
{
    return sectionAbs(lambda);
}

////////////////////////////////////////////////////////////////////

double SpinFlipHydrogenGasMix::opacityAbs(double lambda, const MaterialState* state, const PhotonPacket* /*pp*/) const
{
    if (absorptionRange.contains(lambda))
    {
        double number = state->numberDensity() * _atomicFraction(state);
        if (number > 0.) return crossSection(lambda, _temperature(state)) * number;
    }
    return 0.;
}

////////////////////////////////////////////////////////////////////

double SpinFlipHydrogenGasMix::opacitySca(double /*lambda*/, const MaterialState* /*state*/,
                                          const PhotonPacket* /*pp*/) const
{
    return 0.;
}

////////////////////////////////////////////////////////////////////

double SpinFlipHydrogenGasMix::opacityExt(double lambda, const MaterialState* state, const PhotonPacket* pp) const
{
    return opacityAbs(lambda, state, pp);
}

////////////////////////////////////////////////////////////////////

Array SpinFlipHydrogenGasMix::lineEmissionCenters() const
{
    Array centers(1);
    centers[0] = lambdaSF;
    return centers;
}

////////////////////////////////////////////////////////////////////

Array SpinFlipHydrogenGasMix::lineEmissionMasses() const
{
    Array masses(1);
    masses[0] = Constants::Mproton();
    return masses;
}

////////////////////////////////////////////////////////////////////

Array SpinFlipHydrogenGasMix::lineEmissionSpectrum(const MaterialState* state, const Array& /*Jv*/) const
{
    // calculate the 21 cm luminosity
    constexpr double front = 0.75 * Constants::EinsteinASpinFlip() * Constants::h() * Constants::c() / lambdaSF;
    double L = front * _atomicFraction(state) * state->numberDensity() * state->volume();

    // encapsulate the result in an array
    Array luminosities(1);
    luminosities[0] = L;
    return luminosities;
}

////////////////////////////////////////////////////////////////////

double SpinFlipHydrogenGasMix::indicativeTemperature(const MaterialState* state, const Array& /*Jv*/) const
{
    return _temperature(state);
}

////////////////////////////////////////////////////////////////////
