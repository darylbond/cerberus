#include "MFP_argon.H"
#include "MFP.H"
#include "MFP_state.H"
#include "sol.hpp"

std::string CollisionsArgon::tag = "argon_ionisation";
bool CollisionsArgon::registered = GetActionFactory().Register(CollisionsArgon::tag, ActionBuilder<CollisionsArgon>);

CollisionsArgon::CollisionsArgon(){}
CollisionsArgon::~CollisionsArgon(){}

CollisionsArgon::CollisionsArgon(const int idx, const sol::table &def) : Collisions(idx, def)
{
    // get some shortcuts for the species
    for (const SpeciesInfo& info : species_info) {
        if (info.q > 0) {
            index[2] = info.idx;
        } else if (info.q < 0) {
            index[1] = info.idx;
        } else {
            index[0] = info.idx;
        }
    }

    constexpr Real Av = 6.0221409e23;

    // {Ar, e-, Ar+}
    C = {10.12, 22.59e4, 10.12};
    eta = {1.5, 1.5, 1.5};
    theta = {135300, 135300, 135300};
    CK = {2.9e22, 2.9e22, 2.9e22};
    zeta = {1.5, 1.5, 1.5};
    phi = {183100, 183100, 183100};
    m = {6.628e-26, 9.108e-31, 0.0}; m[2]=m[0]-m[1];
    M = {Av*m[0], Av*m[1], Av*m[2]};
    d = {4e-10, 0, 4e-10};
    Z = {0, -1, 1};

    // do some checks on mass, charge, gamma here!

    return;
}

void CollisionsArgon::calc_update()
{

    // reaction rates
    Array<Real,+SpeciesIndex::NUM> kf, kb, Ke, X, rho;



    for (int s=0; s<+SpeciesIndex::NUM; ++s) {
        SpeciesInfo& info = species_info[index[s]];
        Ke[s] = CK[s]*pow(info.T, zeta[s]);
        kf[s] = C[s]*pow(info.T, eta[s])(theta[s]/info.T+2)*exp(-theta[s]/info.T);
        kb[s] = kf[s]/Ke[s];

        rho[s] = MFP::rho_ref*info.rho; // << note change of units

        X[s] = rho[s]/M[s];

    }

    // calculate the rate of reaction
    Real R = 0.0;
    for (int s=0; s<+SpeciesIndex::NUM; ++s) {
        R += -kf[s]*X[+SpeciesIndex::Ar]*X[s] + kb[s]*X[+SpeciesIndex::Ar_plus]*X[+SpeciesIndex::e_minus]*X[s];
    }

    // mass production terms
    Array<Real,+SpeciesIndex::NUM> w;
    w[+SpeciesIndex::Ar]      =  M[+SpeciesIndex::Ar]*R;
    w[+SpeciesIndex::Ar_plus] = -M[+SpeciesIndex::Ar_plus]*R;
    w[+SpeciesIndex::e_minus] = -M[+SpeciesIndex::e_minus]*R;


    // momentum transfer
    Array<Real,+SpeciesIndex::NUM> Qu;
    Real nu, sigma, r, cr;

    constexpr Real kB = 1.38064852e-23;

    // Argon & electrons
    const Real r_ne = 0.5*(d[+SpeciesIndex::Ar] + d[+SpeciesIndex::e_minus]);
    const Real sigma_ne = PI*r*r;
    const Real T_n = species_info[index[+SpeciesIndex::Ar]].T;
    const Real c2_n = 8*kB*T_n/(PI*m[+SpeciesIndex::Ar]);
    const Real T_e = species_info[index[+SpeciesIndex::Ar]].T;
    const Real c2_e = 8*kB*T_e/(PI*m[+SpeciesIndex::e_minus]);
    nu = rho[+SpeciesIndex::Ar]*sigma*cr

    Qu[+SpeciesIndex::Ar] = rho[+SpeciesIndex::Ar]*nu




    size_t n_species = species_info.size();

    for (int a = 0; a < n_species; ++a) {

        SpeciesInfo& sda = species_info[a];

        for (int b = a+1; b < n_species; ++b) {

            SpeciesInfo& sdb = species_info[b];

            const Real m_ab = (sda.m*sdb.m)/(sda.m + sdb.m);

            Real du2 = 0.0;
            Array<Real,3> du;
            for (int i=0; i<3; ++i) {
                du[i] = sdb.vel[i] - sda.vel[i];
                du2 += du[i]*du[i];
            }

            // collision frequency
            Real nu;
            if ((sda.q != 0) && (sdb.q != 0)) {
                const Real coeff_1 = c1*((sda.q2*sdb.q2*lnC)/(m_ab*sda.m));
                nu = sdb.n*coeff_1*pow(c2*du2 + sda.T/sda.m + sdb.T/sdb.m, -1.5);
            } else {
                const Real coeff_2 = (sdb.m*ccs[sda.idx][sdb.idx])/(sda.m + sdb.m);
                nu = sdb.n*coeff_2*c3*sqrt(sda.T/sda.m + sdb.T/sdb.m);
            }


            const Real S = 0.0;

            // energy exchange
            Real Q = m_ab*sda.n*nu*du2 + 3*sda.rho*nu/(sda.m + sdb.m)*(sdb.T - sda.T);

            // momentum exchange
            Array<Real,3> R;
            for (int n=0; n<3; ++n) {
                R[n] = sda.rho*nu*du[n];
                Q += R[n]*sda.vel[n];
            }

            // contribution to change of density
            sda.delta[+HydroDef::ConsIdx::Density] += S;
            sdb.delta[+HydroDef::ConsIdx::Density] -= S;

            // contribution to change of momentum
            for (int n=0; n<3; ++n) {
                sda.delta[+HydroDef::ConsIdx::Xmom+n] += R[n];
                sdb.delta[+HydroDef::ConsIdx::Xmom+n] -= R[n];
            }

            // contribution to change of energy
            sda.delta[+HydroDef::ConsIdx::Eden] += Q;
            sdb.delta[+HydroDef::ConsIdx::Eden] -= Q;
        }
    }
}
