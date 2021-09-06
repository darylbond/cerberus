#ifndef MFP_HYDRO_VISCOUS_H
#define MFP_HYDRO_VISCOUS_H

#include "MFP.H"
#include "MFP_hydro_defs.H"

class HydroViscous;

class HydroState;
#include "MFP_hydro.H"

//-------------------------------------------------------------------------


class HydroViscous
{
public:
    HydroViscous();
    virtual ~HydroViscous() = default;

    enum DiffusionType {
      Neutral = 0,
      Ion,
      Electron
    };

    virtual int get_type();
    virtual int get_num();
    virtual void update_linked_states();
    virtual Vector<int> get_linked_states() {return {idx};}

    virtual const std::string& get_tag() const = 0;
    virtual const std::map<std::string, Real> get_refs() const {return {};}
    virtual std::string str() const;

    virtual Real get_min_dt(MFP *mfp) const = 0;

    int idx;
    Real cfl;
};

template <typename Derived>
class AbstractHydroViscous : public HydroViscous {
public:

    virtual Real get_min_dt(MFP *mfp) const override
    {
        BL_PROFILE("get_max_speed");
        const HydroState& istate = HydroState::get_state(idx);

        MultiFab& data = mfp->get_new_data(istate.data_idx);

        Real min_dt = std::numeric_limits<Real>::max();

        const Real* dx = mfp->Geom().InvCellSize();
        const Real dx2 = AMREX_D_TERM(dx[0]*dx[0],+dx[1]*dx[1],+dx[2]*dx[2]);

        Array<Real,+HydroDef::ConsIdx::NUM> U;

        for (MFIter mfi(data); mfi.isValid(); ++mfi) {
            const Box& box = mfi.tilebox();
            const Dim3 lo = amrex::lbound(box);
            const Dim3 hi = amrex::ubound(box);

    #ifdef AMREX_USE_EB
            // get the EB data required for later calls
            const FArrayBox& vfrac = istate.eb_data.volfrac[mfi];

            if (vfrac.getType() == FabType::covered) continue;

            Array4<const Real> const& vf4 = vfrac.array();
    #endif

            Array4<Real const> const dat = data.const_array(mfi);

            for     (int k = lo.z; k <= hi.z; ++k) {
                for   (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

    #ifdef AMREX_USE_EB
                        if (vf4(i,j,k) == 0.0) {
                            continue;
                        }
    #endif

                        for (int n=0; n<+HydroDef::ConsIdx::NUM; ++n) {
                            U[n] = dat(i,j,k,n);
                        }

                        const Real dt_ = static_cast<Derived const&>(*this).calc_dt(U, dx2);

                        min_dt = amrex::min(min_dt, dt_);
                    }
                }
            }
        }

        return min_dt;
    }

protected:
   // We make clear AbstractHydroViscous class needs to be inherited
   AbstractHydroViscous() = default;
   AbstractHydroViscous(const AbstractHydroViscous&) = default;
   AbstractHydroViscous(AbstractHydroViscous&&) = default;

};

template <typename D>
std::unique_ptr<HydroViscous> ViscousBuilder(const sol::table& def)
{
    std::unique_ptr<HydroViscous> viscous;

    const int global_idx = def["global_idx"];
    const sol::table& visc = def["viscosity"].get_or(sol::table());

    if (visc.valid()) {

        // load the diffusion coefficients

        if (visc["type"] == D::tag) {
            return std::unique_ptr<D>(new D(global_idx, visc));
        } else {
            return nullptr;
        }
    } else {
        return nullptr;
    }
}

ClassFactory<HydroViscous>& GetViscousFactory();

// ====================================================================================

class Sutherland : public AbstractHydroViscous<Sutherland>
{
public:
    Sutherland();
    Sutherland(const int global_idx, const sol::table &def);
    ~Sutherland();

    Real mu_0, T0, S, Prandtl;

    virtual int get_type() override;
    virtual int get_num() override;
    virtual const std::string& get_tag() const override {return tag;}
    virtual const std::map<std::string, Real> get_refs() const override {
        return {{"mu_0", mu_0},{"T0", T0}, {"S",S}, {"Pr",Prandtl}};
    }

    Real calc_dt(Array<Real,+HydroDef::ConsIdx::NUM>& U, const Real dx2) const;
    void get_coeffs(const Array<Real, +HydroDef::PrimIdx::NUM> &Q, Real &T, Real &mu, Real &kappa);
    static bool valid_state(const int idx);

    static bool registered;
    static std::string tag;
};

// ====================================================================================

class PowerLaw : public AbstractHydroViscous<Sutherland>
{
public:
    PowerLaw();
    PowerLaw(const int global_idx, const sol::table& def);
    ~PowerLaw();

    enum NeutralDiffusiveVariable {
        Temp = 0,
        Mu,
        Kappa,
        NUM_NEUTRAL_DIFF
    };

    Real mu_0, T0, n, Prandtl;

    virtual int get_type() override;
    virtual int get_num() override;
    virtual const std::string& get_tag() const override {return tag;}
    virtual const std::map<std::string, Real> get_refs() const override {
        return {{"mu_0", mu_0},{"T0", T0}, {"n",n}, {"Pr",Prandtl}};
    }

    Real calc_dt(Array<Real,+HydroDef::ConsIdx::NUM>& U, const Real dx2) const;
    void get_coeffs(const Array<Real,+HydroDef::PrimIdx::NUM>& Q, Real &T, Real &mu, Real &kappa);
    static bool valid_state(const int idx);

    static bool registered;
    static std::string tag;
};

// ====================================================================================

class UserDefinedViscosity : public AbstractHydroViscous<Sutherland>
{
public:
    UserDefinedViscosity();
    UserDefinedViscosity(const int global_idx, const sol::table& def);
    ~UserDefinedViscosity();

    enum NeutralDiffusiveVariable {
        Temp = 0,
        Mu,
        Kappa,
        NUM_NEUTRAL_DIFF
    };

    Real mu_0, T0, n, Prandtl;

    virtual int get_type() override;
    virtual int get_num() override;
    virtual const std::string& get_tag() const override {return tag;}
    virtual const std::map<std::string, Real> get_refs() const override {
        return {{"mu_0", mu_0}, {"Pr",Prandtl}};
    }

    Real calc_dt(Array<Real,+HydroDef::ConsIdx::NUM>& U, const Real dx2) const;
    void get_coeffs(const Array<Real,+HydroDef::PrimIdx::NUM>& Q, Real &T, Real &mu, Real &kappa);
    static bool valid_state(const int idx);

    static bool registered;
    static std::string tag;
};

// ====================================================================================

#endif // MFP_HYDRO_VISCOUS_H
