#include "MFP_hydro_refine.H"
#include "MFP_hydro.H"
#include <AMReX_AmrLevel.H>

ClassFactory<Refinement>& GetHydroRefinementFactory() {
    static ClassFactory<Refinement> F;
    return F;
}

std::string HydroGradientRefinement::tag = "hydro_gradient";
bool HydroGradientRefinement::registered = GetHydroRefinementFactory().Register(HydroGradientRefinement::tag, HydroRefinementBuilder<HydroGradientRefinement>);


HydroGradientRefinement::HydroGradientRefinement(){}
HydroGradientRefinement::HydroGradientRefinement(const int global_idx, const sol::table &def)
{

    idx = global_idx;

    max_level = def.get_or("max_level", -1);
    min_value = def.get_or("min_value", -1);

    // primitive variables
    for (int i = 0; i<HydroState::prim_names.size(); ++i) {
        std::string comp = HydroState::prim_names[i];

        if (def[comp].valid()) {
            prim.push_back(std::make_pair(i,def[comp]));
        }
    }

    // conserved variables
    for (int i = 0; i<HydroState::cons_names.size(); ++i) {
        std::string comp = HydroState::cons_names[i];

        if (def[comp].valid()) {
            cons.push_back(std::make_pair(i,def[comp]));
        }
    }
}

void HydroGradientRefinement::get_tags(MFP* mfp, TagBoxArray& tags) const
{

    if ((mfp->get_level() < max_level) || (max_level < 0)) {

        HydroState& istate = HydroState::get_state_global(idx);

        // grab the conservative state
        const int num_grow = 1;
        MultiFab U(mfp->boxArray(), mfp->DistributionMap(), +HydroDef::ConsIdx::NUM, num_grow, MFInfo(),mfp->Factory());

#ifdef AMREX_USE_EB
        EB2::IndexSpace::push(const_cast<EB2::IndexSpace*>(istate.eb2_index));
#endif

        const Real time = mfp->get_cum_time();

        mfp->FillPatch(*mfp, U, num_grow, time, istate.data_idx, 0, +HydroDef::ConsIdx::NUM);

#ifdef AMREX_USE_EB
        auto const& flags = istate.eb_data.flags;
        auto const& volfrac = istate.eb_data.volfrac;
#endif

        FArrayBox Q;
        const Real* dx = mfp->Geom().CellSize();
        const Real* prob_lo = mfp->Geom().ProbLo();

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(U); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();

#ifdef AMREX_USE_EB
            const FArrayBox& vfrac = volfrac[mfi];
            const EBCellFlagFab& flag = flags[mfi];
            FabType flag_type = flag.getType(bx);
#else
            FabType flag_type = FabType::regular;
#endif

            if (flag_type != FabType::covered) {

                // now go through the list of things to check for high gradients
                for (const auto& info : cons) {
                    const int& icomp = info.first; // index
                    const Real& refine_grad = info.second; // threshold value

                    tag_refinement(bx, U[mfi],
               #ifdef AMREX_USE_EB
                                   flags[mfi],
               #endif
                                   tags[mfi],
                                   icomp,
                                   refine_grad);

                }


                // only calculate the primitives if we need them
                if (!prim.empty()) {

                    Box pbx = grow(bx, 1);

                    Q.resize(pbx, +HydroDef::PrimIdx::NUM);

                    // get primitives
                    istate.calc_primitives(pbx,
                                           U[mfi],
                                           Q,
                                           dx,
                                           time,
                                           prob_lo
                       #ifdef AMREX_USE_EB
                                           ,vfrac
                       #endif
                                           );


                    // now go through the list of things to check for high gradients
                    for (const auto& info : cons) {
                        const int& icomp = info.first; // index
                        const Real& refine_grad = info.second; // threshold value

                        tag_refinement(bx, Q,
               #ifdef AMREX_USE_EB
                                       flags[mfi],
               #endif
                                       tags[mfi],
                                       icomp,
                                       refine_grad);

                    }
                }
            }
        }
    }
}

Real refine_criteria(Array<Real,3> S, Real low_val=1e-12)
{
    Real a = S[2] - 2*S[1] + S[0];
    Real b = std::abs(S[2] - S[1]) + std::abs(S[1] - S[0]) + 0.01*(S[2] + 2*S[1] + S[0]);
    if (std::abs(b) <= low_val) {
        return 0.0;
    } else {
        return std::abs(a/b);
    }
}

void HydroGradientRefinement::tag_refinement(const Box& box,
                                             const FArrayBox& src,
                                             #ifdef AMREX_USE_EB
                                             const EBCellFlagFab& flags,
                                             #endif
                                             TagBox& tags,
                                             const int n,
                                             const Real threshold) const
{
    BL_PROFILE("HydroGradientRefinement::tag_refinement");
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array4<Real const> const& src4 = src.array();
    Array4<char> const& tag = tags.array();

    Real val;
    Array<Real,3> S;
    Array<int,3> index;

#ifdef AMREX_USE_EB
    Array4<const EBCellFlag> const& flag4 = flags.array();
#endif

    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

#ifdef AMREX_USE_EB
                if (flag4(i,j,k).isCovered()) continue;
#endif
                if (tag(i,j,k) == TagBox::SET) continue;

                val = 0.0;
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    // fill our stencil
                    index.fill(0);
                    index[d] = 1;
                    S[0] = src4(i-index[0], j-index[1], k-index[2], n);
                    S[1] = src4(i, j, k, n);
                    S[2] = src4(i+index[0], j+index[1], k+index[2], n);
                    // calculate the refinement criteria and compare it to previous value
                    val = std::max(val, refine_criteria(S, min_value));
                }

                // check against threshold and mark for refinement if necessary
                if (val >= threshold) {
                    tag(i,j,k) = TagBox::SET;
                }
            }
        }
    }
}
