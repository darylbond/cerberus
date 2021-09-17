#include "MFP_eulerian.H"
#include "MFP.H"

EulerianState::EulerianState()
{
    num_grow = 0;
}

EulerianState::~EulerianState(){}

void EulerianState::set_reconstruction()
{

    ClassFactory<Reconstruction> rfact = GetReconstructionFactory();

    sol::table state_def = MFP::lua["states"][name];
    state_def["global_idx"] = global_idx;

    std::string rec = state_def["reconstruction"].get_or<std::string>("null");

    state_def["reconstruction"] = rec; // consistency when using default option "null"

    if (is_transported() && rec == "null")
        Abort("Reconstruction option required for state '"+name+"'. Options are "+vec2str(rfact.getKeys()));

    reconstructor = rfact.Build(rec, state_def);

    if (!reconstructor)
        Abort("Invalid reconstruction option '"+rec+"'. Options are "+vec2str(rfact.getKeys()));

    int ng = reconstructor->get_num_grow();
    set_num_grow(ng);

    return;

}

void EulerianState::set_reflux()
{
    if (AMREX_SPACEDIM < 2) {
        reflux = 0;
    } else {
        reflux = MFP::lua["states"][name]["reflux"].get_or(1);
    }
}

void EulerianState::init_from_lua()
{
    BL_PROFILE("EulerianState::init_from_lua");

    //
    // reflux
    //
    set_reflux();

    //
    // get reconstruction
    //

    set_reconstruction();

}

#ifdef AMREX_USE_EB
Real EulerianState::interp2d(Real cym,
                             Real cy0,
                             Real cyp,
                             Real czm,
                             Real cz0,
                             Real czp,
                             Eigen::Matrix3d v)
{
    return czm*(cym*v(0,0) + cy0*v(1,0) + cyp*v(2,0))
            +     cz0*(cym*v(0,1) + cy0*v(1,1) + cyp*v(2,1))
            +     czp*(cym*v(0,2) + cy0*v(1,2) + cyp*v(2,2));
}


void EulerianState::calc_eb_divergence(const Box& box,
                                       const FArrayBox &cons,
                                       Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                       Array<FArrayBox, AMREX_SPACEDIM> &wall_fluxes,
                                       FArrayBox& du,
                                       const EBCellFlagFab& flag,
                                       const FArrayBox& vfrac,
                                       const Array<const FArrayBox *, AMREX_SPACEDIM> &afrac,
                                       const Array<const FArrayBox *, AMREX_SPACEDIM> &fcent,
                                       int as_crse,
                                       int as_fine,
                                       const IArrayBox *rrflag_as_crse,
                                       const IArrayBox &levmsk, FArrayBox *rr_drho_crse,
                                       FArrayBox &dm_as_fine,
                                       const Real *dx,
                                       const Real dt) const
{
    BL_PROFILE("EulerianState::calc_eb_divergence");
    // make sure du is empty
    du.setVal(0.0);

    int N = du.nComp();

    const Dim3 lo = amrex::lbound(du.box());
    const Dim3 hi = amrex::ubound(du.box());

    Array<int,3> index = {0,0,0};
    Array4<Real> const& du4 = du.array();
    Real dxinv;

    Array4<const EBCellFlag> const& flag4 = flag.array();
    Array4<const Real> const& vfrac4 = vfrac.array();

    Array<Array4<const Real>,AMREX_SPACEDIM> wflux4;
    Array<Array4<const Real>,AMREX_SPACEDIM> flux4;
    Array<Array4<const Real>,AMREX_SPACEDIM> afrac4;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        flux4[d] = fluxes[d].array();
        wflux4[d] = wall_fluxes[d].array();
        afrac4[d] = afrac[d]->array();
    }

    Real wflux, vf;

    for (int n=0; n<N; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {

                    const EBCellFlag &cflag = flag4(i,j,k);
                    vf = vfrac4(i,j,k);

                    if (vf <= 0.0) continue;

                    wflux = 0.0;

                    // cycle over dimensions
                    for (int d=0; d<AMREX_SPACEDIM; ++d) {
                        index[d] = 1;
                        dxinv = dt/dx[d];

                        if (cflag.isSingleValued()) {
                            wflux = wflux4[d](i,j,k,n);
                        }

                        Real flux_lo  = flux4[d](i,j,k,n);
                        const Real alpha_lo = afrac4[d](i,j,k);

                        Real flux_hi  = flux4[d](i+index[0], j+index[1], k+index[2], n);
                        const Real alpha_hi = afrac4[d](i+index[0], j+index[1], k+index[2]);

                        flux_lo *= alpha_lo;
                        flux_hi *= alpha_hi;
                        wflux *= (alpha_hi - alpha_lo);

                        du4(i,j,k,n) += dxinv*(flux_lo - flux_hi + wflux)/vf;
                        index[d] = 0;

                    }
                }
            }
        }
    }

    return;
}

bool is_inside(const int i,const int j, const int k, const Dim3 &lo, const Dim3 &hi)
{
    return  i >= lo.x && i <= hi.x &&
            j >= lo.y && j <= hi.y &&
            k >= lo.z && k <= hi.z;
}

void EulerianState::merge_cells(const Box& box,
                                FArrayBox &cons,
                                FArrayBox& du,
                                const EBCellFlagFab& flag,
                                const FArrayBox& vfrac,
                                const Array<const FArrayBox *, AMREX_SPACEDIM> &afrac,
                                int as_fine,
                                FArrayBox &dm_as_fine,
                                const IArrayBox& levmsk) const
{
    BL_PROFILE("EulerianState::merge_cells");
    int nc = du.nComp();

    Array<int,3> index = {0,0,0};

    Array4<Real> const& cons4 = cons.array();
    Array4<Real> const& du4 = du.array();

    Array4<const EBCellFlag> const& flag4 = flag.array();
    Array4<const Real> const& vfrac4 = vfrac.array();
    Array4<Real> const& dm_as_fine4 = dm_as_fine.array();
    Array4<const int> const& levmsk4 = levmsk.array();

    Array<Array4<const Real>,AMREX_SPACEDIM> afrac4;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        afrac4[d] = afrac[d]->array();
    }

    Real vf;

    // the sets of cells that are to be merged
    Vector<Vector<Array<int,3>>> merge;

    // a container for all of the unique cells that are part of the problem and
    // the super-cell that they are a part of
    std::map<Array<int,3>,int> cells;


    int ii, jj, kk;
    int mi, mj, mk;

    Real side_alpha, side_alpha_max;

    Box halo = grow(box,2);
    const Dim3 halo_lo = amrex::lbound(halo);
    const Dim3 halo_hi = amrex::ubound(halo);

    // expand zone over which we work so that each block sees the same modifications
    // to du without having to do inter-block communication
    Box calc = grow(box,1);
    const Dim3 lo = amrex::lbound(calc);
    const Dim3 hi = amrex::ubound(calc);


    for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {

                const EBCellFlag &cflag = flag4(i,j,k);

                vf = vfrac4(i,j,k);

                if (vf <= 0.0) continue;

                if (vf < merge_threshold) {

                    ii = i;
                    jj = j;
                    kk = k;

                    Array<int,3> cell_idx = {ii, jj, kk};

                    // check that this cell isn't already part of another super-cell
                    if (cells.find(cell_idx) != cells.end())
                        continue;

                    Vector<Array<int,3>> super_cell = {cell_idx};
                    bool new_super_cell = true;

                    while (vf < merge_threshold) {
                        // get the fluid fraction of each of the sides
                        // scaled by if it is a viable candidate
                        side_alpha_max = 0.0;
                        for (int d=0; d<AMREX_SPACEDIM; ++d) {
                            index[d] = 1;

                            // lo side

                            if (is_inside(ii-index[0], jj-index[1], kk-index[2], halo_lo, halo_hi)) {
                                side_alpha = afrac4[d](ii,jj,kk);
                                if (side_alpha > side_alpha_max) {
                                    side_alpha_max = side_alpha;
                                    mi = ii-index[0];
                                    mj = jj-index[1];
                                    mk = kk-index[2];

                                }
                            }

                            // hi side

                            if (is_inside(ii+index[0], jj+index[1], kk+index[2], halo_lo, halo_hi)) {
                                side_alpha = afrac4[d](ii+index[0], jj+index[1], kk+index[2]);
                                if (side_alpha > side_alpha_max) {
                                    side_alpha_max = side_alpha;
                                    mi = ii+index[0];
                                    mj = jj+index[1];
                                    mk = kk+index[2];

                                }
                            }

                            index[d] = 0;
                        }

                        ii = mi;
                        jj = mj;
                        kk = mk;

                        Array<int,3> cell_idx = {ii, jj, kk};

                        // check if this cell is part of another super cell
                        if (cells.find(cell_idx) == cells.end()) {
                            vf += vfrac4(ii,jj,kk);
                            super_cell.push_back(cell_idx);
                        } else {
                            int si = cells[cell_idx];

                            // register the cells with a super cell
                            for (const auto &cell : super_cell) {
                                merge[si].push_back(cell);
                                cells[cell] = si;
                            }
                            new_super_cell = false;
                            break;
                        }
                    }

                    // add the merged cell to the list
                    if (new_super_cell) {

                        for (const auto &cell : super_cell) {
                            cells[cell] = merge.size();
                        }

                        merge.push_back(super_cell);
                    }



                } else if (cflag.isSingleValued() && as_fine) {
                    for (int n=0; n<nc; ++n) {
                        dm_as_fine4(i,j,k,n) = 0.0;
                    }
                }
            }
        }
    }

    int nsuper = merge.size();

    // do update

    for (int n=0; n<nc; ++n) {

        for (int mi=0; mi<nsuper; ++mi) {
            const auto& mc = merge[mi];


            // compute the sums
            Real vf_sum   = 0.0;
            Real sum_cons = 0.0;
            Real sum_du   = 0.0;

            for (const auto& cell_idx : mc) {
                vf = vfrac4(cell_idx[0],cell_idx[1], cell_idx[2]);
                vf_sum += vf;

                sum_cons += cons4(cell_idx[0],cell_idx[1], cell_idx[2],n)*vf;
                sum_du += du4(cell_idx[0],cell_idx[1], cell_idx[2],n)*vf;
            }

            sum_cons /= vf_sum;
            sum_du /= vf_sum;

            for (const auto& cell_idx : mc) {

                const int i = cell_idx[0];
                const int j = cell_idx[1];
                const int k = cell_idx[2];

                const Real u_orig = cons4(i,j,k,n);
                const Real u_final = sum_cons + sum_du;
                const Real du_merge = u_final - u_orig;

                du4(i,j,k,n) = du_merge;

                // tell other grids how du has changed
                if (as_fine && (levmsk4(i,j,k) == LevelMask::NotCovered)) {
                    dm_as_fine4(i,j,k,n) = du_merge - du4(cell_idx[0],cell_idx[1], cell_idx[2],n);
                }
            }
        }

    }
}

#endif

void EulerianState::calc_divergence(const Box& box,
                                    Array<FArrayBox, AMREX_SPACEDIM> &fluxes,
                                    FArrayBox& du,
                                    const Real *dx,
                                    const Real dt) const
{
    BL_PROFILE("State::calc_divergence");
    // make sure du is empty
    du.setVal(0.0);

    int N = du.nComp();

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array<int,3> index = {0,0,0};
    Array4<Real> const& du4 = du.array();
    Real dxinv;

    Array<Array4<const Real>,AMREX_SPACEDIM> flux4;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        flux4[d] = fluxes[d].array();
    }

    for (int n=0; n<N; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {

                    // cycle over dimensions
                    for (int d=0; d<AMREX_SPACEDIM; ++d) {
                        index[d] = 1;
                        dxinv = dt/dx[d];

                        const Real& flux_lo  = flux4[d](i,j,k,n);

                        const Real& flux_hi  = flux4[d](i+index[0], j+index[1], k+index[2], n);

                        du4(i,j,k,n) += dxinv*(flux_lo - flux_hi);


                        index[d] = 0;
                    }
                }
            }
        }
    }

    return;
}

void EulerianState::write_info(nlohmann::json &js) const
{

    State::write_info(js);

    // write out stuff that is common to all states

    js["num_grow"] = num_grow;

}
