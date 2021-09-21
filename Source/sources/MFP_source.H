#ifndef SOURCE_H
#define SOURCE_H

#include <AMReX_Vector.H>
#include "sol.hpp"
#include "MFP_factory.H"

class Source;
class MFP;

using namespace amrex;

class Source
{
public:
    Source();
    ~Source();

    int src_idx;
    std::string name;

    virtual void solve(MFP* mfp, const Real dt) const = 0;

    virtual std::string get_tag() const = 0;
};

template <typename D>
std::unique_ptr<Source> SourceBuilder(const sol::table& def)
{
    int idx = def["src_idx"];

    if (def["type"] == D::tag) {
        return std::unique_ptr<D>(new D(idx, def));
    } else {
        return nullptr;
    }
}

ClassFactory<Source> &GetSourceFactory();

#endif // SOURCE_H
