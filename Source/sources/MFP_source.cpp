#include "MFP_source.H"

Source::Source(){}
Source::~Source(){}

ClassFactory<Source>& GetSourceFactory() {
    static ClassFactory<Source> F;
    return F;
}
