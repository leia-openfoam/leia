#include "leiaVersionRegistry.H"
#include "Time.H"
#include "OFstream.H"
#include "Pstream.H"

namespace Foam
{
namespace leia
{

// Function-local static: constructed on first use, so registrars in other
// translation units are safe whatever the static initialisation order.
static DynamicList<Tuple2<word, string>>& registry()
{
    static DynamicList<Tuple2<word, string>> r;
    return r;
}

versionRegistrar::versionRegistrar(const char* library, const char* stamp)
{
    registry().push_back(Tuple2<word, string>(word(library), string(stamp)));
}

const DynamicList<Tuple2<word, string>>& versions()
{
    return registry();
}

void reportVersions(Ostream& os)
{
    for (const auto& v : registry())
    {
        os << "leia library " << v.first() << " : " << v.second().c_str() << nl;
    }
    os << endl;
}

void writeVersions(const Time& runTime)
{
    if (!Pstream::master()) return;
    OFstream os(runTime.globalPath()/"leia.version");
    for (const auto& v : registry())
    {
        os << v.first() << ' ' << v.second().c_str() << nl;
    }
}

} // End namespace leia
} // End namespace Foam
