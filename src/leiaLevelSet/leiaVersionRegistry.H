/*---------------------------------------------------------------------------*\
    leia :: version registry

    Every leia library registers, at load time, its name and the
    `git describe --always[-dirty]` of the source tree it was built from
    (the registrar object lives in the generated leiaStamp_<library>.C that
    etc/leia-stamp.sh writes before each wmake). A solver prints one line per
    loaded library in its startup banner and writes the same lines to
    <case>/leia.version, so a binary that does not match its clone shows up in
    every log and every case, next to the gitCommit column of the workflow.
    MEASURED 2026-09-09: a library about 200 commits ahead of its clone ran for
    two weeks without any log showing it (STATUS.md section 9).
\*---------------------------------------------------------------------------*/

#ifndef leiaVersionRegistry_H
#define leiaVersionRegistry_H

#include "word.H"
#include "string.H"
#include "Tuple2.H"
#include "DynamicList.H"
#include "Ostream.H"

namespace Foam
{
class Time;

namespace leia
{

//- One static object per library, constructed from the generated stamp file.
struct versionRegistrar
{
    versionRegistrar(const char* library, const char* stamp);
};

//- (library, stamp) pairs of every loaded leia library, in load order.
const DynamicList<Tuple2<word, string>>& versions();

//- One banner line per library: "leia library <library> : <stamp>".
void reportVersions(Ostream& os);

//- Write the same lines to <case>/leia.version (master rank only).
void writeVersions(const Time& runTime);

} // End namespace leia
} // End namespace Foam

#endif
