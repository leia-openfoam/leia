/*---------------------------------------------------------------------------*\
    leia :: davof :: davofCurvature

    Curvature of the DAVOF state, run-time selectable by the DAVOF applications
    and living in libleiaDavof: nothing here depends on the level-set
    libraries.

        fvSolution davof
        {
            curvature
            {
                models (normalsOnly hermite);          // evaluated in this order
                normalsOnly { type quadricFit; positionWeight 0; }
                hermite     { type quadricFit; positionWeight 1; }
            }
        }

    Data per interface cell (davofState): the plane-polygon centroid xPlane_c
    (second order) and the area normal n_c = m_c/|m_c| (the order q of the face
    fractions). A model fits a local surface to the ring of these Hermite data
    and reports, AT THE CENTROID, the total curvature

        kappa = kappa_1 + kappa_2 = div(n)     (2/R on a sphere; n out of the
                                                liquid, positive for a droplet),

    the Gaussian curvature K = kappa_1 kappa_2 and the principal curvatures.
    It also answers the OFFSET query for a point x (a cell or a face centre):
    the signed distance d of x from the cell's local surface model along the
    model's normal, positive on the gas side (leia's psi convention), and the
    model's kappa and K at the ray foot. That is what the delivery to cell
    centres and faces (davofCurvatureDelivery.H) needs for the parallel-surface
    closed form kappa_i(d) = kappa_i/(1 + d kappa_i): the interface centroid is
    neither the cell centre nor a face centre, so the value is converted, as
    the level-set solver converts its contour curvature, not extended.

    A model packs itself into a fixed-size scalar list (davofModelPack) so the
    delivery can evaluate the model of a cell on the other side of a processor
    or cyclic boundary after syncTools::swapBoundaryFaceList: an np-4 run gives
    the serial numbers.

    Models: none (this base: nothing filled) and quadricFit
    (davofQuadricCurvature.H).
\*---------------------------------------------------------------------------*/

#ifndef davofCurvature_H
#define davofCurvature_H

#include "fvMesh.H"
#include "volFields.H"
#include "dictionary.H"
#include "FixedList.H"
#include "runTimeSelectionTables.H"
#include "davofState.H"

namespace Foam
{

//- A model packed for exchange across coupled patches (layout: the model's).
typedef FixedList<scalar, 20> davofModelPack;

class davofCurvature
{
protected:

    const word name_;
    const fvMesh& mesh_;
    const davofState& state_;
    const dictionary dict_;

    volScalarField kappaS_;     // kappa.davof.<name>, total curvature at xPlane [1/m]
    volScalarField KS_;         // K.davof.<name>, Gaussian curvature [1/m^2]
    volScalarField kappa1S_;    // kappa1.davof.<name>
    volScalarField kappa2S_;    // kappa2.davof.<name>
    boolList filled_;           // a model exists in this cell
    label nFallback_;           // cells whose model fell back (kappa = K = 0)
    scalar cpuSeconds_;

public:

    TypeName("none");

    declareRunTimeSelectionTable
    (
        autoPtr,
        davofCurvature,
        Dictionary,
        (
            const word& name,
            const fvMesh& mesh,
            const davofState& state,
            const dictionary& dict
        ),
        (name, mesh, state, dict)
    );

    davofCurvature
    (
        const word& name,
        const fvMesh& mesh,
        const davofState& state,
        const dictionary& dict
    );

    //- Select by dict.type (default none)
    static autoPtr<davofCurvature> New
    (
        const word& name,
        const fvMesh& mesh,
        const davofState& state,
        const dictionary& dict
    );

    virtual ~davofCurvature() = default;

    //- Fit the models and fill kappa, K, kappa1, kappa2 on the interface cells.
    virtual void compute() {}

    //- Pack the model of cell c; false when the cell has none.
    virtual bool pack(const label c, davofModelPack& p) const
    {
        return false;
    }

    //- The offset query on a packed model: signed distance d of x from the
    //  model surface along the model normal (positive on the gas side) and the
    //  model's total and Gaussian curvature at the ray foot. False when the
    //  pack holds no model.
    virtual bool offsetPacked
    (
        const davofModelPack& p,
        const point& x,
        scalar& d,
        scalar& kappa,
        scalar& K
    ) const
    {
        return false;
    }

    //- The offset query on the model of cell c.
    bool offset
    (
        const label c,
        const point& x,
        scalar& d,
        scalar& kappa,
        scalar& K
    ) const
    {
        davofModelPack p;
        if (!pack(c, p)) return false;
        return offsetPacked(p, x, d, kappa, K);
    }

    const word& name() const { return name_; }
    const davofState& state() const { return state_; }
    const volScalarField& kappa() const { return kappaS_; }
    const volScalarField& K() const { return KS_; }
    const volScalarField& kappa1() const { return kappa1S_; }
    const volScalarField& kappa2() const { return kappa2S_; }
    const boolList& filled() const { return filled_; }
    label nFallback() const { return nFallback_; }
    scalar cpuSeconds() const { return cpuSeconds_; }

    //- Write the four fields.
    void write() const;
};

} // namespace Foam

#endif

// ************************************************************************* //
