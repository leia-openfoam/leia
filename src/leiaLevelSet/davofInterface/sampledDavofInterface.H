/*---------------------------------------------------------------------------*\
    leia :: davofInterface :: sampledDavofInterface

    A sampledSurface (the `surfaces` function object) that draws the DAVOF
    plane of every interface cell from two registered volVectorFields: the
    area normal pointing out of the liquid (m.davof, |m_c| the interface
    area) and a point on the plane (xPlane.davof, the polygon centroid that
    davofState::planePosition() writes). Any pair of fields with that meaning
    can be named, e.g. a geometricVoF scheme's -normal_ and centre_.

        surfaces
        {
            type            surfaces;
            libs            (leiaDavofInterface);
            surfaceFormat   vtk;
            fields          (alpha.davof);
            surfaces
            (
                davof
                {
                    type      davofInterface;
                    normal    m.davof;        // default
                    centre    xPlane.davof;   // default
                    wispTol   1e-3;           // interface cell: |m_c| > wispTol V_c^(2/3), default 0
                    interpolate false;
                }
            );
        }

    Ported from TwoPhaseFlow's src/postProcessing/interface/sampledInterface
    (Henning Scheufler, DLR, GPL-3); registered as `davofInterface` so it can
    coexist with OpenFOAM's own `interface` (libgeometricVoF) and
    TwoPhaseFlow's `reconstructedInterface`.
\*---------------------------------------------------------------------------*/

#ifndef sampledDavofInterface_H
#define sampledDavofInterface_H

#include "sampledSurface.H"
#include "plicInterfaceSurface.H"

namespace Foam
{

class sampledDavofInterface
:
    public sampledSurface
{
    // Private data

        word normalName_;
        word centreName_;
        scalar wispTol_;

        mutable autoPtr<plicInterfaceSurface> surfPtr_;

        //- Time index at the last rebuild
        mutable label prevTimeIndex_;


    // Private Member Functions

        //- Rebuild the surface when the time index changed
        bool updateGeometry() const;

        template<class Type>
        tmp<Field<Type>> sampleOnFaces
        (
            const interpolation<Type>& sampler
        ) const;

        template<class Type>
        tmp<Field<Type>> sampleOnPoints
        (
            const interpolation<Type>& interpolator
        ) const;


public:

    TypeName("sampledDavofInterface");


    // Constructors

        sampledDavofInterface
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    virtual ~sampledDavofInterface() = default;


    // Member Functions

        const plicInterfaceSurface& surface() const
        {
            return surfPtr_();
        }

        virtual bool needsUpdate() const;
        virtual bool expire();
        virtual bool update();

        virtual const pointField& points() const
        {
            return surface().points();
        }

        virtual const faceList& faces() const
        {
            return surface().surfFaces();
        }

        virtual const labelList& zoneIds() const
        {
            return Foam::emptyLabelList;
        }

        virtual const vectorField& Sf() const
        {
            return surface().Sf();
        }

        virtual const scalarField& magSf() const
        {
            return surface().magSf();
        }

        virtual const vectorField& Cf() const
        {
            return surface().Cf();
        }


        // Sample

            virtual tmp<scalarField> sample
            (
                const interpolation<scalar>& sampler
            ) const;

            virtual tmp<vectorField> sample
            (
                const interpolation<vector>& sampler
            ) const;

            virtual tmp<sphericalTensorField> sample
            (
                const interpolation<sphericalTensor>& sampler
            ) const;

            virtual tmp<symmTensorField> sample
            (
                const interpolation<symmTensor>& sampler
            ) const;

            virtual tmp<tensorField> sample
            (
                const interpolation<tensor>& sampler
            ) const;


        // Interpolate

            virtual tmp<scalarField> interpolate
            (
                const interpolation<scalar>& interpolator
            ) const;

            virtual tmp<vectorField> interpolate
            (
                const interpolation<vector>& interpolator
            ) const;

            virtual tmp<sphericalTensorField> interpolate
            (
                const interpolation<sphericalTensor>& interpolator
            ) const;

            virtual tmp<symmTensorField> interpolate
            (
                const interpolation<symmTensor>& interpolator
            ) const;

            virtual tmp<tensorField> interpolate
            (
                const interpolation<tensor>& interpolator
            ) const;


    // Output

        virtual void print(Ostream& os, int level = 0) const;
};

} // End namespace Foam

#ifdef NoRepository
    #include "sampledDavofInterfaceTemplates.C"
#endif

#endif

// ************************************************************************* //
