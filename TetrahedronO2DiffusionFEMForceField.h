//
// Created by Somers on 1/18/2021.
//

#ifndef SOFAOXYGENDIFFUSION_TETRAHEDRONO2DIFFUSIONFEMFORCEFIELD_H
#define SOFAOXYGENDIFFUSION_TETRAHEDRONO2DIFFUSIONFEMFORCEFIELD_H

#include <SofaSimpleFem/TetrahedronDiffusionFEMForceField.h>

namespace sofa::component::forcefield {

template<class DataTypes>
class TetrahedronO2DiffusionFEMForceField : public TetrahedronDiffusionFEMForceField<DataTypes>
{
public:
    /// Constructor
    TetrahedronO2DiffusionFEMForceField();
    /// Destructor
    virtual ~TetrahedronO2DiffusionFEMForceField();

    /// Vector of densities for each tetrahedron
    Data<sofa::helper::vector<Real> > d_tetraNormoxicDensity;

    /// Vector of densities for each tetrahedron
    Data<sofa::helper::vector<Real> > d_tetraHypoxicDensity;

    Data<Real> d_alphaN;

    Data<Real> d_alphaH;

    Data<Real> d_O2SaturationConstant;
    Data<sofa::helper::vector<Real> > d_tetraUptakeCoefficient;


protected:
    /// Function computing the edge diffusion coefficient from tetrahedral information
    void computeEdgeDiffusionCoefficient();
    /// Forcefield functions for Matrix system. Adding force to global forcefield vector.
    void addForce (const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& dF, const DataVecCoord& dX, const DataVecDeriv& /*v*/);

    sofa::helper::vector<Real> edgeUptakeCoefficient;

};
}  // namespace sofa::component::forcefield

#include "TetrahedronO2DiffusionFEMForceField.inl"
#endif //SOFAOXYGENDIFFUSION_TETRAHEDRONO2DIFFUSIONFEMFORCEFIELD_H
