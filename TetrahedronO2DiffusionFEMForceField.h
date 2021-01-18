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
    TetrahedronO2DiffusionFEMForceField(){}
    /// Destructor
    virtual ~TetrahedronO2DiffusionFEMForceField(){};


};
}  // namespace sofa::component::forcefield

#endif //SOFAOXYGENDIFFUSION_TETRAHEDRONO2DIFFUSIONFEMFORCEFIELD_H
