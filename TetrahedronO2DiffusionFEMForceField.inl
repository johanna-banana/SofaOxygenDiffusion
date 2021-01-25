//
// Created by Somers on 1/18/2021.
//

#include "TetrahedronO2DiffusionFEMForceField.h"
#pragma once
#include <SofaSimpleFem/TetrahedronDiffusionFEMForceField.h>

#include <sofa/core/visual/VisualParams.h>

#include <SofaBaseTopology/TopologyData.inl>
#include <SofaBaseTopology/TopologySubsetData.inl>

#include <sofa/core/behavior/ForceField.inl>
#include <sofa/core/behavior/MultiMatrixAccessor.h>

#include <sofa/helper/AdvancedTimer.h>

namespace sofa::component::forcefield {

using namespace sofa::defaulttype;
using namespace sofa::component::topology;
using namespace core::topology;
using namespace core::objectmodel;
using core::topology::BaseMeshTopology;
typedef BaseMeshTopology::EdgesInTetrahedron EdgesInTetrahedron;

// --------------------------------------------------------------------------------------
// --- constructor
// --------------------------------------------------------------------------------------
template<class DataTypes>
TetrahedronO2DiffusionFEMForceField<DataTypes>::TetrahedronO2DiffusionFEMForceField()
        : TetrahedronDiffusionFEMForceField<DataTypes>::TetrahedronDiffusionFEMForceField(),  // initialize parent class
          d_tetraNormoxicDensity(initData(&d_tetraNormoxicDensity, "tetraNormoxicDensity","Density of normoxic cells for each tetrahedron, by default equal to 0.")),
          d_tetraHypoxicDensity(initData(&d_tetraHypoxicDensity, "tetraHypoxicDensity","Density of hypoxic cells for each tetrahedron, by default equal to 0.")),
          d_alphaN(initData(&d_alphaN, "alphaN","Oxygen uptake coefficient for normoxic cells.")),
          d_alphaH(initData(&d_alphaH, "alphaH","Oxygen uptake coefficient for hypoxic cells.")),
          d_O2SaturationConstant(initData(&d_O2SaturationConstant, "O2SaturationConstant","Oxygen saturation constant for tissue, by default equal to 2.5 mmHg.")),
          d_tetraUptakeCoefficient(initData(&d_tetraUptakeCoefficient, "tetraUptakeCoefficient","Oxygen uptake coefficient for each tetra."))
{
}
// --------------------------------------------------------------------------------------
// --- destructor
// --------------------------------------------------------------------------------------
template<class DataTypes>
TetrahedronO2DiffusionFEMForceField<DataTypes>::~TetrahedronO2DiffusionFEMForceField() = default;

template< class DataTypes>
void TetrahedronO2DiffusionFEMForceField<DataTypes>::computeEdgeDiffusionCoefficient()
{
    this->edgeDiffusionCoefficient.clear();
    this->edgeUptakeCoefficient.clear();
    this->edgeDiffusionCoefficient.resize(this->nbEdges);
    this->edgeUptakeCoefficient.resize(this->nbEdges);

    /* TODO: Here is where the uptake coefficients should be updated. It may be better to 'copy' this function and
             make a version just for the oxygen uptake edge diffusion */

    sofa::Size nbTetra;
    sofa::Index i;
    sofa::Index j,k,l;
    typename DataTypes::Real val1, val2, volume;
    typename DataTypes::Real diff, diff2;
    Vec3 point[4],shapeVector[4];

    for (i=0; i<this->nbEdges; ++i)
    {
        this->edgeDiffusionCoefficient[i] = 0;
        this->edgeUptakeCoefficient[i] = 0;
    }
    nbTetra = this->m_topology->getNbTetrahedra();
    const typename TetrahedronDiffusionFEMForceField<DataTypes>::MechanicalTypes::VecCoord position = this->mechanicalObject->read(core::ConstVecCoordId::position())->getValue();

    typename DataTypes::Real anisotropyRatio = this->d_transverseAnisotropyRatio.getValue();
    bool isotropicDiffusion = (anisotropyRatio == 1.0);

    for (i=0; i<nbTetra; ++i)
    {
        // get a reference on the edge set of the ith added tetrahedron
        const EdgesInTetrahedron &te= this->m_topology->getEdgesInTetrahedron(i);
        //get a reference on the vertex set of the ith added tetrahedron
        const Tetrahedron &t= this->m_topology->getTetrahedron(i);

        // store points
        for(j=0; j<4; ++j)
            point[j]= position[t[j]];

        // compute 6 times the rest volume
        volume = dot(cross(point[1]-point[0], point[2]-point[0]), point[0]-point[3]);

        // store shape vectors
        for(j=0;j<4;++j)
        {
            if ((j%2)==0)
                shapeVector[j] = -cross(point[(j+2)%4] - point[(j+1)%4],point[(j+3)%4] - point[(j+1)%4])/volume;
            else
                shapeVector[j] = cross(point[(j+2)%4] - point[(j+1)%4],point[(j+3)%4] - point[(j+1)%4])/volume;
        }

        diff=(this->d_tetraDiffusionCoefficient.getValue())[i]*fabs(volume)/6;
        diff2=(this->d_tetraUptakeCoefficient.getValue())[i]*fabs(volume)/6;

        // isotropic case
        if (isotropicDiffusion)
        {
            for(j=0;j<6;++j)
            {
                /// local indices of the edge
                k = this->m_topology->getLocalEdgesInTetrahedron(j)[0];
                l = this->m_topology->getLocalEdgesInTetrahedron(j)[1];
                auto dot_product = dot(shapeVector[k],shapeVector[l]);
                val1 = dot_product*diff;
                this->edgeDiffusionCoefficient[te[j]] += val1;
                val2 = dot_product*diff2;
                this->edgeUptakeCoefficient[te[j]] += val2;
            }
        }
            // anisotropic case
        else
        {
            Vec3 direction = this->d_transverseAnisotropyDirectionArray.getValue()[i];
            direction.norm();

            for(j=0;j<6;++j)
            {
                /// local indices of the edge
                k = this->m_topology->getLocalEdgesInTetrahedron(j)[0];
                l = this->m_topology->getLocalEdgesInTetrahedron(j)[1];

                val1= dot(shapeVector[k],shapeVector[l]+direction * ((anisotropyRatio-1)*dot(direction,shapeVector[l])))*diff;
                this->edgeDiffusionCoefficient[te[j]] += val1;
            }
        }
    }
}

template <class DataTypes>
void TetrahedronO2DiffusionFEMForceField<DataTypes>::addForce (const core::MechanicalParams* /*mparams*/, DataVecDeriv& dataf, const DataVecCoord& datax, const DataVecDeriv& /*v*/)
{
    helper::AdvancedTimer::stepBegin("addForceO2Diffusion");

    VecDeriv& f = *dataf.beginEdit();
    const VecCoord& x = datax.getValue();

    sofa::Index v0,v1;

    Coord dp, dp2, x_average;

    for(sofa::Index i=0; i<this->nbEdges; i++ )
    {
        v0=this->m_topology->getEdge(i)[0];
        v1=this->m_topology->getEdge(i)[1];

        //Case 1D Diffusion
        if (this->d_1DDiffusion.getValue())
        {
            dp[0] = (x[v1][0]-x[v0][0]) * this->edgeDiffusionCoefficient[i];

            f[v1][0]+=dp[0];
            f[v0][0]-=dp[0];
        }
            //Case 3D Diffusion
        else
        {
            dp = (x[v1]-x[v0]) * this->edgeDiffusionCoefficient[i];  // calculate normal diffusion: diffusion_coef*delta_T
            x_average = (x[v1] + x[v0]) / 2;    // average the values on each vertex
            Coord o2_coef;
            o2_coef[0] = d_O2SaturationConstant.getValue();
            dp2 = (edgeUptakeCoefficient[i]/(o2_coef[0] + x_average[0])) * x_average[0];  // calculate amount to remove from each vertex:
                                                                                          // (uptake_coefficient/(O2_saturation - average_current_oxygen))*average_current_oxygen
            f[v1]+=dp;  // add normal diffusion vertex 1
            f[v1]-=dp2; // subtract uptake vertex 1
            f[v0]-=dp;  // subtract normal diffusion vertex 0
            f[v0]-=dp2; // subtract uptake vertex 0
        }
    }

    dataf.endEdit() ;
    sofa::helper::AdvancedTimer::stepEnd("addForceO2Diffusion");
}




}  // namespace sofa::component::forcefield