/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU General Public License as published by the Free  *
* Software Foundation; either version 2 of the License, or (at your option)   *
* any later version.                                                          *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    *
* more details.                                                               *
*                                                                             *
* You should have received a copy of the GNU General Public License along     *
* with this program. If not, see <http://www.gnu.org/licenses/>.              *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include "mex.h"
#include "class_handle.hpp"
#include "SofaOxygenDiffusion.h"


#include <sofa/defaulttype/VecTypes.h>

#include <SofaGraphComponent/Gravity.h>
//#include <SofaLoader/MeshGmshLoader.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <SofaBaseMechanics/UniformMass.h>
#include <SofaBoundaryCondition/ConstantForceField.h>
#include <SofaBoundaryCondition/FixedConstraint.h>
//#include <SofaSimpleFem/TetrahedronDiffusionFEMForceField.h>
 //solvers
#include <SofaBaseLinearSolver/CGLinearSolver.h>

#include <sofa/core/objectmodel/Context.h>
#include <sofa/helper/system/FileRepository.h>

#include <SofaSimulationTree/init.h>
#include <SofaSimulationTree/GNode.h>
#include <SofaSimulationTree/TreeSimulation.h>

#include <iostream>
#include <fstream>

#include <SofaBase/initSofaBase.h>

using namespace sofa::simulation::tree;
using namespace sofa::defaulttype;

using sofa::simulation::Node;
using sofa::helper::system::DataRepository;
#include <sofa/helper/system/FileSystem.h>
using sofa::helper::system::FileSystem;

using sofa::component::linearsolver::CGLinearSolver;
using sofa::component::linearsolver::GraphScatteredMatrix;
using sofa::component::linearsolver::GraphScatteredVector;

// mechanical object
using sofa::component::container::MechanicalObject;
using sofa::defaulttype::StdVectorTypes;
using sofa::defaulttype::Vec;
//using sofa::core::loader::MeshLoader;

using sofa::core::Mapping;
using sofa::core::behavior::MechanicalState;
using sofa::core::State;
using sofa::core::objectmodel::New;

using sofa::component::projectiveconstraintset::FixedConstraint;

void SofaOxygenDiffusion::init_simulation()
{
//    glutInit(&argc,argv);
    sofa::simulation::tree::init();
    sofa::component::initSofaBase();
//    sofa::component::initSofaCommon();
//    sofa::component::initSofaGeneral();
//    sofa::component::initSofaMisc();
//    sofa::helper::parse("This is a SOFA application.")
//            (argc,argv);
//    sofa::component::init();
//    sofa::gui::initMain();
//    sofa::gui::GUIManager::Init(argv[0]);
//
//    // The graph root node : gravity already exists in a GNode by default
    m_root_node = New<GNode>();
    m_root_node->setName( "root" );
    m_root_node->setGravity( Vec3dTypes::Coord(0,0,0) );
    m_root_node->setDt(0.02);

    /*
     * Sub nodes: DRE
     */
    GNode::SPtr dreNode = New<GNode>();
    dreNode->setName("DRE");


    GNode::SPtr cylNode = New<GNode>();
    cylNode->setName("Cylinder");

    // solvers
    typedef CGLinearSolver<GraphScatteredMatrix, GraphScatteredVector> CGLinearSolverGraph;
    CGLinearSolverGraph::SPtr cgLinearSolver = New<CGLinearSolverGraph>();
    cgLinearSolver->setName("cgLinearSolver");
    cgLinearSolver->f_maxIter.setValue(25);
    cgLinearSolver->f_tolerance.setValue(1.0e-9);
    cgLinearSolver->f_smallDenominatorThreshold.setValue(1.0e-9);

    // mechanical object
    typedef MechanicalObject< Vec3dTypes > MechanicalObject3d;
    MechanicalObject3d::SPtr mechanicalObject = New<MechanicalObject3d>();
    mechanicalObject->setTranslation(0,0,0);
    mechanicalObject->setRotation(0,0,0);
    mechanicalObject->setScale(1,1,1);

    // diffusion forcefield
//    typedef sofa::component::forcefield::TetrahedronDiffusionFEMForceField<Vec3dTypes> TetraDiffFEM;
//    TetraDiffFEM::SPtr tetra_diffusion_fem = New<TetraDiffFEM>();
//    tetra_diffusion_fem->setName("diffusion_fem");
//    tetra_diffusion_fem->setDiffusionCoefficient(0.6);

//    // fixed constraint
    typedef FixedConstraint< StdVectorTypes<Vec<3,double>,Vec<3,double>,double> > FixedConstraint3d;
    FixedConstraint3d::SPtr fixedConstraints = New<FixedConstraint3d>();
//    fixedConstraints->setName("Box Constraints");
//    fixedConstraints->addConstraint(0);
//    fixedConstraints->addConstraint(1); fixedConstraints->addConstraint(2); fixedConstraints->addConstraint(6); fixedConstraints->addConstraint(12); fixedConstraints->addConstraint(17); fixedConstraints->addConstraint(21); fixedConstraints->addConstraint(22);
//    fixedConstraints->addConstraint(24); fixedConstraints->addConstraint(25); fixedConstraints->addConstraint(26); fixedConstraints->addConstraint(30); fixedConstraints->addConstraint(36); fixedConstraints->addConstraint(41); fixedConstraints->addConstraint(46); fixedConstraints->addConstraint(47);
//    fixedConstraints->addConstraint(50); fixedConstraints->addConstraint(51); fixedConstraints->addConstraint(52); fixedConstraints->addConstraint(56); fixedConstraints->addConstraint(62); fixedConstraints->addConstraint(68); fixedConstraints->addConstraint(73); fixedConstraints->addConstraint(74);
//    fixedConstraints->addConstraint(77); fixedConstraints->addConstraint(78); fixedConstraints->addConstraint(79); fixedConstraints->addConstraint(83); fixedConstraints->addConstraint(89); fixedConstraints->addConstraint(95); fixedConstraints->addConstraint(100); fixedConstraints->addConstraint(101);
//    fixedConstraints->addConstraint(104); fixedConstraints->addConstraint(105); fixedConstraints->addConstraint(106); fixedConstraints->addConstraint(110); fixedConstraints->addConstraint(116); fixedConstraints->addConstraint(122); fixedConstraints->addConstraint(127); fixedConstraints->addConstraint(128);
//    fixedConstraints->addConstraint(131); fixedConstraints->addConstraint(132); fixedConstraints->addConstraint(133); fixedConstraints->addConstraint(137); fixedConstraints->addConstraint(143); fixedConstraints->addConstraint(149); fixedConstraints->addConstraint(154); fixedConstraints->addConstraint(155);
//    fixedConstraints->addConstraint(158); fixedConstraints->addConstraint(159); fixedConstraints->addConstraint(160); fixedConstraints->addConstraint(164); fixedConstraints->addConstraint(170); fixedConstraints->addConstraint(175); fixedConstraints->addConstraint(180); fixedConstraints->addConstraint(181);
//    fixedConstraints->addConstraint(184); fixedConstraints->addConstraint(185); fixedConstraints->addConstraint(186); fixedConstraints->addConstraint(190); fixedConstraints->addConstraint(196); fixedConstraints->addConstraint(201); fixedConstraints->addConstraint(206); fixedConstraints->addConstraint(205);

    m_root_node->addChild(dreNode);

    // Init the scene
    sofa::simulation::tree::getSimulation()->init(m_root_node.get());
    m_simulation = getSimulation();
//    //=======================================
//    // Run the main loop
//    sofa::gui::GUIManager::MainLoop(groot);

//
//    sofa::simulation::tree::cleanup();
}

void SofaOxygenDiffusion::run_simulation()
{
    for(uint16_t i = 0; i<5;i++)
    {
        m_simulation->animate(m_root_node.get());
        std::cout << std::to_string(i) << std::endl;
    }
}