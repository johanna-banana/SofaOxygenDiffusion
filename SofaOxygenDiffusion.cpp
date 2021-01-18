
#include "SofaOxygenDiffusion.h"

#ifdef USE_GUI
#include <sofa/gui/GUIManager.h>
#include <sofa/gui/Main.h>
#endif


#include <sofa/defaulttype/VecTypes.h>

#include <SofaGeneralLoader/MeshGmshLoader.h>
using GmshLoader = sofa::component::loader::MeshGmshLoader;

#include <SofaBaseTopology/MeshTopology.h>
#include <SofaBaseTopology/TetrahedronSetTopologyContainer.h>
using TetraMeshTopology = sofa::component::topology::TetrahedronSetTopologyContainer;
#include <SofaBaseTopology/TetrahedronSetGeometryAlgorithms.h>
using TetraGeoAlgorithm = sofa::component::topology::TetrahedronSetGeometryAlgorithms<sofa::defaulttype::Vec3dTypes>;

#include <SofaBaseMechanics/MechanicalObject.h>
#include <SofaBoundaryCondition/ConstantForceField.h>
#include <SofaBoundaryCondition/FixedConstraint.h>
#include "TetrahedronO2DiffusionFEMForceField.h"

//solvers
#include <SofaBaseLinearSolver/CGLinearSolver.h>
using CGLinearSolver = sofa::component::linearsolver::CGLinearSolver<sofa::component::linearsolver::GraphScatteredMatrix, sofa::component::linearsolver::GraphScatteredVector>;
#include <SofaImplicitOdeSolver/EulerImplicitSolver.h>
using EulerImplicitSolver = sofa::component::odesolver::EulerImplicitSolver;

#include <sofa/core/objectmodel/Context.h>
#include <sofa/helper/system/FileRepository.h>
using sofa::helper::system::DataRepository;

#include <SofaSimulationGraph/SimpleApi.h>
#include <SofaSimulationGraph/DAGNode.h>
#include <SofaSimulation/initSofaSimulation.h>
#include <SofaSimulationGraph/init.h>
#include <SofaSimulationCommon/init.h>
using DAGNode = sofa::simulation::graph::DAGNode;
using sofa::simulation::Node;

#include <iostream>
#include <SofaBase/initSofaBase.h>

// mechanical object
using sofa::component::container::MechanicalObject;
using sofa::defaulttype::StdVectorTypes;
using vec1d = sofa::defaulttype::Vec1dTypes;
using sofa::defaulttype::Vec;

using sofa::core::behavior::MechanicalState;
using sofa::core::State;
using sofa::core::objectmodel::New;

using sofa::component::projectiveconstraintset::FixedConstraint;
using Tag = sofa::core::objectmodel::Tag;
using Data = sofa::core::objectmodel::Data<std::string>;

#ifdef USE_GUI
void SofaOxygenDiffusion::init_gui()
{

}
#endif


void SofaOxygenDiffusion::init_simulation()
{
#ifdef USE_GUI
    init_gui();
#endif
    sofa::component::initSofaBase();
    sofa::initSofaSimulation();

    m_simulation = sofa::simpleapi::createSimulation();
    m_root_node = sofa::simpleapi::createRootNode(m_simulation, "root");
    m_root_node->setName( "root" );
    m_root_node->setGravity( sofa::defaulttype::Vec3dTypes::Coord(0,0,0) );
    m_root_node->setDt(0.02);

    // solvers
    EulerImplicitSolver::SPtr eulerSolver = New<EulerImplicitSolver>();
    eulerSolver->setName("eulerSolver");
    eulerSolver->f_rayleighStiffness.setValue(0.1);
    eulerSolver->f_rayleighMass.setValue(0.1);
    eulerSolver->f_firstOrder.setValue(true);
    CGLinearSolver::SPtr linearSolver = New<CGLinearSolver>();
    linearSolver->setName("linearSolver");
    linearSolver->f_maxIter.setValue(1000);
    linearSolver->f_tolerance.setValue(1.0e-10);

    // Model
    GmshLoader::SPtr mesh = New<GmshLoader>();
    mesh->setFilename("box.msh");
    mesh->doLoad();

    // Topology
    TetraMeshTopology::SPtr topology = New<TetraMeshTopology>();
    topology->setSrc("", mesh.get());
    std::cout << "tetrahedra:" << topology->getNumberOfTetrahedra() << std::endl;
    std::cout << "nodes:" << mesh->d_positions.getValue().size() << std::endl;
    TetraGeoAlgorithm::SPtr geoAlgo = New<TetraGeoAlgorithm>();
    topology->addTag(Tag("mechanics"));
    geoAlgo->addTag(Tag("mechanics"));

    // mechanical object - positions
    typedef MechanicalObject< sofa::defaulttype::Vec3dTypes > MechanicalObject3d;
    MechanicalObject3d::SPtr mechanicalObject = New<MechanicalObject3d>();
    mechanicalObject->setSrc("", mesh.get());
    mechanicalObject->setName("node_positions");
    mechanicalObject->setTranslation(0,0,0);
    mechanicalObject->setRotation(0,0,0);
    mechanicalObject->setScale(1,1,1);
    mechanicalObject->addTag(Tag("mechanics"));

    // mechanical object - oxygen levels
    typedef MechanicalObject<vec1d> MechanicalObject1d;
    MechanicalObject1d::SPtr oxygen = New<MechanicalObject1d>();
    oxygen->resize(mesh->d_positions.getValue().size());
    oxygen->setName("oxygen_levels");
    oxygen->addTag(Tag("oxygen"));

    // diffusion forcefield
    typedef sofa::component::forcefield::TetrahedronO2DiffusionFEMForceField<sofa::defaulttype::Vec1dTypes> TetraDiffFEM;
    TetraDiffFEM::SPtr tetra_diffusion_fem = New<TetraDiffFEM>();
    tetra_diffusion_fem->setName("diffusion_fem");
    tetra_diffusion_fem->l_topology = topology;
    tetra_diffusion_fem->addTag(Tag("oxygen"));
    tetra_diffusion_fem->d_tagMeshMechanics.setValue("mechanics");

    m_root_node->addObject(eulerSolver);
    m_root_node->addObject(linearSolver);
    m_root_node->addObject(mechanicalObject);
    m_root_node->addObject(oxygen);
    m_root_node->addObject(mesh);
    m_root_node->addObject(topology);
    m_root_node->addObject(geoAlgo);
    m_root_node->addObject(tetra_diffusion_fem);

    // Init the scene
    m_simulation->init(m_root_node.get());
}

SofaOxygenDiffusion::~SofaOxygenDiffusion()
{
    sofa::simulation::graph::cleanup();
    sofa::simulation::common::cleanup();
}

void SofaOxygenDiffusion::run_simulation()
{
#ifdef USE_GUI
    sofa::gui::initMain();
    sofa::gui::GUIManager::Init("simple_scene", "qt");
    sofa::gui::GUIManager::createGUI(m_root_node, "SofaOxygenDiffusion");
    sofa::gui::GUIManager::MainLoop(m_root_node);

#endif
    for(uint16_t i = 0; i<5;i++)
    {
        m_simulation->animate(m_root_node.get());
        std::cout << "Time: " << m_root_node->getTime() << "s" << std::endl;
    }
}