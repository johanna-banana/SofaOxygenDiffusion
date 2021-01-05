//
// Created by Somers on 12/23/2020.
//

#ifndef TESTSOFAEXECUTABLE_SOFAOXYGENDIFFUSION_H
#define TESTSOFAEXECUTABLE_SOFAOXYGENDIFFUSION_H

#ifdef MAKE_MEX_FILE
#include "mex.h"
#include "class_handle.hpp"
#endif

#include <sofa/simulation/Simulation.h>
#include <SofaSimulationGraph/SimpleApi.h>
#include <SofaSimulationGraph/DAGNode.h>
using DAGNode = sofa::simulation::graph::DAGNode;

// The class that we are interfacing to
class SofaOxygenDiffusion
{
public:
    SofaOxygenDiffusion() {}
    ~SofaOxygenDiffusion();
    void run_simulation();
    void init_simulation();
private:
    sofa::simulation::Simulation::SPtr m_simulation;
    sofa::simulation::Node::SPtr m_root_node;
};


#ifdef MAKE_MEX_FILE
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Get the command string
    char cmd[64];
    if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
        mexErrMsgTxt("First input should be a command string less than 64 characters long.");

    // New
    if (!strcmp("new", cmd)) {
        // Check parameters
        if (nlhs != 1)
            mexErrMsgTxt("New: One output expected.");
        // Return a handle to a new C++ instance
        plhs[0] = convertPtr2Mat<SofaOxygenDiffusion>(new SofaOxygenDiffusion);
        return;
    }

    // Check there is a second input, which should be the class instance handle
    if (nrhs < 2)
        mexErrMsgTxt("Second input should be a class instance handle.");

    // Delete
    if (!strcmp("delete", cmd)) {
        // Destroy the C++ object
        destroyObject<SofaOxygenDiffusion>(prhs[1]);
        // Warn if other commands were ignored
        if (nlhs != 0 || nrhs != 2)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        return;
    }

    // Get the class instance pointer from the second input
    SofaOxygenDiffusion* example_instance = convertMat2Ptr<SofaOxygenDiffusion>(prhs[1]);

    // Call the various class methods
    // Train
    if (!strcmp("init_simulation", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("Train: Unexpected arguments.");
        // Call the method
        example_instance->init_simulation();
        return;
    }
    // Test
    if (!strcmp("run_simulation", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("Test: Unexpected arguments.");
        // Call the method
        example_instance->run_simulation();
        return;
    }

    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}
#endif

#endif //TESTSOFAEXECUTABLE_SOFAOXYGENDIFFUSION_H
