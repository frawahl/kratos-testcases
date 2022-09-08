from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

import sys
import time

class FluidDynamicsAnalysisDummyForInlet(FluidDynamicsAnalysis):

    def __init__(self,model,project_parameters,flush_frequency=10.0):
        super(FluidDynamicsAnalysisDummyForInlet,self).__init__(model,project_parameters)
        self.flush_frequency = flush_frequency
        self.last_flush = time.time()

    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """

        while self.KeepAdvancingSolutionLoop():
            self.time = self._GetSolver().AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            # self._GetSolver().Predict()
            # is_converged = self._GetSolver().SolveSolutionStep()
            # self.__CheckIfSolveSolutionStepReturnsAValue(is_converged)
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

    def FinalizeSolutionStep(self):
        super(FluidDynamicsAnalysisDummyForInlet,self).FinalizeSolutionStep()

        if self.parallel_type == "OpenMP":
            now = time.time()
            if now - self.last_flush > self.flush_frequency:
                sys.stdout.flush()
                self.last_flush = now

if __name__ == "__main__":

    with open("ProjectParametersCustomOpenMP.json",'r') as parameter_file:
    #with open("ProjectParametersCustomOpenMPRestart.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    print('Before:')
    print(parameters)

    parameters.RemoveValue('output_processes')
    parameters['processes'].RemoveValue('initial_conditions_process_list')
    parameters['processes'].RemoveValue('gravity')

    for bc_cond_proc in parameters['processes']['boundary_conditions_process_list']:
        #Parameters Object name to string and identify python_module name
        val = str(bc_cond_proc['python_module']).split(' ')[-1]
        # get rid of "
        val = val.split('"')[-2]
        
        if val == "impose_windgen_inlet_process":
            inlet_bc_cond_proc = bc_cond_proc.Clone()

    parameters['processes'].RemoveValue('boundary_conditions_process_list')
    parameters['processes'].AddEmptyArray('boundary_conditions_process_list')#,KratosMultiphysics.Parameters('''[]'''))
    parameters['processes']['boundary_conditions_process_list'].Append(inlet_bc_cond_proc)

    parameters['processes'].RemoveValue('auxiliar_process_list')
    new_output_settings = '''
    [{
            "python_module": "line_output_process",
            "kratos_module": "KratosMultiphysics",
            "process_name": "LineOutputProcess",
            "Parameters": {
                "start_point": [-449.999,0.0,0.001],
                "end_point": [-449.999,0.0,225],
                "sampling_points": 56,
                "model_part_name": "FluidModelPart.Parts_fluid",
                "output_file_settings": {
                    "file_name": "line25HupD",
                    "folder_name": "results/ascii_output/25HupD"
                },
                "output_variables": ["VELOCITY_X", "VELOCITY_Y", "VELOCITY_Z" ]
            }
        }]
    '''
    parameters['processes'].AddValue('auxiliar_process_list',KratosMultiphysics.Parameters(new_output_settings))
    
    print('AFTER:')
    print(parameters)

    model = KratosMultiphysics.Model()
    simulation = FluidDynamicsAnalysisDummyForInlet(model,parameters)
    simulation.Run()
