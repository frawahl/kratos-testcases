from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.ConvectionDiffusionApplication

from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis

import sys
import math
import time

def _EvaluateSolutionField(x, y):
    return x**2 + y**2

def _EvaluateFluxField(x, y):
    normal = KM.Vector(2)
    normal[0] = 0.0
    normal[1] = 1.0
    du_dx = 2*x
    du_dy = 2*y
    return (normal[0]*du_dx + normal[1]*du_dy)

def _EvaluateSourceTerm(x, y):
    return -4.0

class ConvectionDiffusionAnalysisWithFlush(ConvectionDiffusionAnalysis):

    def __init__(self,model,project_parameters,flush_frequency=10.0):
        super(ConvectionDiffusionAnalysisWithFlush,self).__init__(model,project_parameters)
        self.flush_frequency = flush_frequency
        self.last_flush = time.time()

    def ModifyInitialGeometry(self):
        # Set the levelset function
        dist_y = 0.750000000001
        	#8.0/9.0 pretty #0.751
        	#0.750000000001 huge values! #0.7500000000001 Kratos::Exception!
        dist_min = 10.0      
        
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            dist = dist_y - node.Y
            node.SetSolutionStepValue(KM.DISTANCE, 0, dist)
                     
            if dist > 0.0 and dist < dist_min:
                dist_min = dist
                
        print("\n\n\n ---------   min dist = " + str(dist_min) + "  ---------\n\n\n")

    def ApplyBoundaryConditions(self):
        super().ApplyBoundaryConditions()

        # Get variables
        conv_diff_settings = self._GetSolver().GetComputingModelPart().ProcessInfo[KM.CONVECTION_DIFFUSION_SETTINGS]
        unknown_variable = conv_diff_settings.GetUnknownVariable()
        flux_variable = conv_diff_settings.GetSurfaceSourceVariable()
        source_variable = conv_diff_settings.GetVolumeSourceVariable()

        # Set BCs in active nodes
        dist_y = 0.750000000001
        tol = 0.0  # 0.2
        for node in self.model.GetModelPart("ThermalModelPart").Nodes:
            if (node.Y < 0.000001 or (node.Y < (dist_y+tol) and node.X < 0.000001) or (node.Y < (dist_y+tol) and node.X > 0.999999)):
                node.Fix(unknown_variable)
                temp = _EvaluateSolutionField(node.X, node.Y)
                node.SetSolutionStepValue(unknown_variable, 0, temp)
       
        # Set source term
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            heat_flux = _EvaluateSourceTerm(node.X, node.Y)
            node.SetSolutionStepValue(source_variable, 0, heat_flux)

        # # Set the analytical solution in the entire domain
        # for node in self._GetSolver().GetComputingModelPart().Nodes:
        #     if node.GetSolutionStepValue(KM.DISTANCE) > 0.0:
        #         node.Fix(unknown_variable)
        #         temp = _EvaluateSolutionField(node.X, node.Y)
        #         node.SetSolutionStepValue(unknown_variable, 0, temp)

        # Set Nitsche penalty coefficient gamma
        self._GetSolver().GetComputingModelPart().ProcessInfo[KM.ConvectionDiffusionApplication.PENALTY_DIRICHLET] = 1e0

    def OutputSolutionStep(self):
        # Plot the nodal error
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            exact = _EvaluateSolutionField(node.X, node.Y)
            obtained = node.GetSolutionStepValue(KM.TEMPERATURE)
            node.SetSolutionStepValue(KM.ConvectionDiffusionApplication.PROJECTED_SCALAR1, 0, exact - obtained)

        super(ConvectionDiffusionAnalysisWithFlush,self).OutputSolutionStep()

    def FinalizeSolutionStep(self):
        super(ConvectionDiffusionAnalysisWithFlush,self).FinalizeSolutionStep()

        if self.parallel_type == "OpenMP":
            now = time.time()
            if now - self.last_flush > self.flush_frequency:
                sys.stdout.flush()
                self.last_flush = now

if __name__ == "__main__":

    error_norm_vect = []
    flux_error_norm_vect = []
    h_vect = [0.2] #, 0.1, 0.05, 0.025, 0.0125, 0.00625]
    for i_mesh in range(len(h_vect)):

        with open("ProjectParameters.json",'r') as parameter_file:
            parameters = KM.Parameters(parameter_file.read())

        aux_problem_name = "rectangle_1_15_str_ref_" + "5" #str(i_mesh)
        parameters["problem_data"]["problem_name"].SetString(aux_problem_name)
        parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(aux_problem_name)
        parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString(aux_problem_name)

        model = KM.Model()
        simulation = ConvectionDiffusionAnalysisWithFlush(model,parameters)
        simulation.Run()

        # Calculate domain error norm
        settings = KM.Parameters(r'''{
            "model_part_name" : "ThermalModelPart",
            "level_set_type" : "continuous",
            "shape_functions" : "standard",
            "distance_variable_name" : "DISTANCE"
        }''')
        error_norm_process = KM.ConvectionDiffusionApplication.GaussPointErrorUtility(model, settings)
        error_norm = error_norm_process.Execute()
        error_norm_vect.append(error_norm)

    print(h_vect)
    print(error_norm_vect)
