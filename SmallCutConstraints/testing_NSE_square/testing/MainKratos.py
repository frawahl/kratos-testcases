from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

import sys
import time
import math
import datetime

class FluidDynamicsAnalysisWithFlush(FluidDynamicsAnalysis):

    def __init__(self,model,project_parameters,flush_frequency=10.0):
        super(FluidDynamicsAnalysisWithFlush,self).__init__(model,project_parameters)
        self.flush_frequency = flush_frequency
        self.last_flush = time.time()

    def FinalizeSolutionStep(self):
        super(FluidDynamicsAnalysisWithFlush,self).FinalizeSolutionStep()

        if self.parallel_type == "OpenMP":
            now = time.time()
            if now - self.last_flush > self.flush_frequency:
                sys.stdout.flush()
                self.last_flush = now

    def ModifyInitialGeometry(self):
        super(FluidDynamicsAnalysisWithFlush,self).ModifyInitialGeometry()

        # Read the cylinder geometry
        square_model_part = self.model.CreateModelPart('SquareModelPart')
        square_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        square_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        KratosMultiphysics.ModelPartIO('square', KratosMultiphysics.ModelPartIO.READ | KratosMultiphysics.ModelPartIO.SKIP_TIMER).ReadModelPart(square_model_part)
        
        KratosMultiphysics.CalculateDistanceToSkinProcess2D(
            self._GetSolver().GetComputingModelPart(),
            square_model_part).Execute()
        
        self.calculate_distances = False

    def ApplyBoundaryConditions(self):
        # Calculate the level-set function
        if self.calculate_distances:
            square_model_part = self.model.GetModelPart('SquareModelPart')
            KratosMultiphysics.CalculateDistanceToSkinProcess2D(
            self._GetSolver().GetComputingModelPart(),
            square_model_part).Execute()
            self.calculate_distances = False

        # Apply the rest of boundary conditions
        super(FluidDynamicsAnalysisWithFlush,self).ApplyBoundaryConditions()


if __name__ == "__main__":

    time_start = datetime.datetime.now().replace(microsecond=0)

    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = FluidDynamicsAnalysisWithFlush(model,parameters)
    simulation.Run()
    
    # print duration of simulation
    time_end = datetime.datetime.now().replace(microsecond=0)
    delta = time_end - time_start
    hours, remainder = divmod(delta.seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    print("\n DURATION " + str(int(hours)) + " h, " + str(int(minutes)) + " min, " + str(int(seconds)) + " sec\n")
