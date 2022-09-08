from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

import sys
import time
import math
import datetime

class FluidDynamicsAnalysisWithFlush(FluidDynamicsAnalysis):

    def __init__(self,model,project_parameters,min_dist_pos,min_dist_neg,flush_frequency=10.0):
        super(FluidDynamicsAnalysisWithFlush,self).__init__(model,project_parameters)
        self.min_dist_pos = min_dist_pos
        self.min_dist_neg = min_dist_neg
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
        #cylinder_model_part = self.model.CreateModelPart('CylinderModelPart')
        #cylinder_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        #cylinder_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        #KratosMultiphysics.ModelPartIO('circle', KratosMultiphysics.ModelPartIO.READ | KratosMultiphysics.ModelPartIO.SKIP_TIMER).ReadModelPart(cylinder_model_part)
        
        #KratosMultiphysics.CalculateDistanceToSkinProcess2D(
        #    self._GetSolver().GetComputingModelPart(),
        #    cylinder_model_part).Execute()
    

    def ModifyAfterSolverInitialize(self):
        super(FluidDynamicsAnalysisWithFlush,self).ModifyAfterSolverInitialize()
        
        # Set the continuous distance field
        r = 0.1 - 1e-10
        c = [0.0,0.0]
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            d = ((node.X - c[0])**2 + (node.Y - c[1])**2)**0.5 - r
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, d)
            if d > 0:
                if d < self.min_dist_pos:
                   self.min_dist_pos = d
            else:
                if d > self.min_dist_neg:
                    self.min_dist_neg = d

        # Set the elemental distance field
        for element in self._GetSolver().GetComputingModelPart().Elements:
            i = 0
            elem_dist = KratosMultiphysics.Vector(3)
            for node in element.GetNodes():
                elem_dist[i] = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                i += 1
            element.SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES, elem_dist)



if __name__ == "__main__":

    time_start = datetime.datetime.now().replace(microsecond=0)
    
    with open("ProjectParameters.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = FluidDynamicsAnalysisWithFlush(model,parameters, 10.0, -10.0)
    simulation.Run()
    
    # print duration of simulation
    time_end = datetime.datetime.now().replace(microsecond=0)
    delta = time_end - time_start
    hours, remainder = divmod(delta.seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    print("\n DURATION: " + str(int(hours)) + " h, " + str(int(minutes)) + " min, " + str(int(seconds)) + " sec\n")
    
    # print minimal distances
    print("\n POSITIVE AND NEGATIVE MINIMAL DISTANCES: " + str(simulation.min_dist_neg) + ", " + str(simulation.min_dist_pos) + "\n")
    
