'''An inlet boundary condition process for KratosMultiphysics

license: license.txt
'''

'''
Example of usage for user-defined exponential profile:
the variable umean is initially taken from the *.h5 file
thereafter overwritten by this provided value
the standard value for z_ref should be 10 m
{
    "python_module" : "windgen_inlet_process",
    "Parameters"    : {
        "inlet_model_part_name" : "FluidModelPart.ModelPartName",
        "ramp_time"   : 30.0,
        "write_to_file": true,
        "print_to_screen": true,
        "output_file_settings": {
            "file_name": "inlet_position",
            "folder_name": "results/ascii_output"
        },
        "y0"              : -450.0,
        "z0"              : 0.0,
        "wind_filename"   : "<filename>.h5'
        "sc_fctr_mean"    : 1.01,
        "sc_fctr_fluct"   : 0.96,
        "profile_type"    : "ud_exponential",
        "umean"           : 40.0,
        "z_ref"           : 180,
        "gamma"           : 0.23
    }
}
'''

'''
Example of usage for EC-based exponential profile:
the variable umean is taken from the *.h5 file
the velicity profile is generated according to DIN EN 1991-1-4/NA 2012-12 Tab. NA.B.2
the standard value for z_ref should be 10 m
{
    "python_module" : "windgen_inlet_process",
    "Parameters"    : {
        "inlet_model_part_name" : "FluidModelPart.ModelPartName",
        "ramp_time"   : 30.0,
        "write_to_file": true,
        "print_to_screen": true,
        "output_file_settings": {
            "file_name": "inlet_position",
            "folder_name": "results/ascii_output"
        },
        "y0"              : -450.0,
        "z0"              : 0.0,
        "wind_filename"   : "<filename>.h5'
        "sc_fctr_mean"    : 1.01,
        "sc_fctr_fluct"   : 0.96,
        "profile_type"    : "ec_exponential",
        "z_ref"           : 290.0,
        "terrain_category": "II"
    }
}
'''

'''
Example of usage for logarithmic profile:
the variables: log_z0, z and umean are taken from the *.h5 file
{
    "python_module" : "windgen_inlet_process",
    "Parameters"    : {
        "inlet_model_part_name" : "FluidModelPart.ModelPartName",
        "ramp_time"   : 30.0,
        "write_to_file": true,
        "print_to_screen": true,
        "output_file_settings": {
            "file_name": "inlet_position",
            "folder_name": "results/ascii_output"
        },
        "y0"              : -450.0,
        "z0"              : 0.0,
        "wind_filename"   : "<filename>.h5"
        "sc_fctr_mean"    : 1.01,
        "sc_fctr_fluct"   : 0.96,
        "profile_type"    : "logarithmic"
    }
}
'''

__all__ = ['Factory', 'ImposeWindGenInletProcess']


import h5py
import os

from collections import namedtuple, Mapping
from math import isclose, floor, ceil, log

import KratosMultiphysics as Kratos
from KratosMultiphysics import VELOCITY_X, VELOCITY_Y, VELOCITY_Z
from KratosMultiphysics import TIME, DELTA_TIME
from KratosMultiphysics import Model, Logger
from wind_characteristics import WindCharacteristics
# other imports
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility
from KratosMultiphysics import Parameters as KM_Parameters
from KratosMultiphysics import IS_RESTARTED as KM_IS_RESTARTED

class Parameters(Mapping):

    def __init__(self, kratos_parameters):
        self._kratos_parameters = kratos_parameters

    def __getitem__(self, key):
        param = self._kratos_parameters[key]
        if param.IsDouble():
            value = param.GetDouble()
        elif param.IsInt():
            value = param.GetInt()
        elif param.IsString():
            value = param.GetString()
        elif param.IsBool():
            value = param.GetBool()
        else:
            value = param
        return value

    def __iter__(self):
        yield from self._kratos_parameters.keys()

    def __len__(self):
        return self._kratos_parameters.size()


def Factory(settings, Model):
    return ImposeWindGenInletProcess(Model, Parameters(settings['Parameters']))


Extent = namedtuple('Extent', ['lower', 'upper'])


class RegularGrid1D:

    @property
    def lower_bound(self):
        return self.extent.lower

    @property
    def upper_bound(self):
        return self.extent.upper

    def __init__(self, start_pos, length, size):
        self.extent = Extent(start_pos, start_pos+length)
        self.size = size
        self.step_size = (self.upper_bound-self.lower_bound) / (self.size-1)

    def __getitem__(self, index):
        return self.lower_bound + self.step_size*index

    def __len__(self):
        return self.size

    def floor_index(self, coord):
        if isclose(coord, self.lower_bound, abs_tol=1e-4):
            # Guard against negative index due to floating point representation.
            return 0
        local_coord = coord - self.lower_bound
        return int(floor(local_coord / self.step_size))

    def ceil_index(self, coord):
        if isclose(coord, self.upper_bound, abs_tol=1e-4):
            return len(self) - 1
        local_coord = coord - self.lower_bound
        return int(ceil(local_coord / self.step_size))

    def has(self, coord):
        return (self.floor_index(coord) >= 0
                and self.ceil_index(coord) < len(self))


IndexSpan = namedtuple('IndexSpan', ['begin', 'end'])


class SubGrid1D:

    @property
    def lower_bound(self):
        return self.grid[self.span.begin]

    @property
    def upper_bound(self):
        return self.grid[self.span.end]

    @property
    def step_size(self):
        return self.grid.step_size

    @property
    def gslice(self):
        return slice(self.span.begin, self.span.end+1)

    def __init__(self, grid, start_pos, end_pos):
        self.grid = grid
        start_index = grid.floor_index(start_pos)
        end_index = grid.ceil_index(end_pos)
        end_index = max(end_index, start_index+1)
        self.span = IndexSpan(start_index, end_index)

    def __getitem__(self, index):
        return self.grid[self.span.begin + index]

    def __len__(self):
        return self.span.end - self.span.begin + 1

    def has(self, coord):
        return (self.floor_index(coord) >= 0
                and self.ceil_index(coord) < len(self))

    def floor_index(self, coord):
        if isclose(coord, self.lower_bound, abs_tol=1e-4):
            return 0
        return self.grid.floor_index(coord) - self.span.begin

    def ceil_index(self, coord):
        if isclose(coord, self.upper_bound, abs_tol=1e-4):
            return len(self) - 1
        return self.grid.ceil_index(coord) - self.span.begin


Grid3D = namedtuple('Grid3D', ['x', 'y', 'z'])
Grid2D = namedtuple('Grid2D', ['x', 'z'])


class InletPanel3D:

    def __init__(self, grid, data):
        self.grid = grid
        self.data = data
        self.dy = self.grid.y.step_size
        self.dz = self.grid.z.step_size
        self.y0 = self.grid.y.lower_bound
        self.z0 = self.grid.z.lower_bound

    def update(self, pos):
        i = self.grid.x.floor_index(pos)
        tx = (pos-self.grid.x[i]) / self.grid.x.step_size
        data_0 = self.data[i, self.grid.y.gslice, self.grid.z.gslice]
        data_1 = self.data[i+1, self.grid.y.gslice, self.grid.z.gslice]
        self.cut_data = (1.0 - tx) * data_0 + tx * data_1

    def interpolate(self, node):
        # xi and eta are local mesh cell coordinates in the interval [-1, 1].
        xi = ((node.Y - self.y0) % self.dy) / self.dy
        eta = ((node.Z - self.z0) % self.dz) / self.dz
        xi = 2.0 * (xi-0.5)
        eta = 2.0 * (eta-0.5)
        # Interpolate using bilinear shape functions.
        weights = (
            0.25 * (1.0-xi) * (1.0-eta),
            0.25 * (1.0+xi) * (1.0-eta),
            0.25 * (1.0+xi) * (1.0+eta),
            0.25 * (1.0-xi) * (1.0+eta)
        )
        j0 = self.grid.y.floor_index(node.Y)
        j1 = self.grid.y.ceil_index(node.Y)
        k0 = self.grid.z.floor_index(node.Z)
        k1 = self.grid.z.ceil_index(node.Z)
        return (
            weights[0] * self.cut_data[j0, k0]
            + weights[1] * self.cut_data[j1, k0]
            + weights[2] * self.cut_data[j1, k1]
            + weights[3] * self.cut_data[j0, k1]
        )


class InletPanel2D:

    def __init__(self, grid, data):
        self.grid = grid
        self.data = data
        self.dz = self.grid.z.step_size
        self.z0 = self.grid.z.lower_bound

    def update(self, pos):
        i = self.grid.x.floor_index(pos)
        tx = (pos-self.grid.x[i]) / self.grid.x.step_size
        data_0 = self.data[i, self.grid.z.gslice]
        data_1 = self.data[i+1, self.grid.z.gslice]
        self.cut_data = (1.0 - tx) * data_0 + tx * data_1

    def interpolate(self, node):
        # xi and eta are local mesh cell coordinates in the interval [-1, 1].
        xi = ((node.Z - self.z0) % self.dz) / self.dz
        xi = 2.0 * (xi-0.5)
        # Interpolate using linear shape functions.
        weights = (
            0.5 * (1.0-xi),
            0.5 * (1.0+xi)
        )
        k = self.grid.z.floor_index(node.Z)
        return (
            weights[0] * self.cut_data[k]
            + weights[1] * self.cut_data[k+1]
        )


def weak_min(nodes, key): return key(min(nodes, key=key))


def weak_max(nodes, key): return key(max(nodes, key=key))


def get_extent(nodes, key):
    lower_bound = weak_min(nodes, key=key)
    upper_bound = weak_max(nodes, key=key)
    return Extent(lower_bound, upper_bound)

# ------------------------------------------------------------------
class LogMeanProfile:

    def __init__(self, friction_velocity, log_z0, bulk_wind_speed,
                 key=lambda node: node.Z):
        self.friction_velocity = friction_velocity
        self.log_z0 = log_z0
        self.bulk_wind_speed = bulk_wind_speed
        self.key = key

    def mean_wind_speed(self, node, von_karman_const=0.41):
        return (1/von_karman_const * self.friction_velocity
                * log((self.key(node) + self.log_z0) / self.log_z0))

    @staticmethod
    def wind_speed_profile(z, friction_velocity, log_z0):
        return (1/0.41 * friction_velocity
                * log((z + log_z0) / log_z0))

class EcExpMeanProfile:

    def __init__(self, ubasic, z, z_ref, terrain_category, key=lambda node: node.Z):
        self.key = key
        self.ubasic = ubasic
        self.bulk_wind_speed = ubasic
        self.z_ref = z_ref
        self.terrain_category = terrain_category

    def mean_wind_speed(self, node):
        available_ec_coefs = {
            'I': {'alpha': 1.18, 'beta': 0.12},
            'II': {'alpha': 1.00, 'beta': 0.16},
            'III': {'alpha': 0.89, 'beta': 0.19},
            'IV': {'alpha': 0.56, 'beta': 0.30}}

        if self.terrain_category not in available_ec_coefs.keys():
            raise Exception('Terraing category ', self.terrain_category,
                            ' not possible for terrain_category. Possible values are : ',
                            ', '.join(available_ec_coefs.keys()))

        return available_ec_coefs[self.terrain_category]['alpha'] * self.ubasic * (self.key(node)/self.z_ref)**available_ec_coefs[self.terrain_category]['alpha']


class UdExpMeanProfile:

    def __init__(self, umean, z, z_ref, gamma, key=lambda node: node.Z):
        self.key = key
        self.umean = umean
        self.bulk_wind_speed = umean
        self.z_ref = z_ref
        self.gamma = gamma

    def mean_wind_speed(self, node):
        return self.umean * (self.key(node)/self.z_ref)**self.gamma

    @staticmethod
    def wind_speed_profile(z, u_mean, z_ref, gamma):
        return u_mean * (z/z_ref)**gamma
# ------------------------------------------------------------------

class WindCharacteristicsMeanProfile:
    def __init__(
        self,
        process_info,
        WRF_h5_file,
        WRF_delta_time,
        interpolation_max_height,
        bulk_wind_speed,
        wind_velocity_profile_method,
        key_x=lambda node: node.X,
        key_z=lambda node: node.Z):

        self.__wind_characteristics = WindCharacteristics(
            WRF_h5_file,
            WRF_delta_time,
            wind_velocity_profile_method,
            interpolation_max_height)

        self.__wind_velocity_profile_method = wind_velocity_profile_method

        self.bulk_wind_speed = bulk_wind_speed
        self.__process_info = process_info
        self.__key_x = key_x
        self.__key_z = key_z

    def mean_wind_speed(self, node):
        x = self.__key_x(node)
        coeffs, _, _, _, _ = self.__wind_characteristics.Calculate(x, self.__process_info[Kratos.TIME])
        z = self.__key_z(node)
        return self.__wind_velocity_profile_method(z, *coeffs)

# ------------------------------------------------------------------


class ImposeWindGenInletProcess:

    @property
    def inlet_nodes(self):
        return self.model_part.Nodes

    def __init__(self, Model, settings):

        # setting user-defined values as attributes
        for name, value in settings.items():
            if hasattr(self, name):
                # will not be the case, here for robustness
                Logger.PrintInfo('ImposeWindGenInletProcess',
                                 'Already has attribute ', name,
                                 ' with value of: ', getattr(self, name), '\n',
                                 'not setting with value: ', value)
            else:
                Logger.PrintInfo('ImposeWindGenInletProcess',
                                 'Setting attribute ', name, ' with value: ', value)
                setattr(self, name, value)

        # set value from file
        # if attribute exists from user input above,
        # values from file will not be considered
        with self.OpenFile() as file_:
            for key, value in file_.items():
                if key in ['lx', 'ly', 'lz', 'log_z0', 'z', 'umean']:
                    if hasattr(self, key):
                        Logger.PrintInfo('ImposeWindGenInletProcess',
                                         'Already has attribute ', key,
                                         ' with value of: ', getattr(
                                             self, key), '\n',
                                         'not setting value: ', value[0])
                    else:
                        Logger.PrintInfo('ImposeWindGenInletProcess',
                                         'Setting attribute ', key, ' with value: ', value[0])
                        setattr(self, key, value[0])

        self.model_part = Model[self.model_part_name]

        if (self.model_part.GetCommunicator().MyPID() == 0):

            # writing input position to file
            if self.write_to_file and hasattr(self, 'output_file_settings'):
                file_handler_params = KM_Parameters(settings["output_file_settings"])

                file_header = self._GetFileHeader()
                self.output_file = TimeBasedAsciiFileWriterUtility(self.model_part,
                                                                   file_handler_params, file_header).file

            if self.model_part.ProcessInfo[KM_IS_RESTARTED]:
                with open(self.output_file.name, 'r') as f:
                    lines = f.read().splitlines()
                    last_line = lines[-1].split(' ')

                    last_inlet_time = float(last_line[0])
                    last_inlet_pos =  float(last_line[1])

                self.inlet_position = last_inlet_pos

                Logger.PrintInfo('ImposeWindGenInletProcess',
                    'Overwriting inlet_position: ', self.inlet_position, '\n'
                    'with: ', last_inlet_pos,
                    ' because of restart from time: ', last_inlet_time)

        if self.model_part.ProcessInfo[KM_IS_RESTARTED]:
            self.inlet_position = self.model_part.GetCommunicator().GetDataCommunicator().Broadcast(self.inlet_position, 0)

        if self.profile_type == 'logarithmic':
            self.mean_profile = self.CreateLogMeanProfile()

            Logger.PrintInfo('ImposeWindGenInletProcess',
                             'applying logarithmic profile with: ',
                             'u_mean =', self.umean,
                             ', log_z0 = ', self.log_z0)

        elif self.profile_type == 'ec_exponential':
            if not hasattr(self, 'ubasic'):
                self.ubasic = self.umean
                Logger.PrintInfo('ImposeWindGenInletProcess',
                                 'Does not have attribute ubasic', '\n',
                                 'taking it equal to umean: ', self.umean)

            self.mean_profile = self.CreateEcExpMeanProfile()
            Logger.PrintInfo('ImposeWindGenInletProcess',
                             'applying EC-based exponential profile for terrain category', self.terrain_category,
                             'with: z_ref = ', self.z_ref,
                             ', ubasic = ', self.ubasic)

        elif self.profile_type == 'ud_exponential':
            self.mean_profile = self.CreateUdExpMeanProfile()
            Logger.PrintInfo('ImposeWindGenInletProcess',
                             'applying user-defined exponential profile',
                             'with: z_ref = ', self.z_ref,
                             ', umean = ', self.umean,
                             ', gamma = ', self.gamma)

        elif self.profile_type == "wrf_logarithmic":
            self.mean_profile = self.CreateWrfBasedMeanProfile(LogMeanProfile.wind_speed_profile)
            Logger.PrintInfo('ImposeWindGenInletProcess',
                             'applying WRF based mean with LogMeanProfile')
        elif self.profile_type == "wrf_ud_exponential":
            self.mean_profile = self.CreateWrfBasedMeanProfile(UdExpMeanProfile.wind_speed_profile)
            Logger.PrintInfo('ImposeWindGenInletProcess',
                             'applying WRF based mean with UdExpMeanProfile')
        else:
            raise Exception('Value ', self.profile_type, ' not possible for profile_type.',
                            'Possible values are : ud_exponential, ec_exponential, logarithmic')

        if len(self.inlet_nodes) > 0:
            # In MPI we pin the file if the process has nodes on the inlet.
            # The mappers are only valid during the lifetime of the file.
            self.file_ = self.OpenFile()
            self.mappers = self.CreateMappers(self.file_)

    def ExecuteInitialize(self):
        for node in self.inlet_nodes:
            for var, _ in self.mappers:
                node.Fix(var)

    def ExecuteFinalize(self):
        if (self.model_part.GetCommunicator().MyPID() == 0):
            self.output_file.close()

    def ExecuteInitializeSolutionStep(self):
        self.UpdateInletPosition()

        if len(self.inlet_nodes) > 0:
            self.AssignVelocity()
            self.ApplyRamp()

    def OpenFile(self):
        return h5py.File(self.wind_filename, 'r')

    def _GetFileHeader(self):
        header = '# WindGen inlet for model part ' + self.model_part_name + '\n'
        header += '# Time InletPosition\n'
        return header

    def _PrintToScreen(self, result_msg):
        Logger.PrintInfo(
            "ImposeWindGenInletProcess", "Current time: " + result_msg)

    # ------------------------------------------------------------------
    def CreateLogMeanProfile(self, key=lambda node: node.Z, von_karman_const=0.41):
        z = self.z
        lz = self.lz
        log_z0 = self.log_z0
        umean = self.umean
        bulk_wind_speed = umean * (log(lz/log_z0) - 1.0) / log(z/log_z0)
        friction_velocity = umean * von_karman_const / \
            log((z + log_z0) / log_z0)
        return LogMeanProfile(friction_velocity, log_z0, bulk_wind_speed)

    def CreateEcExpMeanProfile(self, key=lambda node: node.Z):
        ubasic = self.ubasic
        z = self.z
        z_ref = self.z_ref
        terrain_category = self.terrain_category
        return EcExpMeanProfile(ubasic, z, z_ref, terrain_category)

    def CreateUdExpMeanProfile(self, key=lambda node: node.Z):
        umean = self.umean
        z = self.z
        z_ref = self.z_ref
        gamma = self.gamma
        return UdExpMeanProfile(umean, z, z_ref, gamma)

    def CreateWrfBasedMeanProfile(self, method, key=lambda node: node.Z):
        return WindCharacteristicsMeanProfile(self.model_part.ProcessInfo, self.wrf_h5_filename, self.wrf_delta_time, self.max_height, self.umean, method, lambda node: (node.Y + 450.0))

    # ------------------------------------------------------------------

    def CreateMappers(self, file_):
        if len(file_['u'].shape) == 2:
            mappers = self.Create2DMappers(file_)
        else:
            mappers = self.Create3DMappers(file_)
        return mappers

    def Create2DMappers(self, file_):
        nx, nz = file_['u'].shape
        x_grid = RegularGrid1D(0., self.lx, nx)
        z_grid = RegularGrid1D(self.z0, self.lz, nz)
        z_extent = get_extent(self.inlet_nodes, key=lambda node: node.Z)
        z_subgrid = SubGrid1D(z_grid, z_extent.lower, z_extent.upper)
        grid = Grid2D(x_grid, z_subgrid)
        # mappers become invalid after file_ is destructed.
        mappers = ((VELOCITY_X, InletPanel2D(grid, file_['u'])),
                   (VELOCITY_Z, InletPanel2D(grid, file_['w'])))
        return mappers

    def Create3DMappers(self, file_):
        nx, ny, nz = file_['u'].shape
        x_grid = RegularGrid1D(0., self.lx, nx)
        y_grid = RegularGrid1D(self.y0, self.ly, ny)
        z_grid = RegularGrid1D(self.z0, self.lz, nz)
        y_extent = get_extent(self.inlet_nodes, key=lambda node: node.Y)
        z_extent = get_extent(self.inlet_nodes, key=lambda node: node.Z)
        y_subgrid = SubGrid1D(y_grid, y_extent.lower, y_extent.upper)
        z_subgrid = SubGrid1D(z_grid, z_extent.lower, z_extent.upper)
        grid = Grid3D(x_grid, y_subgrid, z_subgrid)
        # mappers become invalid after file_ is destructed.
        mappers = ((VELOCITY_X, InletPanel3D(grid, file_['u'])),
                   (VELOCITY_Y, InletPanel3D(grid, file_['v'])),
                   (VELOCITY_Z, InletPanel3D(grid, file_['w'])))
        return mappers

    def UpdateInletPosition(self):
        time = self.model_part.ProcessInfo[TIME]
        dt = self.model_part.ProcessInfo[DELTA_TIME]
        self.inlet_position += dt * self.mean_profile.bulk_wind_speed

        if (self.model_part.GetCommunicator().MyPID() == 0):

            res_labels = ["time: ", "inlet_position: "]
            output_vals = [str(time), str(self.inlet_position)]

            if (self.print_to_screen):
                result_msg = 'WinGen Inlet for model part ' + \
                    self.model_part_name + '\n'
                result_msg += ', '.join([a+b for a,
                                            b in zip(res_labels, output_vals)])
                self._PrintToScreen(result_msg)

            if (self.write_to_file):
                self.output_file.write(' '.join(output_vals) + "\n")

        # Use the property that the wind data is periodic to cycle the domain.
        if self.inlet_position >= self.lx:
            self.inlet_position -= self.lx

    def AssignVelocity(self):
        for var, mapper in self.mappers:
            mapper.update(self.inlet_position)
            if var == VELOCITY_X:
                for node in self.inlet_nodes:
                    vel = self.mean_profile.mean_wind_speed(
                        node) * self.sc_fctr_mean
                    vel += mapper.interpolate(node) * self.sc_fctr_fluct
                    node.SetSolutionStepValue(var, vel)
            else:
                for node in self.inlet_nodes:
                    vel = mapper.interpolate(node) * self.sc_fctr_fluct
                    node.SetSolutionStepValue(var, vel)

    def ApplyRamp(self):
        time = self.model_part.ProcessInfo[TIME]
        if time < self.ramp_time:
            scal = time / self.ramp_time
            for node in self.inlet_nodes:
                for var, _ in self.mappers:
                    vel = node.GetSolutionStepValue(var)
                    node.SetSolutionStepValue(var, scal * vel)

    def Check(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass
