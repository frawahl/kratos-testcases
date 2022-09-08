from math import fabs
import h5py
import numpy as np
import math
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

class WindCharacteristics:
    def __init__(self, WRF_h5_file, delta_time, f_velocity_profile, interpolation_max_height):
        self.__file = h5py.File(WRF_h5_file, "r")

        self.__f_velocity_profile = f_velocity_profile
        self.__interpolation_max_height = interpolation_max_height
        self.__y = self.__file["z"]
        self.__u = self.__file["u"]
        self.__v = self.__file["v"]
        self.__w = self.__file["w"]
        self.__latitudes = self.__file["lat"]
        self.__longitudes = self.__file["lon"]
        self.__topography = self.__file["hgt"]

        self.__TranslateGlobalCoordinatesToLocalCoordinates()
        self.__time_steps = np.arange(0, self.__y.shape[0]) * delta_time

    def Calculate(self, x, t):
        h, v, h_raw, v_raw = self.__InterpolatedWindData(x, t)
        h = h[h <= self.__interpolation_max_height]
        v = v[:len(h)]
        popt, _ = curve_fit(self.__f_velocity_profile, h, v)

        return popt, h, v, h_raw, v_raw

    def __InterpolatedWindData(self, x, t):
        time_index, time_ratio = WindCharacteristics.__InterpolatationData(t, self.__time_steps)
        x_index, x_ratio = WindCharacteristics.__InterpolatationData(x, self.__x)

        # time linear interpolation
        y_t = WindCharacteristics.__InterpolateValue(self.__y, time_index, time_ratio)
        u_t = WindCharacteristics.__InterpolateValue(self.__u, time_index, time_ratio)
        v_t = WindCharacteristics.__InterpolateValue(self.__v, time_index, time_ratio)
        w_t = WindCharacteristics.__InterpolateValue(self.__w, time_index, time_ratio)

        def get_velocity_and_h(ratio):
            y_0 = WindCharacteristics.__InterpolateValue(self.__topography, x_index, ratio)
            h = WindCharacteristics.__InterpolateValue(y_t, x_index, ratio) - y_0
            u = WindCharacteristics.__InterpolateValue(u_t, x_index, ratio)
            v = WindCharacteristics.__InterpolateValue(v_t, x_index, ratio)
            w = WindCharacteristics.__InterpolateValue(w_t, x_index, ratio)
            velocity = np.sqrt(u**2 + v**2 + w**2)

            return h, velocity

        # space linear interpolation
        h_0, velocity_0 = get_velocity_and_h(0.0)
        h_1, velocity_1 = get_velocity_and_h(1.0)
        h, velocity = WindCharacteristics.__GenerateInterpolatedCurve(h_0, velocity_0, h_1, velocity_1, x_ratio)

        return h, velocity, [h_0, h_1], [velocity_0, velocity_1]

    def __TranslateGlobalCoordinatesToLocalCoordinates(self):
        """Translates longitudes and latitudes to local x, y coordinates
        """

        latitudes = np.array(self.__latitudes) * math.pi / 180.0
        latitudes_1 = latitudes[:-1]
        latitudes_2 = latitudes[1:]
        longitude_diff = np.diff(np.array(self.__longitudes) * math.pi / 180.0)

        # convert longitudes and latitudes to x,y coordinates
        # earth radius
        r = 6.371e+6
        self.__x = r * 2*np.arcsin(np.sqrt((np.sin((latitudes_1-latitudes_2)/2))**2+np.cos(latitudes_1)*np.cos(latitudes_2)*(np.sin((longitude_diff)/2))**2))
        self.__x = np.insert(self.__x, 0, 0.0)
        self.__x = np.cumsum(self.__x)

    @staticmethod
    def __InterpolatationData(value, vector):
        index = int(np.where(vector == vector[vector <= value][-1])[0][0])
        if (index >= vector.shape[0]):
            raise Exception("Extrapolation required which is not supported.")
        ratio = (value - vector[index]) / (vector[index + 1] - vector[index])

        return index, ratio

    @staticmethod
    def __InterpolateValue(vector, index, ratio):
        return vector[index] + ratio * (vector[index + 1] - vector[index])

    @staticmethod
    def __GenerateInterpolatedCurve(x_0, y_0, x_1, y_1, ratio):
        n = x_0.shape[0]
        f_0_interpolated = interp1d(x_0, y_0)
        f_1_interpolated = interp1d(x_1, y_1)
        x_max = min(max(x_0), max(x_1))
        x_min = max(min(x_0), min(x_1))

        x = np.linspace(x_min, x_max, n)
        y_0_new = f_0_interpolated(x)
        y_1_new = f_1_interpolated(x)

        y = y_0_new + ratio * (y_1_new - y_0_new)

        return x, y

if __name__=="__main__":
    import matplotlib.pyplot as plt

    def WindLogarithmicVelocityProfile(z, a, z_0):
        return a * np.log((z + z_0) / z_0 )

    def WindPolynomialVelocityProfile(z, a, b, c, d, e, f):
        return a + b * z + c * z**2 + d * z**3 + e * z**4 + f * z**5

    # box size
    width = 3000.0
    height = 1000.0 # 1km
    width_dx = 1.0
    wrf_time = 5000.0
    velocity_profile = WindLogarithmicVelocityProfile
    n_dx = int(width / width_dx)

    wind_characteristics = WindCharacteristics("Testfile_Daniel.h5", 0.8, velocity_profile, height)
    x_space = np.linspace(0.0, width, n_dx)
    y_space = np.linspace(0.0, height, n_dx)
    v_space = np.zeros([n_dx, n_dx])
    for i, x in enumerate(x_space):
        coeffs, _, _, _, _ = wind_characteristics.Calculate(x, wrf_time)
        v_space[i, :] = velocity_profile(y_space, *coeffs)

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(x_space, y_space, v_space,cmap='viridis', edgecolor='none')
    ax.set_xlabel("width [m]")
    ax.set_ylabel("height [m]")
    ax.set_zlabel("velocity [m/s]")
    plt.grid(True)
    # plt.legend(loc = "lower right")
    plt.show()
