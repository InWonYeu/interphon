import numpy as np
from InterPhon import error
from InterPhon.util import MatrixLike


class ThermalProperty(object):
    """
    ThermalProperty class to analyze the thermal properties determined by phonon dispersion.
    Its instance variables contain an instance of :class:`core.PostProcess` class storing the information on eigen-frequency and eigen-mode.
    The given information is further refined to predict vibrational entropy and free energy.
    The thermal properties are written to external data and graphic files using the :class:`analysis.ThermalProperty.write`
    and :class:`analysis.ThermalProperty.plot` methods, respectively.

    :param process: Instance of PostProcess class
    :type process: :class:`core.PostProcess`
    :param temp: Temperatures to be studied, defaults to range(0, 1000, 10)
    :type temp: MatrixLike
    """
    def __init__(self, process,
                 temp: MatrixLike = range(0, 1000, 10)):
        """
        Constructor of ThermalProperty class.
        """
        self.process = process
        self.temp = np.array(temp)
        self.free_energy = np.zeros(self.temp.shape)
        self.entropy = np.zeros(self.temp.shape)
        # self.heat_capacity = np.zeros(self.temp.shape)

    def set(self):
        """
        Set the vibrational entropy and free energy in the range of temperatures.
        """
        kb = 1.38 * 10 ** (-23) / (1.602 * 10 ** (-19))
        h = 6.626 * 10 ** (-34) / (1.602 * 10 ** (-19))

        is_imaginary = False
        check_zero = np.isin(self.temp, 0).nonzero()[0]
        if check_zero.shape[0] == 0:
            for eig_freqs in self.process.w_q:
                for eig_freq in eig_freqs:
                    if eig_freq < 0:
                        is_imaginary = True
                    else:
                        eig_freq *= 10 ** 12

                        # To solve the overflow encountered in exp
                        self.free_energy += (1.0 / 2.0 * h * eig_freq
                                             + kb * self.temp * np.log(1 - np.exp(-h * eig_freq / (kb * self.temp)))) \
                                            / len(self.process.k_points) \
                                            / (self.process.w_q.shape[1] / 3)

                        self.entropy += (1 / (2 * self.temp) * h * eig_freq \
                                        * np.cosh(h * eig_freq / (2 * kb * self.temp)) \
                                        / np.sinh(h * eig_freq / (2 * kb * self.temp)) \
                                        - kb * np.log(2 * np.sinh(h * eig_freq / (2 * kb * self.temp)))) \
                                        / len(self.process.k_points) \
                                        / (self.process.w_q.shape[1] / 3)

                        # self.heat_capacity += kb * np.power(h * eig_freq / (kb * self.temp), 2) \
                        #                       * np.exp(h * eig_freq / (kb * self.temp)) \
                        #                       / np.power(np.exp(h * eig_freq / (kb * self.temp)) - 1, 2) \
                        #                       / len(self.process.k_points) \
                        #                       / (self.process.w_q.shape[1] / 3)
            if is_imaginary:
                print("\nCaution: ", error.Thermal_Imaginary_Frequency())

        elif check_zero.shape[0] == 1:
            zero_index = check_zero[0]
            tmp_temp = np.delete(self.temp, zero_index)
            tmp_free_energy = np.zeros(tmp_temp.shape)
            zero_point_energy = 0.0
            tmp_entropy = np.zeros(tmp_temp.shape)

            for eig_freqs in self.process.w_q:
                for eig_freq in eig_freqs:
                    if eig_freq < 0:
                        is_imaginary = True
                    else:
                        eig_freq *= 10 ** 12

                        # To solve the overflow encountered in exp
                        tmp_free_energy += (1.0 / 2.0 * h * eig_freq
                                            + kb * tmp_temp * np.log(1 - np.exp(-h * eig_freq / (kb * tmp_temp)))) \
                                           / len(self.process.k_points) \
                                           / (self.process.w_q.shape[1] / 3)

                        zero_point_energy += 1.0 / 2.0 * h * eig_freq \
                                             / len(self.process.k_points) / (self.process.w_q.shape[1] / 3)

                        tmp_entropy += (1 / (2 * tmp_temp) * h * eig_freq
                                        * np.cosh(h * eig_freq / (2 * kb * tmp_temp))
                                        / np.sinh(h * eig_freq / (2 * kb * tmp_temp))
                                        - kb * np.log(2 * np.sinh(h * eig_freq / (2 * kb * tmp_temp)))) \
                                       / len(self.process.k_points) \
                                       / (self.process.w_q.shape[1] / 3)

            self.free_energy = np.insert(tmp_free_energy, zero_index, zero_point_energy)
            self.entropy = np.insert(tmp_entropy, zero_index, 0.0)

            if is_imaginary:
                print("\nCaution: ", error.Thermal_Imaginary_Frequency())

    def write(self, out_folder='.'):
        """
        Write thermal properties in the file name of **thermal_properties.dat**.

        :param out_folder: Folder path for **thermal_properties.dat** to be stored, defaults to .
        :type out_folder: str
        """
        with open(out_folder + '/thermal_properties.dat', 'w') as outfile:
            comment = "Thermal Properties / atom"
            outfile.write("%s" % comment + '\n')
            outfile.write("%s" %
                          '    Temperature (K)' +
                          '    Free_energy (eV/atom)' +
                          '    Entropy (meV/K/atom)' + '\n')

            for x, y1, y2 in zip(self.temp, self.free_energy, self.entropy):
                line = ' %20.9f ' % x + ' %20.9f ' % y1 + ' %20.9f ' % (y2 * 1000)
                outfile.write("%s" % line + '\n')

    def plot(self, legend_location='best'):
        """
        Plot thermal properties in the file name of **thermal_properties.png**.

        :param legend_location: Location of DOS legend, defaults to best
        :type legend_location: str
        """
        from matplotlib import pyplot as plt

        # plt.rcParams["font.family"] = "Arial"
        # plt.rcParams["figure.figsize"] = (13, 9)
        # plt.rcParams['lines.linewidth'] = 4
        # # plt.rcParams['lines.color'] = 'r'
        # plt.rcParams['axes.grid'] = False
        # plt.rc("font", size=22)  # controls default text sizes
        # # plt.rc("axes", titlesize=22)  # fontsize of the axes title
        # plt.rc("axes", labelsize=22)  # fontsize of the x and y labels
        # plt.rc("xtick", labelsize=22)  # fontsize of the tick labels
        # plt.rc("ytick", labelsize=22)  # fontsize of the tick labels

        font_x = {'family': 'Arial', 'size': 36, 'color': 'black', 'weight': 'bold'}
        font_y = {'family': 'Arial', 'size': 36, 'color': 'black', 'weight': 'bold'}
        font_tick = {'size': 36, 'color': 'black'}
        font_legend = {'family': 'Arial', 'size': 24, 'color': 'black', 'weight': 'bold'}

        self.process.num_figure = self.process.num_figure + 1
        fig = plt.figure(self.process.num_figure, figsize=(13, 9))
        ax = fig.subplots()

        x_min = self.temp.min()
        x_max = self.temp.max()

        ax.plot(self.temp, self.free_energy, label='Free energy (eV/atom)', linewidth=4)
        ax.plot(self.temp, self.entropy * 1000, label='Entropy (meV/K/atom)', linewidth=4)
        # ax.plot(self.temp, self.heat_capacity * 1000, label='Heat capacity (meV/K/atom)', linewidth=4)

        ax.set_xlabel('Temperature (K)', fontdict=font_x)
        labels = [item for item in ax.get_xticks()]
        ax.set_xlim(x_min, x_max)
        xtick_labels = [round(float(label)) for label in labels]
        ax.set_xticks(xtick_labels)
        ax.set_xticklabels(labels=xtick_labels, fontdict=font_tick)
        ax.set_xlim(x_min, x_max)

        ax.set_ylabel('Thermal Properties', fontdict=font_y)
        labels = [item for item in ax.get_yticks()]
        ytick_labels = [round(float(label), 2) for label in labels]
        ax.set_yticklabels(labels=ytick_labels, fontdict=font_tick)

        ax.legend(loc=legend_location, fontsize='xx-large')

        fig.tight_layout()
        plt.savefig('thermal_properties.png', dpi=300, format='png', bbox_inches='tight')
        plt.show()
