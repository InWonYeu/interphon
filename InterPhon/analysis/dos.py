import numpy as np
from InterPhon.util import tetrahedron_1d, tetrahedron_2d, tetrahedron_3d


class DOS(object):
    def __init__(self, process, sigma: float = 0.1, num_dos: int = 200):
        self.process = process
        self.sigma = sigma
        self.num_dos = num_dos

        _minimum_freq = np.min(self.process.w_q)
        _maximum_freq = np.max(self.process.w_q)

        self._freq = np.arange(_minimum_freq - 2, _maximum_freq + 2,
                               (_maximum_freq - _minimum_freq + 4) / self.num_dos)
        self._pdos = np.zeros((len(self.process.unit_cell.xyz_true), self._freq.shape[0]))
        self._tdos = np.empty((self._freq.shape[0],))

    @property
    def freq(self):
        return self._freq

    @property
    def pdos(self):
        return self._pdos

    @property
    def tdos(self):
        return self._tdos

    def set(self):
        if self.sigma == 0.0:
            # Linear Tetrahedron Method for Brillouin zone integration
            _ind_pbc = self.process.user_arg.periodicity.nonzero()[0]
            if _ind_pbc.shape[0] == 0:
                self._freq = self.process.w_q[0, :]
                self._tdos = np.ones(self._freq.shape[0])

            elif _ind_pbc.shape[0] == 1:
                self._pdos = tetrahedron_1d(self._freq,
                                            self._pdos,
                                            self.process.k_points,
                                            self.process.w_q,
                                            self.process.v_q,
                                            self.process.auto_k_points,
                                            _ind_pbc)
                self._tdos = self._pdos.sum(axis=0)

            elif _ind_pbc.shape[0] == 2:
                self._pdos = tetrahedron_2d(self._freq,
                                            self._pdos,
                                            self.process.k_points,
                                            self.process.w_q,
                                            self.process.v_q,
                                            self.process.auto_k_points,
                                            _ind_pbc)
                self._tdos = self._pdos.sum(axis=0)

            elif _ind_pbc.shape[0] == 3:
                self._pdos = tetrahedron_3d(self._freq,
                                            self._pdos,
                                            self.process.k_points,
                                            self.process.w_q,
                                            self.process.v_q,
                                            self.process.auto_k_points,
                                            _ind_pbc)
                self._tdos = self._pdos.sum(axis=0)

        else:
            # Gaussian Smearing Method for Brillouin zone integration
            for eig_freqs, eig_modes in zip(self.process.w_q, self.process.v_q):
                for ind_freq, eig_freq in enumerate(eig_freqs):
                    for ind_mode, _ in enumerate(self.process.unit_cell.xyz_true):
                        self._pdos[ind_mode, :] = self._pdos[ind_mode, :] \
                                        + 1 / (self.sigma * np.sqrt(2 * np.pi)) \
                                        * np.exp(- (self._freq - eig_freq) ** 2 / (2 * self.sigma ** 2)) \
                                        / len(self.process.k_points) \
                                        * (abs(eig_modes[ind_freq, ind_mode]) ** 2)

            self._tdos = self._pdos.sum(axis=0)

    def write(self, out_folder='.'):
        with open(out_folder + '/total_dos.dat', 'w') as outfile:
            if self.sigma == 0.0:
                comment = "Total phonon DOS by Linear Tetrahedron Method"
                outfile.write("%s" % comment + '\n')

            elif self.sigma != 0.0:
                comment = "Total phonon DOS by Gaussian Smearing with Sigma = %f" % self.sigma
                outfile.write("%s" % comment + '\n')

            outfile.write("%s" %
                          '    Frequency (THz)' +
                          '    Total_Density_of_State' + '\n')

            for x, y in zip(self._freq, self._tdos):
                line = ' %16.9f ' % x + ' %16.9f ' % y
                outfile.write("%s" % line + '\n')

        with open(out_folder + '/projected_dos.dat', 'w') as outfile:
            if self.sigma == 0.0:
                comment = "Projected phonon DOS by Linear Tetrahedron Method"
                outfile.write("%s" % comment + '\n')

            elif self.sigma != 0.0:
                comment = "Projected phonon DOS by Gaussian Smearing with Sigma = %f" % self.sigma
                outfile.write("%s" % comment + '\n')

            line = "%s" % '    Frequency (THz)'
            for xyz_true in self.process.unit_cell.xyz_true:
                line += '   pdos-atom-{0:0>3}  '.format(xyz_true)
            outfile.write(line + '\n')

            for x_ind, x in enumerate(self._freq):
                line = ' %16.9f ' % x
                for ind, _ in enumerate(self.process.unit_cell.xyz_true):
                    line += ' %16.9f ' % self._pdos[ind, x_ind]
                outfile.write(line + '\n')

    def plot(self,
             plot_total=True,
             atoms=None,
             elimit=None,
             color='tab:orange',
             option='plain',
             orientation='horizontal',
             label_position='left',
             legends=None,
             legend_location='best',
             bulk_dos=None,
             bulk_option='fill',
             bulk_legend=None,
             proportion=1):

        from matplotlib import pyplot as plt

        font_x = {'family': 'Arial', 'size': 36, 'color': 'black', 'weight': 'bold'}
        font_y = {'family': 'Arial', 'size': 36, 'color': 'black', 'weight': 'bold'}
        font_tick = {'size': 36, 'color': 'black'}
        font_legend = {'family': 'Arial', 'size': 24, 'color': 'black', 'weight': 'bold'}

        if orientation == 'horizontal':
            self.process.num_figure = self.process.num_figure + 1
            fig = plt.figure(self.process.num_figure, figsize=(13, 9))
            ax = fig.subplots()

            if elimit is None:
                x_min = np.floor(self.process.w_q.min()) - 1
                x_max = np.ceil(self.process.w_q.max()) + 1
            else:
                x_min = elimit[0]
                x_max = elimit[1]

            if option == 'plain':
                if bulk_dos is not None:
                    if bulk_legend is not None:
                        label = bulk_legend
                    else:
                        label = 'bulk dos'

                    if bulk_option == 'line':
                        ax.plot(bulk_dos[:, 0], bulk_dos[:, 1] * proportion,
                                color='black', linestyle=':', linewidth=4, label=label)
                    else:
                        ax.fill(bulk_dos[:, 0], bulk_dos[:, 1] * proportion, 'black', alpha=0.3, label=label)
                        # ax.fill(bulk_dos[0][:, 0], bulk_dos[0][:, 1] * proportion, 'tab:purple', alpha=0.5, label=label[0])
                        # ax.fill(bulk_dos[1][:, 0], bulk_dos[1][:, 1] * proportion, 'black', alpha=0.5, label=label[1])

                ax.plot(self._freq, self._tdos, color=color, label='total dos', linewidth=4)
                # ax.plot(self._freq, self._tdos, color=color, label='total dos', linewidth=4)

                ax.set_xlabel('Frequency (THz)', fontdict=font_x)
                ax.set_xlim(x_min, x_max)
                labels = [item for item in ax.get_xticks()]
                xtick_labels = [round(float(label), 2) for label in labels]
                ax.set_xticks(xtick_labels)
                ax.set_xticklabels(labels=xtick_labels, fontdict=font_tick)
                ax.set_xlim(x_min, x_max)

                ax.set_ylabel('DOS (a.u.)', fontdict=font_y)
                ax.set_yticks([])
                ax.set_ylim(0, None)

                ax.legend(loc=legend_location, fontsize='xx-large')

                fig.tight_layout()
                plt.savefig('dos.png', dpi=300, format='png', bbox_inches='tight')
                plt.show()

            elif option == 'line':
                if bulk_dos is not None:
                    if bulk_legend is not None:
                        label = bulk_legend
                    else:
                        label = 'bulk dos'

                    if bulk_option == 'line':
                        ax.plot(bulk_dos[:, 0], bulk_dos[:, 1] * proportion,
                                color='black', linestyle=':', linewidth=4, label=label)
                    else:
                        ax.fill(bulk_dos[:, 0], bulk_dos[:, 1] * proportion, 'black', alpha=0.3, label=label)
                        # ax.fill(bulk_dos[0][:, 0], bulk_dos[0][:, 1] * proportion, 'tab:purple', alpha=0.5, label=label[0])
                        # ax.fill(bulk_dos[1][:, 0], bulk_dos[1][:, 1] * proportion, 'black', alpha=0.5, label=label[1])

                if plot_total:
                    ax.plot(self._freq, self._tdos,
                            color=color,
                            label='total dos',
                            linewidth=4)
                    # ax.plot(self._freq, self._tdos,
                    #         color=color,
                    #         label='total dos',
                    #         linewidth=4)

                __pdos = np.zeros((len(atoms), self._freq.shape[0]))
                ind_modes = []
                try:
                    for ind_set, atom_set in enumerate(atoms):
                        ind_modes.append([ind for ind, val in enumerate(self.process.unit_cell.xyz_true)
                                          if val in atom_set])
                        for ind_mode in ind_modes[ind_set]:
                            __pdos[ind_set, :] += self._pdos[ind_mode, :]
                except TypeError:
                    print("\nFail to identify the index of atoms='{0}'".format(atoms))
                    print("Please designate the 'atoms' parameter for which DOS will be projected\n")
                    raise

                for ind_set, _pdos_ in enumerate(__pdos):
                    if legends is None:
                        label = 'atom-{0}'.format(set(atoms[ind_set]))
                    else:
                        label = legends[ind_set]

                    ax.plot(self._freq, _pdos_,
                            linestyle='--',
                            linewidth=4,
                            label=label)

                # color_set = ['tab:purple', 'black']
                # for ind_set, _pdos_ in enumerate(__pdos):
                #     if legends is None:
                #         label = 'atom-{0}'.format(set(atoms[ind_set]))
                #     else:
                #         label = legends[ind_set]
                #
                #     ax.plot(self._freq, _pdos_,
                #             linestyle='--',
                #             linewidth=4,
                #             label=label, color=color_set[ind_set])

                ax.set_xlabel('Frequency (THz)', fontdict=font_x)
                ax.set_xlim(x_min, x_max)
                labels = [item for item in ax.get_xticks()]
                xtick_labels = [round(float(label), 2) for label in labels]
                ax.set_xticks(xtick_labels)
                ax.set_xticklabels(labels=xtick_labels, fontdict=font_tick)
                ax.set_xlim(x_min, x_max)

                ax.set_ylabel('DOS (a.u.)', fontdict=font_y)
                ax.set_yticks([])
                ax.set_ylim(0, None)

                ax.legend(loc=legend_location, fontsize='xx-large')

                fig.tight_layout()
                plt.savefig('dos.png', dpi=300, format='png', bbox_inches='tight')
                plt.show()

            elif option == 'stack':
                if bulk_dos is not None:
                    if bulk_legend is not None:
                        label = bulk_legend
                    else:
                        label = 'bulk dos'

                    if bulk_option == 'line':
                        ax.plot(bulk_dos[:, 0], bulk_dos[:, 1] * proportion,
                                color='black', linestyle=':', linewidth=4, label=label)
                        # ax.plot(bulk_dos[0][:, 0], bulk_dos[0][:, 1] * proportion,
                        #         color='tab:purple', linestyle=':', linewidth=4, label=label[0])
                        # ax.plot(bulk_dos[1][:, 0], bulk_dos[1][:, 1] * proportion,
                        #         color='black', linestyle=':', linewidth=4, label=label[1])
                    else:
                        ax.fill(bulk_dos[:, 0], bulk_dos[:, 1] * proportion, 'black', alpha=0.3, label=label)

                if plot_total:
                    ax.plot(self._freq, self._tdos,
                            color=color,
                            label='total dos',
                            linewidth=4)

                __pdos = np.zeros((len(atoms), self._freq.shape[0]))
                ind_modes = []
                try:
                    for ind_set, atom_set in enumerate(atoms):
                        ind_modes.append([ind for ind, val in enumerate(self.process.unit_cell.xyz_true)
                                          if val in atom_set])
                        for ind_mode in ind_modes[ind_set]:
                            __pdos[ind_set, :] += self._pdos[ind_mode, :]
                except TypeError:
                    print("\nFail to identify the index of atoms='{0}'".format(atoms))
                    print("Please designate the 'atoms' parameter for which DOS will be projected\n")
                    raise

                if legends is None:
                    labels = []
                    for atom in atoms:
                        labels.append('atom-{0}'.format(set(atom)))
                else:
                    labels = legends

                ax.stackplot(self._freq, __pdos,
                             linewidth=1,
                             labels=labels)

                # color_set = ['tab:purple', 'black']
                # ax.stackplot(self._freq, __pdos,
                #              linewidth=1,
                #              labels=labels, colors=color_set)

                ax.set_xlabel('Frequency (THz)', fontdict=font_x)
                ax.set_xlim(x_min, x_max)
                labels = [item for item in ax.get_xticks()]
                xtick_labels = [round(float(label), 2) for label in labels]
                ax.set_xticks(xtick_labels)
                ax.set_xticklabels(labels=xtick_labels, fontdict=font_tick)
                ax.set_xlim(x_min, x_max)

                ax.set_ylabel('DOS (a.u.)', fontdict=font_y)
                ax.set_yticks([])
                ax.set_ylim(0, None)

                ax.legend(loc=legend_location, fontsize='xx-large')

                fig.tight_layout()
                plt.savefig('dos.png', dpi=300, format='png', bbox_inches='tight')
                plt.show()

        elif orientation == 'vertical':
            self.process.num_figure = self.process.num_figure + 1
            fig = plt.figure(self.process.num_figure, (6, 9))
            ax = fig.subplots()

            if elimit is None:
                y_min = np.floor(self.process.w_q.min()) - 1
                y_max = np.ceil(self.process.w_q.max()) + 1
            else:
                y_min = elimit[0]
                y_max = elimit[1]

            if option == 'plain':
                if bulk_dos is not None:
                    if bulk_legend is not None:
                        label = bulk_legend
                    else:
                        label = 'bulk dos'

                    if bulk_option == 'line':
                        ax.plot(bulk_dos[:, 1] * proportion, bulk_dos[:, 0],
                                color='black', linestyle=':', linewidth=4, label=label)
                    else:
                        ax.fill(bulk_dos[:, 1] * proportion, bulk_dos[:, 0], 'black', alpha=0.3, label=label)

                if label_position == 'right':
                    ax.yaxis.tick_right()
                    ax.yaxis.set_ticks_position('both')
                    ax.yaxis.set_label_position("right")

                ax.plot(self._tdos, self._freq, color=color, label='total dos', linewidth=4)

                ax.set_ylabel('Frequency (THz)', fontdict=font_x)
                ax.set_ylim(y_min, y_max)
                labels = [item for item in ax.get_yticks()]
                ytick_labels = [round(float(label), 2) for label in labels]
                ax.set_yticks(ytick_labels)
                ax.set_yticklabels(labels=ytick_labels, fontdict=font_tick)
                ax.set_ylim(y_min, y_max)

                ax.set_xlabel('DOS (a.u.)', fontdict=font_y)
                ax.set_xticks([])
                ax.set_xlim(0, None)

                ax.legend(loc=legend_location, fontsize='xx-large')

                fig.tight_layout()
                plt.savefig('dos.png', dpi=300, format='png', bbox_inches='tight')
                plt.show()

            elif option == 'line':
                if bulk_dos is not None:
                    if bulk_legend is not None:
                        label = bulk_legend
                    else:
                        label = 'bulk dos'

                    if bulk_option == 'line':
                        ax.plot(bulk_dos[:, 1] * proportion, bulk_dos[:, 0],
                                color='black', linestyle=':', linewidth=4, label=label)
                    else:
                        ax.fill(bulk_dos[:, 1] * proportion, bulk_dos[:, 0], 'black', alpha=0.3, label=label)

                if plot_total:
                    ax.plot(self._tdos, self._freq,
                            color=color,
                            label='total dos',
                            linewidth=4)

                if label_position == 'right':
                    ax.yaxis.tick_right()
                    ax.yaxis.set_ticks_position('both')
                    ax.yaxis.set_label_position("right")

                __pdos = np.zeros((len(atoms), self._freq.shape[0]))
                ind_modes = []
                try:
                    for ind_set, atom_set in enumerate(atoms):
                        ind_modes.append([ind for ind, val in enumerate(self.process.unit_cell.xyz_true)
                                          if val in atom_set])
                        for ind_mode in ind_modes[ind_set]:
                            __pdos[ind_set, :] += self._pdos[ind_mode, :]
                except TypeError:
                    print("\nFail to identify the index of atoms='{0}'".format(atoms))
                    print("Please designate the 'atoms' parameter for which DOS will be projected\n")
                    raise

                for ind_set, _pdos_ in enumerate(__pdos):
                    if legends is None:
                        label = 'atom-{0}'.format(set(atoms[ind_set]))
                    else:
                        label = legends[ind_set]

                    ax.plot(_pdos_, self._freq,
                            linestyle='--',
                            linewidth=4,
                            label=label)

                ax.set_ylabel('Frequency (THz)', fontdict=font_x)
                ax.set_ylim(y_min, y_max)
                labels = [item for item in ax.get_yticks()]
                ytick_labels = [round(float(label), 2) for label in labels]
                ax.set_yticks(ytick_labels)
                ax.set_yticklabels(labels=ytick_labels, fontdict=font_tick)
                ax.set_ylim(y_min, y_max)

                ax.set_xlabel('DOS (a.u.)', fontdict=font_y)
                ax.set_xticks([])
                ax.set_xlim(0, None)

                ax.legend(loc=legend_location, fontsize='xx-large')

                fig.tight_layout()
                plt.savefig('dos.png', dpi=300, format='png', bbox_inches='tight')
                plt.show()

            elif option == 'stack':
                if bulk_dos is not None:
                    if bulk_legend is not None:
                        label = bulk_legend
                    else:
                        label = 'bulk dos'

                    if bulk_option == 'line':
                        ax.plot(bulk_dos[:, 1] * proportion, bulk_dos[:, 0],
                                color='black', linestyle=':', linewidth=4, label=label)
                    else:
                        ax.fill(bulk_dos[:, 1] * proportion, bulk_dos[:, 0], 'black', alpha=0.3, label=label)

                if plot_total:
                    ax.plot(self._tdos, self._freq,
                            color=color,
                            label='total dos',
                            linewidth=4)

                if label_position == 'right':
                    ax.yaxis.tick_right()
                    ax.yaxis.set_ticks_position('both')
                    ax.yaxis.set_label_position("right")

                __pdos = np.zeros((len(atoms), self._freq.shape[0]))
                ind_modes = []
                try:
                    for ind_set, atom_set in enumerate(atoms):
                        ind_modes.append([ind for ind, val in enumerate(self.process.unit_cell.xyz_true)
                                          if val in atom_set])
                        for ind_mode in ind_modes[ind_set]:
                            __pdos[ind_set, :] += self._pdos[ind_mode, :]
                except TypeError:
                    print("\nFail to identify the index of atoms='{0}'".format(atoms))
                    print("Please designate the 'atoms' parameter for which DOS will be projected\n")
                    raise

                labels = []
                x_data = np.cumsum(__pdos, axis=0)
                for ind, atom in enumerate(atoms):
                    if legends is None:
                        labels.append('atom-{0}'.format(set(atom)))
                    else:
                        labels.append(legends[ind])

                    ax.fill_betweenx(self._freq, x_data[ind, :], linewidth=1, label=labels[ind], zorder=-ind)

                ax.set_ylabel('Frequency (THz)', fontdict=font_x)
                ax.set_ylim(y_min, y_max)
                labels = [item for item in ax.get_yticks()]
                ytick_labels = [round(float(label), 2) for label in labels]
                ax.set_yticks(ytick_labels)
                ax.set_yticklabels(labels=ytick_labels, fontdict=font_tick)
                ax.set_ylim(y_min, y_max)

                ax.set_xlabel('DOS (a.u.)', fontdict=font_y)
                ax.set_xticks([])
                ax.set_xlim(0, None)

                ax.legend(loc=legend_location, fontsize='xx-large')

                fig.tight_layout()
                plt.savefig('dos.png', dpi=300, format='png', bbox_inches='tight')
                plt.show()
