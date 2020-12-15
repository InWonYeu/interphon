import numpy as np


class Band(object):
    def __init__(self, process):
        self.process = process
        self.k_points_length = np.zeros(len(self.process.k_points))
        self.projected_w = np.empty((len(self.process.k_points),
                                     len(self.process.unit_cell.xyz_true),
                                     len(self.process.unit_cell.xyz_true)), dtype=float)
        self.ind_high_sym = [0]

    def set(self):
        for ind, kpt in enumerate(self.process.k_points[1:], 1):
            self.k_points_length[ind] = self.k_points_length[ind - 1] \
                                        + np.sqrt(np.dot(self.process.k_points[ind] - self.process.k_points[ind - 1],
                                                         self.process.k_points[ind] - self.process.k_points[ind - 1]))

            if np.allclose(self.process.k_points[ind], self.process.k_points[ind - 1]):
                self.ind_high_sym.append(ind)
        self.ind_high_sym.append(ind)

        for ind_kpt, eig_freqs in enumerate(self.process.w_q):
            for ind_freq, _ in enumerate(eig_freqs):
                for ind_mode, _ in enumerate(self.process.unit_cell.xyz_true):
                    self.projected_w[ind_kpt, ind_freq, ind_mode] = \
                        1 * (abs(self.process.v_q[ind_kpt, ind_freq, ind_mode]) ** 2)

    def write(self, out_folder='.'):
        with open(out_folder + '/band.dat', 'w') as outfile:
            comment = "Phonon Band"
            outfile.write("%s" % comment + '\n')
            outfile.write("%s" %
                          '    K_Points_Path' +
                          '    Frequency (THz)' + '\n')

            for x, y_set in zip(self.k_points_length, self.process.w_q):
                line = ' %16.9f ' % x
                for y in y_set:
                    line = line + ' %16.9f ' % y
                outfile.write("%s" % line + '\n')

    def plot(self,
             k_labels=[],
             atoms=None,
             elimit=None,
             color='tab:orange',
             option='plain',
             cmap='jet',
             colorbar_label=None,
             colorbar_location='right'):

        from matplotlib import pyplot as plt

        # plt.rc("font", size=22)  # controls default text sizes
        # plt.rc('axes', linewidth=2)

        font_y = {'family': 'Arial', 'size': 36, 'color': 'black', 'weight': 'bold'}
        font_tick = {'size': 36, 'color': 'black'}
        font_bar = {'size': 24, 'color': 'black'}
        font_bar_title = {'family': 'Arial', 'size': 30, 'color': 'black', 'weight': 'bold'}

        if elimit is None:
            y_min = np.floor(self.process.w_q.min()) - 1
            y_max = np.ceil(self.process.w_q.max()) + 1
        else:
            y_min = elimit[0]
            y_max = elimit[1]

        _k_labels = []
        _ind_pbc = self.process.user_arg.periodicity.nonzero()[0]

        if _ind_pbc.shape[0] == 2:
            for k_label in k_labels:
                if k_label == 'G' or k_label == 'g':
                    _k_labels.append(r'$\overline{\Gamma}$')
                else:
                    _k_labels.append(r'$\overline{' + k_label + '}$')

        else:
            for k_label in k_labels:
                if k_label == 'G' or k_label == 'g':
                    _k_labels.append(r'$\Gamma$')
                else:
                    _k_labels.append(r'$' + k_label + '$')

        self.process.num_figure = self.process.num_figure + 1
        fig = plt.figure(self.process.num_figure, figsize=(13, 9))

        if option == 'plain':
            ax = fig.subplots()
            ax.plot(self.k_points_length, self.process.w_q, color=color, linewidth=2)

            # ax.set_xlabel('K-points', fontdict=font_x)
            ax.set_xlim(self.k_points_length[0], self.k_points_length[-1])
            ax.set_xticks(self.k_points_length[self.ind_high_sym])
            ax.set_xticklabels(_k_labels, fontdict=font_tick)

            ax.set_ylabel('Frequency (THz)', fontdict=font_y)
            ax.set_ylim(y_min, y_max)
            labels = [item for item in ax.get_yticks()]
            ytick_labels = [round(float(label), 2) for label in labels]
            ax.set_yticks(ytick_labels)
            ax.set_yticklabels(labels=ytick_labels, fontdict=font_tick)
            ax.set_ylim(y_min, y_max)

            ax.grid(True)

            fig.tight_layout()
            plt.savefig('band.png', dpi=300, format='png', bbox_inches='tight')
            plt.show()

        elif option == 'projection':
            from matplotlib.collections import LineCollection

            _projected_w = np.zeros((len(self.process.k_points),
                                     len(self.process.unit_cell.xyz_true),
                                     1))
            try:
                ind_modes = [ind for ind, val in enumerate(self.process.unit_cell.xyz_true) if val in atoms]
                for ind_mode in ind_modes:
                    _projected_w[:, :, 0] += self.projected_w[:, :, ind_mode]
            except TypeError:
                print("\nFail to identify the index of atoms='{0}'".format(atoms))
                print("Please designate the 'atoms' parameter for which band will be projected\n")
                raise

            if colorbar_location == 'right':
                ax = fig.subplots()
                for ind_freq in range(self.process.w_q.shape[1]):
                    _x = self.k_points_length
                    _y = self.process.w_q[:, ind_freq]

                    points = np.array([_x, _y]).T.reshape(-1, 1, 2)
                    segments = np.concatenate([points[:-1], points[1:]], axis=1)

                    norm = plt.Normalize(vmin=0, vmax=1)
                    lc = LineCollection(segments, cmap=cmap, norm=norm)
                    lc.set_array(_projected_w[:, ind_freq, 0])
                    lc.set_linewidth(2)
                    line = ax.add_collection(lc)

                cbar = fig.colorbar(line, ax=ax, extend='max', ticks=[0, 0.2, 0.4, 0.6, 0.8, 1])
                cbar.ax.set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1], fontdict=font_bar)

                if colorbar_label is None:
                    label = 'Contribution of atom-{0}'.format(set(atoms))
                else:
                    label = 'Contribution of {0}'.format(colorbar_label)

                cbar.set_label(label, fontdict=font_bar_title)

            elif colorbar_location == 'bottom':
                ax, cax = fig.subplots(nrows=2, gridspec_kw={"height_ratios": [1, 0.05]})
                for ind_freq in range(self.process.w_q.shape[1]):
                    _x = self.k_points_length
                    _y = self.process.w_q[:, ind_freq]

                    points = np.array([_x, _y]).T.reshape(-1, 1, 2)
                    segments = np.concatenate([points[:-1], points[1:]], axis=1)

                    norm = plt.Normalize(vmin=0, vmax=1)
                    lc = LineCollection(segments, cmap=cmap, norm=norm)
                    lc.set_array(_projected_w[:, ind_freq, 0])
                    lc.set_linewidth(2)
                    line = ax.add_collection(lc)

                cbar = fig.colorbar(line, cax=cax, extend='max', ticks=[0, 0.2, 0.4, 0.6, 0.8, 1],
                                    orientation="horizontal")
                cbar.cax.set_xticklabels([0, 0.2, 0.4, 0.6, 0.8, 1], fontdict=font_bar)

                if colorbar_label is None:
                    label = 'Contribution of atom-{0}'.format(set(atoms))
                else:
                    label = 'Contribution of {0}'.format(colorbar_label)

                cbar.set_label(label, fontdict=font_bar_title)

            # ax.set_xlabel('K-points', fontdict=font_x)
            ax.set_xlim(self.k_points_length[0], self.k_points_length[-1])
            ax.set_xticks(self.k_points_length[self.ind_high_sym])
            ax.set_xticklabels(_k_labels, fontdict=font_tick)

            ax.set_ylabel('Frequency (THz)', fontdict=font_y)
            ax.set_ylim(y_min, y_max)
            labels = [item for item in ax.get_yticks()]
            ytick_labels = [round(float(label), 2) for label in labels]
            ax.set_yticks(ytick_labels)
            ax.set_yticklabels(labels=ytick_labels, fontdict=font_tick)
            ax.set_ylim(y_min, y_max)

            ax.grid(True)

            fig.tight_layout()
            plt.savefig('band.png', dpi=300, format='png', bbox_inches='tight')
            plt.show()

    def plot_with_dos(self,
                      k_labels=[],
                      band_atoms=None,
                      elimit=None,
                      band_color='tab:orange',
                      band_option='plain',
                      cmap='jet',
                      colorbar_label=None,
                      dos_object=None,
                      dos_atoms=None,
                      dos_color='tab:orange',
                      dos_option='plain',
                      dos_legends=None,
                      legend_location='upper right',
                      bulk_dos=None,
                      bulk_option='fill',
                      bulk_legend=None,
                      proportion=1):

        from matplotlib import pyplot as plt

        font_x = {'family': 'Arial', 'size': 36, 'color': 'black', 'weight': 'bold'}
        font_y = {'family': 'Arial', 'size': 36, 'color': 'black', 'weight': 'bold'}
        font_tick = {'size': 30, 'color': 'black'}

        self.process.num_figure = self.process.num_figure + 1
        fig = plt.figure(self.process.num_figure, (14, 9))

        if elimit is None:
            y_min = np.floor(self.process.w_q.min()) - 1
            y_max = np.ceil(self.process.w_q.max()) + 1
        else:
            y_min = elimit[0]
            y_max = elimit[1]

        _k_labels = []
        _ind_pbc = self.process.user_arg.periodicity.nonzero()[0]

        if _ind_pbc.shape[0] == 2:
            for k_label in k_labels:
                if k_label == 'G' or k_label == 'g':
                    _k_labels.append(r'$\overline{\Gamma}$')
                else:
                    _k_labels.append(r'$\overline{' + k_label + '}$')

        else:
            for k_label in k_labels:
                if k_label == 'G' or k_label == 'g':
                    _k_labels.append(r'$\Gamma$')
                else:
                    _k_labels.append(r'$' + k_label + '$')

        # Band
        # ax_dos = plt.subplot(121)
        if band_option == 'plain':
            (ax_band, ax_dos) = fig.subplots(nrows=1, ncols=2,
                                             gridspec_kw={"width_ratios": [7, 3]},
                                             sharey=True)

            ax_band.plot(self.k_points_length, self.process.w_q, color=band_color, linewidth=2)

            ax_band.yaxis.set_ticks_position('both')

            # ax_band.set_xlabel('K-points', fontdict=font_x)
            ax_band.set_xlim(self.k_points_length[0], self.k_points_length[-1])
            ax_band.set_xticks(self.k_points_length[self.ind_high_sym])
            ax_band.set_xticklabels(_k_labels, fontdict=font_tick)

            ax_band.set_ylabel('Frequency (THz)', fontdict=font_y)
            ax_band.set_ylim(y_min, y_max)
            labels = [item for item in ax_band.get_yticks()]
            ytick_labels = [round(float(label), 2) for label in labels]
            ax_band.set_yticks(ytick_labels)
            ax_band.set_yticklabels(labels=ytick_labels, fontdict=font_tick)
            ax_band.set_ylim(y_min, y_max)

            ax_band.grid(True)

        elif band_option == 'projection':
            from matplotlib.collections import LineCollection
            from mpl_toolkits.axes_grid1.inset_locator import inset_axes

            font_y = {'family': 'Arial', 'size': 36, 'color': 'black', 'weight': 'bold'}
            font_tick = {'size': 36, 'color': 'black'}
            font_bar = {'size': 24, 'color': 'black'}
            font_bar_title = {'family': 'Arial', 'size': 30, 'color': 'black', 'weight': 'bold'}

            (ax_band, ax_dos) = fig.subplots(nrows=1, ncols=2,
                                             gridspec_kw={"width_ratios": [7, 3]},
                                             sharey=True)

            _projected_w = np.zeros((len(self.process.k_points),
                                     len(self.process.unit_cell.xyz_true),
                                     1))
            try:
                ind_modes = [ind for ind, val in enumerate(self.process.unit_cell.xyz_true) if val in band_atoms]
                for ind_mode in ind_modes:
                    _projected_w[:, :, 0] += self.projected_w[:, :, ind_mode]
            except TypeError:
                print("\nFail to identify the index of atoms='{0}'".format(band_atoms))
                print("Please designate the 'atoms' parameter for which band will be projected\n")
                raise

            for ind_freq in range(self.process.w_q.shape[1]):
                _x = self.k_points_length
                _y = self.process.w_q[:, ind_freq]

                points = np.array([_x, _y]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)

                norm = plt.Normalize(vmin=0, vmax=1)

                lc = LineCollection(segments, cmap=cmap, norm=norm)
                lc.set_array(_projected_w[:, ind_freq, 0])
                lc.set_linewidth(2)

                line = ax_band.add_collection(lc)

            # fig.subplots_adjust(bottom=0.1)
            # cbar_ax = fig.add_axes([0.1, 0.02, 0.6, 0.06])
            fig.tight_layout()
            cbaxes = inset_axes(ax_band, width="70%", height="4%", loc='lower center',
                                bbox_to_anchor=(0.0, 0.04, 1, 1), bbox_transform=ax_band.transAxes)  # upper: 0.88, lower: 0.04
            cbar = fig.colorbar(line, cax=cbaxes, orientation="horizontal",
                                ticks=[0, 0.2, 0.4, 0.6, 0.8, 1], extend='max')
            cbar.ax.tick_params(direction='out')
            cbar.ax.set_xticklabels([0, 0.2, 0.4, 0.6, 0.8, 1], fontdict=font_bar)

            if colorbar_label is None:
                label = 'Contribution of atom-{0}'.format(set(band_atoms))
            else:
                label = 'Contribution of {0}'.format(colorbar_label)

            cbaxes.set_title(label, fontdict=font_bar_title)

            ax_band.yaxis.set_ticks_position('both')

            # ax_band.set_xlabel('K-points', fontdict=font_x)
            ax_band.set_xlim(self.k_points_length[0], self.k_points_length[-1])
            ax_band.set_xticks(self.k_points_length[self.ind_high_sym])
            ax_band.set_xticklabels(_k_labels, fontdict=font_tick)

            ax_band.set_ylabel('Frequency (THz)', fontdict=font_y)
            ax_band.set_ylim(y_min, y_max)
            labels = [item for item in ax_band.get_yticks()]
            ytick_labels = [round(float(label), 2) for label in labels]
            ax_band.set_yticks(ytick_labels)
            ax_band.set_yticklabels(labels=ytick_labels, fontdict=font_tick)
            ax_band.set_ylim(y_min, y_max)

            ax_band.grid(True)

        # DOS
        # ax_dos = plt.subplot(122, sharey=ax_band)
        fig.subplots_adjust(wspace=0.1)

        if dos_option == 'plain':
            if bulk_dos is not None:
                if bulk_legend is not None:
                    label = bulk_legend
                else:
                    label = 'bulk_dos'

                if bulk_option == 'line':
                    ax_dos.plot(bulk_dos[:, 1] * proportion, bulk_dos[:, 0],
                                color='black', linestyle=':', linewidth=4, label=label)
                else:
                    ax_dos.fill(bulk_dos[:, 1] * proportion, bulk_dos[:, 0], 'black', alpha=0.3, label=label)

            ax_dos.plot(dos_object.tdos, dos_object.freq, color=dos_color, label='total dos', linewidth=4)
            ax_dos.yaxis.set_ticks_position('both')

            ax_dos.set_xlabel('DOS (a.u.)', fontdict=font_y)
            ax_dos.set_xticks([])
            ax_dos.set_xlim(0, None)

            ax_dos.legend(loc=legend_location, fontsize='xx-large')

        elif dos_option == 'line':
            if bulk_dos is not None:
                if bulk_legend is not None:
                    label = bulk_legend
                else:
                    label = 'bulk_dos'

                if bulk_option == 'line':
                    ax_dos.plot(bulk_dos[:, 1] * proportion, bulk_dos[:, 0],
                                color='black', linestyle=':', linewidth=4, label=label)
                else:
                    ax_dos.fill(bulk_dos[:, 1] * proportion, bulk_dos[:, 0], 'black', alpha=0.3, label=label)

            ax_dos.plot(dos_object.tdos, dos_object.freq, color=dos_color, label='total dos', linewidth=4)
            ax_dos.yaxis.set_ticks_position('both')

            __pdos = np.zeros((len(dos_atoms), dos_object.freq.shape[0]))
            ind_modes = []
            try:
                for ind_set, atom_set in enumerate(dos_atoms):
                    ind_modes.append([ind for ind, val in enumerate(self.process.unit_cell.xyz_true)
                                      if val in atom_set])
                    for ind_mode in ind_modes[ind_set]:
                        __pdos[ind_set, :] += dos_object.pdos[ind_mode, :]
            except TypeError:
                print("\nFail to identify the index of atoms='{0}'".format(dos_atoms))
                print("Please designate the 'atoms' parameter for which DOS will be projected\n")
                raise

            for ind_set, _pdos_ in enumerate(__pdos):
                if dos_legends is None:
                    label = 'atom-{0}'.format(set(dos_atoms[ind_set]))
                else:
                    label = dos_legends[ind_set]

                ax_dos.plot(_pdos_, dos_object.freq, linestyle='--', linewidth=4,
                            label=label)

            ax_dos.set_xlabel('DOS (a.u.)', fontdict=font_y)
            ax_dos.set_xticks([])
            ax_dos.set_xlim(0, None)

            ax_dos.legend(loc=legend_location, fontsize='xx-large')

        elif dos_option == 'stack':
            if bulk_dos is not None:
                if bulk_legend is not None:
                    label = bulk_legend
                else:
                    label = 'bulk_dos'

                if bulk_option == 'line':
                    ax_dos.plot(bulk_dos[:, 1] * proportion, bulk_dos[:, 0],
                                color='black', linestyle=':', linewidth=4, label=label)
                else:
                    ax_dos.fill(bulk_dos[:, 1] * proportion, bulk_dos[:, 0], 'black', alpha=0.3, label=label)

            ax_dos.plot(dos_object.tdos, dos_object.freq, color=dos_color, label='total dos', linewidth=4)
            ax_dos.yaxis.set_ticks_position('both')

            __pdos = np.zeros((len(dos_atoms), dos_object.freq.shape[0]))
            ind_modes = []
            try:
                for ind_set, atom_set in enumerate(dos_atoms):
                    ind_modes.append([ind for ind, val in enumerate(self.process.unit_cell.xyz_true)
                                      if val in atom_set])
                    for ind_mode in ind_modes[ind_set]:
                        __pdos[ind_set, :] += dos_object.pdos[ind_mode, :]
            except TypeError:
                print("\nFail to identify the index of atoms='{0}'".format(dos_atoms))
                print("Please designate the 'atoms' parameter for which DOS will be projected\n")
                raise

            labels = []
            x_data = np.cumsum(__pdos, axis=0)
            for ind, atom in enumerate(dos_atoms):
                if dos_legends is None:
                    labels.append('atom-{0}'.format(set(atom)))
                else:
                    labels.append(dos_legends[ind])

                ax_dos.fill_betweenx(dos_object.freq, x_data[ind, :], linewidth=1, label=labels[ind], zorder=-ind)

            ax_dos.set_xlabel('DOS (a.u.)', fontdict=font_y)
            ax_dos.set_xticks([])
            ax_dos.set_xlim(0, None)

            ax_dos.legend(loc=legend_location, fontsize='xx-large')

        # fig.tight_layout()
        plt.savefig('band_dos.png', dpi=300, format='png', bbox_inches='tight')
        plt.show()
