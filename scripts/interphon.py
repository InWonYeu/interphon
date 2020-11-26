import os
import click
import glob
import yaml
try:
    from InterPhon.core import PreProcess, PostProcess
except ImportError:  # parent class of ModuleNotFoundError
    print("\nInterPhon package should be in one of the following directories: \n 1) current folder \n "
          "2) standard library modules \n 3) third party modules \n")
    raise


def iter_lattice_yaml(lattice_matrix):
    for ind, val in enumerate(lattice_matrix, 1):
        lattice = ('a{0}'.format(ind), '{0:>20.16f}, {1:>20.16f}, {2:20.16f}'.format(val[0], val[1], val[2]).strip())
        yield lattice


def iter_atom_yaml(unit_cell):
    for ind, atom_type in enumerate(unit_cell.atom_type):
        val = unit_cell.atom_cart[ind]
        atom = {'index': ind,
                'type': atom_type,
                'position': '{0:>20.16f}, {1:>20.16f}, {2:20.16f}'.format(val[0], val[1], val[2]).strip(),
                'selection': True if ind in unit_cell.atom_true else False}
        yield atom


@click.command()
@click.argument('argument_file',
                type=click.Path(exists=True),
                required=False)
@click.option('--preprocess/--postprocess', '-pre/-post', 'process',
              default=True,
              required=True,
              show_default=True)
@click.option('--dft_code', '-dft', 'dft',
              default='vasp',
              type=click.Choice(['vasp', 'espresso', 'aims']),
              help='DFT code name.',
              required=True,
              show_default=True)
@click.option('--displacement', '-disp',
              default='0.01',
              type=click.STRING,
              help='Displacement length (unit: Angst).',
              required=True,
              show_default=True)
@click.option('--enlargement', '-enlarge',
              default='1 1 1',
              type=click.STRING,
              help='Extension ratio along each a, b, c lattice.',
              required=True,
              show_default=True)
@click.option('--periodicity', '-pbc',
              default='1 1 0',
              type=click.STRING,
              help='Periodic (True) or not (False) along each a, b, c lattice direction.',
              required=True,
              show_default=True)
@click.option('--unitcell', '-c',
              # prompt='unit cell file',
              type=click.Path(exists=True),
              help='Unit cell file path.')
@click.option('--supercell', '-sc',
              type=click.Path(exists=True),
              help='Super cell file.')
@click.option('--forces', '-fc',
              type=click.STRING,
              help='Force files.')
# Options for DOS write and plot
@click.option('--density_of_state', '-dos', 'dos', is_flag=True,
              default=False,
              required=True,
              show_default=True)
@click.option('--kpoint_dos', '-kdos',
              type=click.Path(exists=True),
              help='K-point file for DOS.')
@click.option('--sigma', '-sig',
              default=0.1,
              type=click.FLOAT,
              help='Sigma of gaussian smearing (0.0: tetrahedron method).',
              required=True,
              show_default=True)
@click.option('--number_dos', '-ndos', 'num_dos',
              default=200,
              type=click.INT,
              help='The number of DOS points.',
              required=True,
              show_default=True)
@click.option('--projection_atom_dos', '-atom_dos', 'atom_dos',
              type=click.STRING,
              help='The Index of atoms to be projected.')
@click.option('--projection_legend_dos', '-legend_dos', 'legend_dos',
              type=click.STRING,
              help='Legends for the projected atoms.')
@click.option('--energy_limit', '-elimit', 'elimit',
              type=click.STRING,
              help='Energy (eV) limitation of DOS and Band plot.')
@click.option('--tdos_color_dos', '-color_dos', 'color_dos',
              default='tab:orange',
              type=click.STRING,
              help='The color of total DOS line.',
              show_default=True)
@click.option('--projection_option_dos', '-option_dos', 'option_dos',
              default='plain',
              type=click.Choice(['plain', 'line', 'stack']),
              help='Option for DOS projection plot.',
              show_default=True)
@click.option('--image_orientation_dos', '-orientation_dos', 'orientation_dos',
              default='horizontal',
              type=click.Choice(['horizontal', 'vertical']),
              help='Orientation of DOS plot.',
              show_default=True)
@click.option('--legend_location_dos', '-legend_loc_dos', 'legend_loc_dos',
              default='best',
              type=click.Choice(['best', 'upper right', 'upper left', 'lower left', 'lower right', 'right',
                                 'center left', 'center right', 'lower center', 'upper center', 'center']),
              help='Location of DOS legend.',
              show_default=True)
# Options for Thermal_Property write and plot
@click.option('--thermal_property', '-thermal', 'thermal', is_flag=True,
              default=False,
              required=True,
              show_default=True)
@click.option('--temperature_minimum', '-tmin', 'tmin',
              default=0,
              type=click.INT,
              help='Temperature minimum (unit: K).',
              required=True,
              show_default=True)
@click.option('--temperature_maximum', '-tmax', 'tmax',
              default=1000,
              type=click.INT,
              help='Temperature maximum (unit: K).',
              required=True,
              show_default=True)
@click.option('--temperature_step', '-tstep', 'tstep',
              default=10,
              type=click.INT,
              help='Temperature step (unit: K).',
              required=True,
              show_default=True)
# Options for Band write and plot
@click.option('--phonon_band', '-band', 'band', is_flag=True,
              default=False,
              required=True,
              show_default=True)
@click.option('--kpoint_band', '-kband',
              type=click.Path(exists=True),
              help='K-point file for Band.')
@click.option('--kpoint_label_band', '-k_label_band', 'k_label_band',
              type=click.STRING,
              help='The label of high-symmetry k-points.')
@click.option('--projection_atom_band', '-atom_band', 'atom_band',
              type=click.STRING,
              help='The Index of atoms to be projected.')
@click.option('--total_color_band', '-color_band', 'color_band',
              default='tab:orange',
              type=click.STRING,
              help='The color of Band line.',
              show_default=True)
@click.option('--projection_option_band', '-option_band', 'option_band',
              default='plain',
              type=click.Choice(['plain', 'projection']),
              help='Option for Band projection plot.',
              show_default=True)
@click.option('--colorbar_label_band', '-bar_label_band', 'bar_label_band',
              type=click.STRING,
              help='The label of colorbar for projection plot.')
@click.option('--colorbar_location_band', '-bar_loc_band', 'bar_loc_band',
              default='right',
              type=click.Choice(['right', 'bottom']),
              help='Location of colorbar.',
              show_default=True)
# Options for Mode write and plot
@click.option('--phonon_mode', '-mode', 'mode', is_flag=True,
              default=False,
              required=True,
              show_default=True)
@click.option('--index_mode', '-ind_mode', 'ind_mode',
              default='0',
              type=click.STRING,
              help='The index of phonon mode(kpoint, index).',
              show_default=True)
@click.option('--k_point_mode', '-kpt_mode', 'kpt_mode',
              default='0.0 0.0 0.0',
              type=click.STRING,
              help='The K-point of phonon mode(kpoint, index).',
              show_default=True)
def main(argument_file, process, dft, displacement, enlargement, periodicity,
         unitcell, supercell, forces, kpoint_dos,
         dos, sigma, num_dos, atom_dos, legend_dos, elimit, color_dos, option_dos, orientation_dos, legend_loc_dos,
         thermal, tmin, tmax, tstep,
         band, kpoint_band, k_label_band, atom_band, color_band, option_band, bar_label_band, bar_loc_band,
         mode, ind_mode, kpt_mode):
    if argument_file is not None:
        with open(argument_file, 'r') as infile:
            lines = infile.readlines()

        arg_dict = {}
        for line in lines:
            key = line.strip().split('=')[0]
            value = line.strip().split('=')[1]
            arg_dict[key.strip().lower()] = value.strip()

        for key, value in arg_dict.items():
            if key in ('dft', 'dft_code'):
                if value in ('vasp', 'espresso', 'aims'):
                    dft = value
                else:
                    raise Exception('invalid choice: {0}. (choose from vasp, espresso, aims)'.format(value))
            elif key in ('displacement', 'disp'):
                displacement = value
            elif key in ('enlargement', 'enlarge'):
                enlargement = value
            elif key in ('periodicity', 'pbc'):
                periodicity = value
            elif key in ('unitcell', 'c'):
                unitcell = value
            elif key in ('supercell', 'sc'):
                supercell = value
            elif key in ('forces', 'fc'):
                forces = value
            elif key in ('kpoint_dos', 'kdos'):
                kpoint_dos = value
            elif key in ('sigma', 'sig'):
                sigma = float(value)
            elif key in ('number_dos', 'ndos'):
                num_dos = int(value)
            elif key in ('projection_atom_dos', 'atom_dos'):
                atom_dos = value
            elif key in ('projection_legend_dos', 'legend_dos'):
                legend_dos = value
            elif key in ('energy_limit', 'elimit'):
                elimit = value
            elif key in ('tdos_color_dos', 'color_dos'):
                color_dos = value
            elif key in ('projection_option_dos', 'option_dos'):
                if value in ('plain', 'line', 'stack'):
                    option_dos = value
                else:
                    raise Exception('invalid choice: {0}. (choose from plain, line, stack)'.format(value))
            elif key in ('image_orientation_dos', 'orientation_dos'):
                if value in ('horizontal', 'vertical'):
                    orientation_dos = value
                else:
                    raise Exception('invalid choice: {0}. (choose from horizontal, vertical)'.format(value))
            elif key in ('legend_location_dos', 'legend_loc_dos'):
                if value in ('best', 'upper right', 'upper left', 'lower left', 'lower right', 'right',
                             'center left', 'center right', 'lower center', 'upper center', 'center'):
                    legend_loc_dos = value
                else:
                    raise Exception('invalid choice: {0}. (choose from best, upper right, upper left, '
                                    'lower left, lower right, right, center left, center right, lower center, '
                                    'upper center, center)'.format(value))
            elif key in ('temperature_minimum', 'tmin'):
                tmin = int(value)
            elif key in ('temperature_maximum', 'tmax'):
                tmax = int(value)
            elif key in ('temperature_step', 'tstep'):
                tstep = int(value)
            elif key in ('kpoint_band', 'kband'):
                kpoint_band = value
            elif key in ('kpoint_label_band', 'k_label_band'):
                k_label_band = value
            elif key in ('projection_atom_band', 'atom_band'):
                atom_band = value
            elif key in ('total_color_band', 'color_band'):
                color_band = value
            elif key in ('projection_option_band', 'option_band'):
                if value in ('plain', 'projection'):
                    option_band = value
                else:
                    raise Exception('invalid choice: {0}. (choose from plain, projection)'.format(value))
            elif key in ('colorbar_label_band', 'bar_label_band'):
                bar_label_band = value
            elif key in ('colorbar_location_band', 'bar_loc_band'):
                if value in ('right', 'bottom'):
                    bar_loc_band = value
                else:
                    raise Exception('invalid choice: {0}. (choose from right, bottom)'.format(value))
            elif key in ('index_mode', 'ind_mode'):
                ind_mode = value
            elif key in ('k_point_mode', 'kpt_mode'):
                kpt_mode = value

    # check the number of arguments
    if process:
        _process = 'pre-process'
        print('\n#########################################')
        print('\tChecking pre-arguments...')
        print('#########################################')

        user_args = {'dft_code': dft,
                     'displacement': displacement,
                     'enlargement': enlargement,
                     'periodicity': periodicity}

        files = {'unit_cell_file': os.path.abspath(unitcell)}
        working_dir = os.path.dirname(files.get('unit_cell_file'))

        print('\n[User Arguments for Pre-Process]')
        for key, value in user_args.items():
            print('{0} = {1}'.format(key.title(), value))

        print('\n[Files for Pre-Process]')
        for key, value in files.items():
            if isinstance(value, str):
                value = os.path.basename(value)
            print('{0} = {1}'.format(key.title(), value))

    else:
        _process = 'post-process'
        print('\n#########################################')
        print('\tChecking post-arguments...')
        print('#########################################')

        user_args = {'dft_code': dft,
                     'displacement': displacement,
                     'enlargement': enlargement,
                     'periodicity': periodicity}

        files = {'unit_cell_file': os.path.abspath(unitcell)}

        dos_args = {'flag': False}
        thermal_args = {'flag': False}
        band_args = {'flag': False}
        mode_args = {'flag': False}

        # Check and fill the files dictionary
        # Supercell file
        if supercell is None:
            raise Exception('Supercell file should be given')
        else:
            files['super_cell_file'] = os.path.abspath(supercell)

        # K-point file
        if kpoint_dos is None:
            if kpoint_band is None:
                raise Exception('K-points file should be given')
            else:
                files['k_point_file_band'] = os.path.abspath(kpoint_band)  # only Band
        else:
            if kpoint_band is None:
                files['k_point_file_dos'] = os.path.abspath(kpoint_dos)  # only DOS
            else:
                files['k_point_file_dos'] = os.path.abspath(kpoint_dos)
                files['k_point_file_band'] = os.path.abspath(kpoint_band)  # Both DOS and Band

        # Force file
        if forces is None:
            raise Exception('Force files should be given')
        else:
            force_path = glob.glob(forces)
            force_path.sort()
            if not force_path:
                raise Exception('Path "{0}" does not exist.'.format(forces))
            else:
                files['force_file'] = force_path

        working_dir = os.path.dirname(files.get('unit_cell_file'))

        print('\n[User Arguments for Post-Process]')
        for key, value in user_args.items():
            print('{0} = {1}'.format(key.title(), value))

        print('\n[Files for Post-Process]')
        for key, value in files.items():
            if isinstance(value, str):
                value = os.path.basename(value)
            print('{0} = {1}'.format(key.title(), value))

        print('\n[User Arguments for Phonon DOS]')
        if dos:
            # Check and fill the dos_args dictionary
            dos_args['flag'] = dos
            dos_args['sigma'] = float(sigma)
            dos_args['num_dos'] = int(num_dos)
            dos_args['color'] = color_dos
            dos_args['option'] = option_dos
            dos_args['orientation'] = orientation_dos
            dos_args['legend_loc'] = legend_loc_dos
            if atom_dos is None:
                if option_dos != 'plain':
                    print('Caution:')
                    print('"atom_dos" should be given, in the Index of selected atoms, to plot projected DOS.')
                    print('"option_dos" is changed from "{0}" to "plain".'.format(option_dos))
                    dos_args['option'] = 'plain'
                dos_args['atom'] = None
            else:
                _atom_dos = atom_dos.strip().split(',')
                atom_dos = []
                for atom in _atom_dos:
                    atom_dos.append([int(value) for value in atom.strip().split()])
                dos_args['atom'] = atom_dos

            if legend_dos is None:
                # if dos_args.get('option') != 'plain':
                #     legend_dos = []
                #     for atom_set in dos_args.get('atom'):
                #         legend_dos.append('group-' + str(atom_set))
                #     dos_args['legend'] = legend_dos
                # else:
                dos_args['legend'] = None
            else:
                _legend_dos = legend_dos.strip().split(',')
                legend_dos = [value.strip() for value in _legend_dos]
                if dos_args.get('option') == 'plain':
                    dos_args['legend'] = None
                else:
                    if len(legend_dos) == len(dos_args.get('atom')):
                        dos_args['legend'] = legend_dos
                    else:
                        raise Exception('The length of "{0}" should be matched with that of "{1}".'.format(legend_dos,
                                                                                                           atom_dos))

            if elimit is None:
                dos_args['elimit'] = None
            else:
                _elimit = elimit.strip().split()
                elimit = [int(value) for value in _elimit]
                if len(elimit) == 2:
                    dos_args['elimit'] = elimit
                else:
                    raise Exception('"elimit" should be given in the form, "e_min e_max".')
        print('DOS Arguments:\n', dos_args)

        print('\n[User Arguments for Thermal Property]')
        if thermal:
            # Check and fill the thermal_args dictionary
            thermal_args['flag'] = thermal
            thermal_args['temp'] = range(tmin, tmax, tstep)
            thermal_args['legend_loc'] = dos_args['legend_loc']
        print('Thermal Arguments:\n', thermal_args)

        print('\n[User Arguments for Phonon Band]')
        if band:
            # Check and fill the band_args dictionary
            band_args['flag'] = band
            band_args['color'] = color_band
            band_args['option'] = option_band
            band_args['bar_loc'] = bar_loc_band
            band_args['bar_label'] = bar_label_band
            if k_label_band is None:
                band_args['k_label'] = []
            else:
                band_args['k_label'] = [label for label in k_label_band.strip().split()]

            if atom_band is None:
                if option_dos != 'plain':
                    print('Caution:')
                    print('"atom_band" should be given, in the Index of selected atoms, to plot projected Band.')
                    print('"option_band" is changed from "{0}" to "plain".'.format(option_dos))
                    band_args['option'] = 'plain'
                band_args['atom'] = None
            else:
                _atom_band = atom_band.strip().split()
                atom_band = [int(value) for value in _atom_band]
                band_args['atom'] = atom_band

            if elimit is None:
                band_args['elimit'] = None
            elif dos_args.get('elimit') is not None:
                band_args['elimit'] = dos_args['elimit']
            else:
                _elimit = elimit.strip().split()
                elimit = [int(value) for value in _elimit]
                if len(elimit) == 2:
                    band_args['elimit'] = elimit
                else:
                    raise Exception('"elimit" should be given in the form, "e_min e_max".')
        print('Band Arguments:\n', band_args)

        print('\n[User Arguments for Phonon Mode]')
        if mode:
            # Check and fill the mode_args dictionary
            mode_args['flag'] = mode
            mode_args['index'] = [int(ind) for ind in ind_mode.strip().split()]
            k_point = [float(kpt) for kpt in kpt_mode.strip().split()]
            if len(k_point) == 3:
                mode_args['k_point'] = k_point
            else:
                raise Exception('The length of "k_point" should be 3.')
        print('Mode Arguments:\n', mode_args)

        information_post_process = {'user_arguments': user_args,
                                    'files': files,
                                    'dos_arguments': dos_args,
                                    'thermal_arguments': thermal_args,
                                    'band_arguments': band_args,
                                    'mode_arguments': mode_args}

    # start process
    if _process == 'pre-process':
        print('\n#########################################')
        print('\tThis is pre-process...')
        print('#########################################')

        # define process
        print('\n>>>>>> Defining pre-process for displacement...')

        pre = PreProcess()

        # define user arguments
        pre.set_user_arg(user_args)
        _pre_user_arg = [{'displacement': pre.user_arg.displacement},
                         {'enlargement': ', '.join([str(_) for _ in pre.user_arg.enlargement])},
                         {'periodicity': ', '.join([str(_) for _ in pre.user_arg.periodicity])}]

        # define unit_cell
        pre.set_unit_cell(in_file=files.get('unit_cell_file'),
                          code_name=user_args.get('dft_code'))
        _pre_unit_cell = [{'lattice_matrix': dict(iter_lattice_yaml(pre.unit_cell.lattice_matrix))},
                          {'num_atom': int(pre.unit_cell.num_atom.sum())},
                          {'selected_atom_index': ', '.join([str(_) for _ in pre.unit_cell.atom_true])},
                          {'atoms': list(iter_atom_yaml(pre.unit_cell))}]

        # make a super_cell
        print('Writing supercell...')
        pre.set_super_cell(out_file=working_dir + '/SUPERCELL',
                           comment='Supercell',
                           write_file=True,
                           code_name=user_args.get('dft_code'))

        # Build a set of super_cells with displacements
        print('Writing displaced supercells...')
        pre.write_displace_cell(out_file=files.get('unit_cell_file'),
                                code_name=user_args.get('dft_code'))

        # Record this pre-process
        serialized_yaml_pre_process = [{'unit_cell_file': os.path.basename(files.get('unit_cell_file'))},
                                       {'dft_code_name': user_args.get('dft_code')},
                                       {'user_arguments': _pre_user_arg},
                                       {'unit_cell': _pre_unit_cell}]

        with open('pre_process.yaml', 'w') as outfile:
            yaml.dump(serialized_yaml_pre_process, outfile)

        print('Done.')

    else:
        print('\n#########################################')
        print('\tThis is post-process...')
        print('#########################################')

        if dos_args.get('flag', False):
            # define process
            print('\n>>>>>> Defining process for DOS...')
            post = PostProcess(in_file_unit_cell=files.get('unit_cell_file'),
                               in_file_super_cell=files.get('super_cell_file'))

            # define user arguments
            post.set_user_arg(dict_args=user_args)
            print('Index of selected atoms:\n', post.unit_cell.atom_true)

            # define reciprocal lattice
            post.set_reciprocal_lattice()
            print('Reciprocal lattice:\n', post.reciprocal_matrix)

            # construct Born-von Karman force constants
            print('Setting force constants...')
            post.set_force_constant(force_files=files.get('force_file'),
                                    code_name=user_args.get('dft_code'))

            # set k-points
            print('Setting k-points...')
            if not files.get('k_point_file_dos', False):
                raise Exception('K-points file for DOS should be given')
            else:
                post.set_k_points(k_file=files.get('k_point_file_dos'))

            # construct Dynamical matrix(q)
            print('Constructing dynamical matrix(q) and Evaluating phonon...')
            post.eval_phonon()

            print('DOS analysis is in progress...')
            from InterPhon.analysis import DOS
            post.dos = DOS(process=post, sigma=dos_args.get('sigma'), num_dos=dos_args.get('num_dos'))
            post.dos.set()
            post.dos.write(out_folder=working_dir)
            post.dos.plot(atoms=dos_args.get('atom'),
                          elimit=dos_args.get('elimit'),
                          color=dos_args.get('color'),
                          option=dos_args.get('option'),
                          orientation=dos_args.get('orientation'),
                          legends=dos_args.get('legend'),
                          legend_location=dos_args.get('legend_loc'),)

            if thermal_args.get('flag', False):
                print('Thermal analysis is in progress...')
                from InterPhon.analysis import ThermalProperty
                post.thermal = ThermalProperty(process=post, temp=thermal_args.get('temp'))
                post.thermal.set()
                post.thermal.write(out_folder=working_dir)
                post.thermal.plot(legend_location=thermal_args.get('legend_loc'))

        if band_args.get('flag', False):
            # define process
            print('\n>>>>>> Defining process for Band...')
            if dos_args.get('flag'):
                from copy import deepcopy
                post_band = deepcopy(post)
            else:
                post_band = PostProcess(in_file_unit_cell=files.get('unit_cell_file'),
                                        in_file_super_cell=files.get('super_cell_file'))

            # define user arguments
            post_band.set_user_arg(dict_args=user_args)
            if not dos_args.get('flag', False):
                print('Index of selected atoms:\n', post_band.unit_cell.atom_true)

            # define reciprocal lattice
            post_band.set_reciprocal_lattice()
            if not dos_args.get('flag', False):
                print('Reciprocal lattice:\n', post_band.reciprocal_matrix)

            # construct Born-von Karman force constants
            print('Setting force constants...')
            post_band.set_force_constant(force_files=files.get('force_file'),
                                         code_name=user_args.get('dft_code'))

            # set k-points
            print('Setting k-points...')
            if not files.get('k_point_file_band', False):
                raise Exception('K-points file for Band should be given')
            else:
                post_band.set_k_points(k_file=files.get('k_point_file_band'))

            # construct Dynamical matrix(q)
            print('Constructing dynamical matrix(q) and Evaluating phonon...')
            post_band.eval_phonon()

            print('Band analysis is in progress...')
            from InterPhon.analysis import Band
            post_band.band = Band(process=post_band)
            post_band.band.set()
            post_band.band.write(out_folder=working_dir)
            post_band.band.plot(k_labels=band_args.get('k_label'),
                                atoms=band_args.get('atom'),
                                elimit=band_args.get('elimit'),
                                color=band_args.get('color'),
                                option=band_args.get('option'),
                                colorbar_label=band_args.get('bar_label'),
                                colorbar_location=band_args.get('bar_loc'))

            if dos_args.get('flag', False):
                # from InterPhon.analysis import Band
                # post_band.band = Band(post_band)
                # post_band.band.set()
                # post_band.band.write(out_folder=working_dir)
                post_band.band.plot_with_dos(k_labels=band_args.get('k_label'),
                                             band_atoms=band_args.get('atom'),
                                             elimit=band_args.get('elimit'),
                                             band_color=band_args.get('color'),
                                             band_option=band_args.get('option'),
                                             colorbar_label=band_args.get('bar_label'),
                                             dos_object=post.dos,
                                             dos_atoms=dos_args.get('atom'),
                                             dos_color=dos_args.get('color'),
                                             dos_option=dos_args.get('option'),
                                             dos_legends=dos_args.get('legend'),
                                             legend_location=dos_args.get('legend_loc'),)

            if mode_args.get('flag', False):
                print('Mode analysis is in progress...')
                from InterPhon.analysis import Mode
                post_band.mode = Mode(process=post_band)
                post_band.mode.set(mode_inds=(0,), k_point=(0.0, 0.0, 0.0))
                post_band.mode.write(out_folder=working_dir)
                post_band.mode.plot(out_folder=working_dir,
                                    unit_cell=files.get('unit_cell_file'),
                                    code_name=user_args.get('dft_code'))

        print('Done.')


if __name__ == '__main__':
    main()
