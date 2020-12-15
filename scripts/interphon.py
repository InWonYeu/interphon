import os
import glob
import click
import yaml
import numpy as np
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


def check_file_order(process, unit_cell_file, force_file, dft_code):
    assert len(process.unit_cell.atom_true * 6) == len(force_file), \
        "The number of force files is not consistent with pre-process"

    from InterPhon.core import UnitCell
    _dis_super_cell = UnitCell()
    _dis_super_position = process.super_cell.atom_cart.copy()
    _current_position = process.super_cell.atom_cart.copy()

    _enlarge = 1
    for ind, value in enumerate(process.user_arg.periodicity):
        if value:
            _enlarge = _enlarge * process.user_arg.enlargement[ind]

    for i, ind_T in enumerate(process.unit_cell.atom_true):
        for j, displace in enumerate(np.eye(3, dtype=float)):
            # Forward Displacement
            _dis_super_position[_enlarge * ind_T, 0:3] = \
                _current_position[_enlarge * ind_T, 0:3] + process.user_arg.displacement * displace

            _dis_super_cell.read_unit_cell(os.path.dirname(force_file[6 * i + 2 * j]) + '/' + unit_cell_file,
                                           code_name=dft_code)

            assert np.allclose(_dis_super_position, _dis_super_cell.atom_cart), \
                "The force files are not consistent with the {0} of pre-process".format(
                    unit_cell_file + '-' + '{{{0:0>4}..{1:0>4}}}'.format(1, len(process.unit_cell.atom_true) * 6))

            # Backward Displacement
            _dis_super_position[_enlarge * ind_T, 0:3] = \
                _current_position[_enlarge * ind_T, 0:3] - process.user_arg.displacement * displace

            _dis_super_cell.read_unit_cell(os.path.dirname(force_file[6 * i + 2 * j + 1]) + '/' + unit_cell_file,
                                           code_name=dft_code)

            assert np.allclose(_dis_super_position, _dis_super_cell.atom_cart), \
                "The force files are not consistent with the {0} of pre-process".format(
                    unit_cell_file + '-' + '{{{0:0>4}..{1:0>4}}}'.format(1, len(process.unit_cell.atom_true) * 6))


@click.command()
@click.argument('force_files',
                nargs=-1,
                # type=click.File('rb'),
                required=False)
@click.option('--option_file', '-option',
              type=click.Path(exists=True),
              help='File of option collections.')
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
              default='0.02',
              type=click.STRING,
              help='Displacement length (unit: Angst).',
              required=True,
              show_default=True)
@click.option('--enlargement', '-enlarge',
              default='2 2 1',
              type=click.STRING,
              help='Extension ratio along each a, b, c lattice direction.',
              required=True,
              show_default=True)
@click.option('--periodicity', '-pbc',
              default='1 1 0',
              type=click.STRING,
              help='Periodic (True or 1) or not (False or 0) along each a, b, c lattice direction.',
              required=True,
              show_default=True)
@click.option('--unitcell', '-c',
              default='POSCAR',
              # prompt='unit cell file',
              type=click.Path(exists=True),
              help='Unit cell file.')
@click.option('--supercell', '-sc',
              type=click.Path(exists=True),
              help='Supercell file.')
# Options for DOS write and plot
@click.option('--density_of_state', '-dos', 'dos', is_flag=True,
              default=False,
              required=True,
              show_default=True,
              help='Flag to DOS.')
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
              help='The Index of atoms to be projected in DOS plot.')
@click.option('--projection_legend_dos', '-legend_dos', 'legend_dos',
              type=click.STRING,
              help='Legends for the projected atoms.')
@click.option('--energy_limit', '-elimit', 'elimit',
              type=click.STRING,
              help='Energy (THz) limitation of DOS and Band plot.')
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
              show_default=True,
              help='Flag to thermal property.')
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
              show_default=True,
              help='Flag to phonon band')
@click.option('--kpoint_band', '-kband',
              type=click.Path(exists=True),
              help='K-point file for Band.')
@click.option('--kpoint_label_band', '-k_label_band', 'k_label_band',
              type=click.STRING,
              help='The label of high-symmetry k-points.')
@click.option('--projection_atom_band', '-atom_band', 'atom_band',
              type=click.STRING,
              help='The Index of atoms to be projected in Band plot.')
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
              show_default=True,
              help='Flag to phonon mode')
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
def main(force_files, option_file, process, dft, displacement, enlargement, periodicity,
         unitcell, supercell, kpoint_dos,
         dos, sigma, num_dos, atom_dos, legend_dos, elimit, color_dos, option_dos, orientation_dos, legend_loc_dos,
         thermal, tmin, tmax, tstep,
         band, kpoint_band, k_label_band, atom_band, color_band, option_band, bar_label_band, bar_loc_band,
         mode, ind_mode, kpt_mode):
    if option_file is not None:
        with open(option_file, 'r') as infile:
            lines = infile.readlines()

        arg_dict = {}
        for line in lines:
            if line.strip().split('#')[0]:
                _line = line.strip().split('#')[0]
                key = _line.strip().split('=')[0]
                value = _line.strip().split('=')[1]
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

    if force_files:
        process = False

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

        with open('pre_process.yaml') as yaml_file:
            pre_record = yaml.load(yaml_file, Loader=yaml.FullLoader)

        # Unit cell file
        if unitcell == 'POSCAR':  # default
            unitcell = pre_record[0].get('unit_cell_file')

        # Supercell file
        if supercell is None:
            supercell = pre_record[1].get('supercell_file')

        if dft == 'vasp':  # default
            dft = pre_record[3].get('dft_code')

        if displacement == '0.02':  # default
            displacement = pre_record[4].get('user_arguments')[0].get('displacement')

        if enlargement == '2 2 1':  # default
            enlargement = pre_record[4].get('user_arguments')[1].get('enlargement')

        if periodicity == '1 1 0':  # default
            periodicity = pre_record[4].get('user_arguments')[2].get('periodicity')

        user_args = {'dft_code': dft,
                     'displacement': displacement,
                     'enlargement': enlargement,
                     'periodicity': periodicity}

        dos_args = {'flag': dos}
        thermal_args = {'flag': thermal}
        band_args = {'flag': band}
        mode_args = {'flag': mode}

        # Check and fill the files dictionary
        files = {'unit_cell_file': os.path.abspath(unitcell),
                 'super_cell_file': os.path.abspath(supercell)}

        # K-point file
        if kpoint_dos is None:
            if kpoint_band is None:
                raise Exception('K-points file should be given')
            else:
                band_args['flag'] = True
                files['k_point_file_band'] = os.path.abspath(kpoint_band)  # only Band
        else:
            if kpoint_band is None:
                dos_args['flag'] = True
                files['k_point_file_dos'] = os.path.abspath(kpoint_dos)  # only DOS
            else:
                dos_args['flag'] = True
                files['k_point_file_dos'] = os.path.abspath(kpoint_dos)
                band_args['flag'] = True
                files['k_point_file_band'] = os.path.abspath(kpoint_band)  # Both DOS and Band

        # Force file
        if force_files is None:
            raise Exception('Force files should be given')
        else:
            if len(force_files) > 1:
                files['force_file'] = force_files
            else:
                force_path = glob.glob(force_files[0])
                force_path.sort()
                if not force_path:
                    raise Exception('Path "{0}" does not exist.'.format(force_files))
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
        if dos_args['flag']:
            # Check and fill the dos_args dictionary
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
        if thermal_args['flag']:
            # Check and fill the thermal_args dictionary
            thermal_args['temp'] = range(tmin, tmax, tstep)
            thermal_args['legend_loc'] = dos_args['legend_loc']
        print('Thermal Arguments:\n', thermal_args)

        print('\n[User Arguments for Phonon Band]')
        if band_args['flag']:
            # Check and fill the band_args dictionary
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
        if mode_args['flag']:
            # Check and fill the mode_args dictionary
            mode_args['index'] = [int(ind) for ind in ind_mode.strip().split()]
            k_point = [float(kpt) for kpt in kpt_mode.strip().split()]
            if len(k_point) == 3:
                mode_args['k_point'] = k_point
            else:
                raise Exception('The length of "k_point" should be 3.')
        print('Mode Arguments:\n', mode_args)

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
                         {'enlargement': ' '.join([str(_) for _ in pre.user_arg.enlargement])},
                         {'periodicity': ' '.join([str(_) for _ in pre.user_arg.periodicity])}]

        # define unit_cell
        pre.set_unit_cell(in_file=files.get('unit_cell_file'),
                          code_name=user_args.get('dft_code'))
        _pre_unit_cell = [{'lattice_matrix': dict(iter_lattice_yaml(pre.unit_cell.lattice_matrix))},
                          {'num_atom': int(pre.unit_cell.num_atom.sum())},
                          {'selected_atom_index': ', '.join([str(_) for _ in pre.unit_cell.atom_true])},
                          {'atoms': list(iter_atom_yaml(pre.unit_cell))}]

        # make a super_cell
        supercell_file = 'SUPERCELL'
        print('Writing supercell... ---> {0}'.format(supercell_file))
        pre.set_super_cell(out_file=working_dir + '/' + supercell_file,
                           comment='Supercell',
                           write_file=True,
                           code_name=user_args.get('dft_code'))

        # Build a set of super_cells with displacements
        _displaced_supercell = os.path.basename(files.get('unit_cell_file')) + '-' \
                               + '{{{0:0>4}..{1:0>4}}}'.format(1, len(pre.unit_cell.atom_true) * 6)
        print('Writing displaced supercells... ---> {0}'.format(_displaced_supercell))
        pre.write_displace_cell(out_file=files.get('unit_cell_file'),
                                code_name=user_args.get('dft_code'))

        # Record this pre-process
        serialized_yaml_pre_process = [{'unit_cell_file': os.path.basename(files.get('unit_cell_file'))},
                                       {'supercell_file': supercell_file},
                                       {'displaced_supercell_files': _displaced_supercell},
                                       {'dft_code': user_args.get('dft_code')},
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
                               in_file_super_cell=files.get('super_cell_file'),
                               code_name=user_args.get('dft_code'))

            # define user arguments
            post.set_user_arg(dict_args=user_args)
            _post_user_arg = [{'displacement': post.user_arg.displacement},
                              {'enlargement': ' '.join([str(_) for _ in post.user_arg.enlargement])},
                              {'periodicity': ' '.join([str(_) for _ in post.user_arg.periodicity])}]
            print('Index of selected atoms:\n', post.unit_cell.atom_true)

            # define reciprocal lattice
            post.set_reciprocal_lattice()
            print('Reciprocal lattice:\n', post.reciprocal_matrix)

            # construct Born-von Karman force constants
            print('Setting force constants...')
            check_file_order(post,
                             os.path.basename(files.get('unit_cell_file')),
                             files.get('force_file'),
                             user_args.get('dft_code'))
            post.set_force_constant(force_files=files.get('force_file'),
                                    code_name=user_args.get('dft_code'))

            # set k-points
            print('Setting k-points from {0}...'.format(os.path.basename(files.get('k_point_file_dos'))))
            if not files.get('k_point_file_dos', False):
                raise Exception('K-points file for DOS should be given')
            else:
                post.set_k_points(k_file=files.get('k_point_file_dos'))

            # construct Dynamical matrix(q)
            print('Constructing dynamical matrix(q) and Evaluating phonon...')
            post.eval_phonon()

            print('DOS analysis is in progress... ---> total_dos.dat and projected_dos.dat')
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
                print('Thermal analysis is in progress... ---> thermal_properties.dat')
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
                                        in_file_super_cell=files.get('super_cell_file'),
                                        code_name=user_args.get('dft_code'))

            # define user arguments
            post_band.set_user_arg(dict_args=user_args)
            _post_user_arg = [{'displacement': post_band.user_arg.displacement},
                              {'enlargement': ' '.join([str(_) for _ in post_band.user_arg.enlargement])},
                              {'periodicity': ' '.join([str(_) for _ in post_band.user_arg.periodicity])}]
            if not dos_args.get('flag', False):
                print('Index of selected atoms:\n', post_band.unit_cell.atom_true)

            # define reciprocal lattice
            post_band.set_reciprocal_lattice()
            if not dos_args.get('flag', False):
                print('Reciprocal lattice:\n', post_band.reciprocal_matrix)

            # construct Born-von Karman force constants
            print('Setting force constants...')
            check_file_order(post_band,
                             os.path.basename(files.get('unit_cell_file')),
                             files.get('force_file'),
                             user_args.get('dft_code'))
            post_band.set_force_constant(force_files=files.get('force_file'),
                                         code_name=user_args.get('dft_code'))

            # set k-points
            print('Setting k-points from {0}...'.format(os.path.basename(files.get('k_point_file_band'))))
            if not files.get('k_point_file_band', False):
                raise Exception('K-points file for Band should be given')
            else:
                post_band.set_k_points(k_file=files.get('k_point_file_band'))

            # construct Dynamical matrix(q)
            print('Constructing dynamical matrix(q) and Evaluating phonon...')
            post_band.eval_phonon()

            print('Band analysis is in progress... ---> band.dat')
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
                print('Mode analysis is in progress... ---> XDATCAR_phonon_[mode_index]_[k_point]')
                from InterPhon.analysis import Mode
                post_band.mode = Mode(process=post_band)
                post_band.mode.set(mode_inds=mode_args.get('index'), k_point=mode_args.get('k_point'))
                post_band.mode.write(out_folder=working_dir)
                post_band.mode.plot(out_folder=working_dir,
                                    unit_cell=files.get('unit_cell_file'),
                                    code_name=user_args.get('dft_code'))

        # Record this post-process
        _post_files = []
        for key, value in files.items():
            if isinstance(value, str):
                value = os.path.basename(value)
            _post_files.append({key: value})

        _dos_args = []
        for key, value in dos_args.items():
            if (not isinstance(value, bool)) and (not isinstance(value, int)) and (not isinstance(value, float)):
                if not isinstance(value, list):
                    value = str(value)
            _dos_args.append({key: value})

        _thermal_args = []
        for key, value in thermal_args.items():
            if (not isinstance(value, bool)) and (not isinstance(value, int)) and (not isinstance(value, float)):
                if not isinstance(value, list):
                    value = str(value)
            _thermal_args.append({key: value})

        _band_args = []
        for key, value in band_args.items():
            if (not isinstance(value, bool)) and (not isinstance(value, int)) and (not isinstance(value, float)):
                if not isinstance(value, list):
                    value = str(value)
            _band_args.append({key: value})

        _mode_args = []
        for key, value in mode_args.items():
            if (not isinstance(value, bool)) and (not isinstance(value, int)) and (not isinstance(value, float)):
                if not isinstance(value, list):
                    value = str(value)
            _mode_args.append({key: value})

        serialized_yaml_post_process = [{'files': _post_files},
                                        {'dft_code': user_args.get('dft_code')},
                                        {'user_arguments': _post_user_arg},
                                        {'dos_arguments': _dos_args},
                                        {'thermal_arguments': _thermal_args},
                                        {'band_arguments': _band_args},
                                        {'mode_arguments': mode_args}]

        with open('post_process.yaml', 'w') as outfile:
            yaml.dump(serialized_yaml_post_process, outfile)

        print('Done.')


if __name__ == '__main__':
    main()
