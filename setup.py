import sys
from setuptools import setup


setup_requires = []

install_requires = ['Click', 'numpy', 'matplotlib', 'PyYAML']
# Optional: (include Atomic Simulation Environment (ASE) package to use the visualization of vibration mode)
# install_requires.append('ase')

packages_interphon = ['InterPhon',
                      'InterPhon.core',
                      'InterPhon.error',
                      'InterPhon.inout',
                      'InterPhon.util',
                      'InterPhon.analysis', ]

scripts_interphon = ['scripts/interphon.py', ]


if __name__ == '__main__':

    assert sys.version_info >= (3, 0), 'python>=3 is required'
    
    with open('InterPhon/__init__.py', 'r') as init_file:
        for line in init_file:
            if "__version__" in line:
                version = line.split()[2].strip('\"')
                break

    setup(name='InterPhon',
          version=version,
          description='A Python Package for Ab initio Interface Phonon Calculations within a 3D Electronic Structure Framework',
          url='https://github.com/inwonyeu/interphon',
          author='In Won Yeu',
          author_email='yeuiw@kist.re.kr',
          license='LGPLv2.1',
          packages=packages_interphon,
          install_requires=install_requires,  # The package written here will be installed with the current package.
          python_requires='>=3',
          setup_requires=setup_requires,
          # scripts=scripts_interphon,
          entry_points={'console_scripts': ['interphon = InterPhon.interphon:main', ], },
          )
