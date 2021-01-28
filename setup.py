from setuptools import setup, Extension

setup(
    name = 'spec1d-package',
    version = '3.8',
    description = '1d dispersion spectrograph extraction software.',
    url = None,
    author = 'Samuel Gill',
    author_email = 'samuel.gill@warwick.ac.uk',
    license = 'GNU',
    packages=['spec1d'],
    package_dir={'spec1d': 'src/spec1d'},
    #package_data={'spec1d': ['data/gr5_-4_ref_spectra.dat', 'data/gr5_-4_ref_lines.dat']},
    package_data={'spec1d': ['data/*.dat']},
    scripts=['Utils/spupnic', 'Utils/spupnic_halpha' , 'Utils/check_calibration'],
    install_requires=['numba']
)