from setuptools import setup, find_packages

setup(
    name='lateral-tools',
    version='0.3',
    description='Utilities for seismic and lateral load engineering',
    author='andrew.mole@gmail.com',
    license=None,
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    install_requires=[
        'matplotlib>=3.0',
        'numpy',
        'scipy',
        'pandas',
        'openpyxl',
    ],
    python_requires='>=3.9',
)