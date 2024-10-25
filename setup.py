from setuptools import setup

setup(
    name='obsprep',
    version='0.1',
    description='MITgcm preprocessing observational data module',
    author='Matthew Goldberg, Shoshana Reich, Caeli Griffin',
    author_email='matthew.goldberg10@gmail.com',
    packages=['obsprep'],
    include_package_data=True,
    test_suite='tests',
    install_requires=[
        'ecco_v4_py',
        'xarray',
        'scipy',
        'numpy',
    ],
    tests_require=[
        'pytest',
    ],
)
