from setuptools import setup, find_packages

setup(
    name='biocal',
    version='2.0.1',
    description='Using sedimentological priors for more accurate calibration of 14C determinations from bioturbated sediment archives.',
    author='Bryan C. Lougheed',
    author_email='bryan.lougheed@outlook.com',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy',
        'numba'
    ],
)
