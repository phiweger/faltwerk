from setuptools import setup, find_packages
# https://packaging.python.org/tutorials/packaging-projects/
# https://docs.pytest.org/en/latest/goodpractices.html


with open('README.md', 'r') as fh:
    long_description = fh.read()


setup(
    name='faltwerk',
    version='0.2.32',
    author='Adrian Viehweger',
    author_email='adrian.viehweger@medizin.uni-leipzig.de',
    description='Spatial analysis of protein structures',
    long_description=long_description,
    long_description_content_type='text/markdown',
    # url='',
    license='BSD 3-clause',
    
    # next 2 lines for pytest
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],

    install_requires=[
        # 'screed',
        # 'tqdm',
        # 'py3dmol',
        # 'markov-clustering',
        # 'hdbscan',
        # 'pdb-tools',
        # 'altair',
        # 'vega_datasets',
        # 'networkx',
        # 'numpy',
        # 'pandas',
        # 'matplotlib',
        # 'pytest',
        # 'libpysal',
        # 'spreg',
        # 'esda',
        # 'geopandas',
    ],

    # https://click.palletsprojects.com/en/7.x/setuptools/#testing-the-script
    # https://click.palletsprojects.com/en/7.x/setuptools/#scripts-in-packages
    # packages=find_packages(),
    # entry_points='''
    #     [console_scripts]
    #     nanotext=nanotext.__main__:cli
    # ''',
    packages=['faltwerk'],
    include_package_data=True,
    # package_data={'faltwerk': ['faltwerk/data/ligands/*.tsv']},

    # data_files=[
    #     ('interacdome', [
    #         'data/ligands/InteracDome_v0.3-confident.tsv',
    #         'data/ligands/InteracDome_v0.3-representableNR.tsv',
    #         'data/ligands/InteracDome_v0.3-representable.tsv',
    #         ]
    #     )
    # ],
)