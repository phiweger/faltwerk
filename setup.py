from setuptools import setup#, find_packages
# https://packaging.python.org/tutorials/packaging-projects/
# https://docs.pytest.org/en/latest/goodpractices.html


with open('README.md', 'r') as fh:
    long_description = fh.read()


setup(
    name='faltwerk',
    version='0.3',
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
        'biopython==1.79',
        'screed==1.0.5',
        'py3Dmol==1.8.1',
        'markov-clustering==0.0.6.dev0',
        'hdbscan==0.8.28',
        'networkx==2.6.3',
        'numpy==1.21.6',
        'pandas==1.3.5',
        'matplotlib==3.5.3',
        'pytest==7.1.2',
        'libpysal==4.6.2',
        'esda==2.4.3',
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