try:
    from setuptools import setup, find_packages
except (ImportError, ModuleNotFoundError) as e:
    print("Please install setuptools and try again.\npip install setuptools")
    import sys
    sys.exit(1)
else:
    setup(
        name='eharmonize',
        author="Alyssa H. Zhu",
        description="eharmonize (ENIGMA harmonization) is a python-based tool for harmonizing outputs of the ENIGMA-DTI pipeline to a provided reference",
        version='0.0.0',
        packages=find_packages(exclude=["data"]),
        include_package_data=True,
        package_data={"": ["data/*"]},
        install_requires=[
            'click',
            'pandas',
            'numpy',
            'seaborn',
            'neuroHarmonize @ git+ssh://git@github.com/ahzhu/neuroHarmonize.git'
        ],
        scripts=[
            'bin/model2csv'
            ],
        entry_points={
            'console_scripts': [
                'eharmonize = eharmonize.scripts.cli:eharmonize'
                ]
            }
    )
