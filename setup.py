from setuptools import setup, find_packages

setup(
    name = 'cosmo_peewee',
    version = '0.0.1',
    description = 'COSMO using the peewee ORM language.',
    author = 'Mees Fix',
    author_email = 'mfix@stsci.edu',
    keywords = ['astronomy'],
    classifiers = ['Programming Language :: Python',
                   'Programming Language :: Python :: 3',
                   'Development Status :: 1 - Planning',
                   'Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering :: Astronomy',
                   'Topic :: Scientific/Engineering :: Physics',
                   'Topic :: Software Development :: Libraries :: Python Modules'],
    packages = find_packages(),
    requires = ['numpy', 'scipy', 'astropy', 'matplotlib'],
    entry_points = {'console_scripts': ['cm_ingest=cosmo_peewee.database:ingest_all',
                                        'cm_monitors=cosmo_peewee.database:run_monitors'
                                        ],
    },
    install_requires = ['setuptools',
                        'numpy>=1.11.1',
                        'astropy>=1.0.1',
                        'peewee>=2.8.5',
                        'pymysql',
                        'matplotlib',
                        'scipy',
                        'fitsio',
                        'psutil'
                        ]
    )
