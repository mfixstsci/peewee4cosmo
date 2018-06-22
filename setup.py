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
                                        'cm_monitors=cosmo_peewee.database:run_monitors',
                                        'drop_tables=cosmo_peewee.database:drop_tables',
                                        'cosmo_webapp=cosmo_peewee.webapp.main:run',
                                        'cosmo_webapp_debug=cosmo_peewee.webapp.main:run_debug',
					                    'cosmo_dark_period=cosmo_peewee.dark.dark_periodicity:main',
                                        'cosmo_gsagtab_by_date=cosmo_peewee.cci.gainmap_sagged_pixel_overplotter:gsagtab_plot_by_date',
                                        'cosmo_gsagtab_residual_plot=cosmo_peewee.cci.gainmap_sagged_pixel_overplotter:compare_and_plot_gsagtable_data_entry',
                                        'cosmo_gsagtab_creator=cosmo_peewee.cci.gainsag:make_gsagtab_db_entry',
                                        'cosmo_data_by_rootname=cosmo_peewee.database.utils:gather_all_data_rootname',
                                        'cosmo_gainmap_gsag_plot=cosmo_peewee.cci.gainmap_sagged_pixel_overplotter:plot_gainmap_and_gsagtab_by_hv_entry',
                                        ],
    },
    install_requires = ['setuptools',
                        'numpy>=1.11.1',
                        'astropy>=1.0.1',
                        'peewee==2.10.2',
                        'pymysql',
                        'matplotlib',
                        'scipy',
                        'psutil'
                        'networkx'
                        'fitsio'
                        'msgpack'
                        ]
    )
