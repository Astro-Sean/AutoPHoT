import setuptools

setuptools.setup(
    name="autophot",
    version="0.0.1",
    author="Sean Brennan",
    author_email="sean.brennan2@ucdconnect.ie",
    description="Automated Photometry of Transients",
    url='https://github.com/Astro-Sean/autophot',
    packages=setuptools.find_packages(),
    python_requires='>=3.5, <4',
    package_data={'': ['databases/*.yml','example/example.fits']
                  },
    classifiers=[
        'Development Status :: 3 - Alpha',

        # WHo project is intended for
        'Intended Audience :: Selected collabortors',
        'Topic :: Software Development :: Build Tools',

        # License
        'License :: CC BY-SA 4.0',

        # Suitable Languages
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        ],

    project_urls={
        'Bug Reports': 'https://github.com/Astro-Sean/autophot/issues',
        'Source': 'https://github.com/Astro-Sean/autophot',
        'Homepage':'https://sn.ie'
        }
)
