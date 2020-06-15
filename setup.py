import setuptools


# Utility function to read the README file.
with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name="autophot",
    version="0.1.1",
    author="Sean Brennan",
    author_email="sean.brennan2@ucdconnect.ie",
    description="Automated Photometry of Transients",
    long_description=long_description,
    url='https://github.com/Astro-Sean/autophot',
    packages=setuptools.find_packages(),
    zip_safe=False,
    python_requires='>=3.5, <4',
    package_data={'': ['databases/*.yml','example/example.fits']
                  },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Selected collabortors',
        'Topic :: Software Development :: Build Tools',
        'License :: CC BY-SA 4.0',
        'Programming Language :: Python :: 3.6',
        ],

    project_urls={
        'Bug Reports': 'https://github.com/Astro-Sean/autophot/issues',
        'Source': 'https://github.com/Astro-Sean/autophot',
        'Homepage':'https://sn.ie'
        }
)
