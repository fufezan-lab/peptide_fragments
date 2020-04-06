#!/usr/bin/env python

"""The setup script."""

from setuptools import find_packages, setup

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('requirements.txt') as requirement_file:
    requirements = requirement_file.readlines()

with open('peptide_fragmentor/version.txt') as version_file:
    version = version_file.read()

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest>=3', ]


setup(
    author="Manuel Kösters, Christian Fufezan",
    author_email='christian@fufezan.net',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Library to calculate fragment ions of peptides",
    install_requires=requirements,
    license="MIT license",
    long_description=readme,
    include_package_data=True,
    name='peptide_fragments',
    packages=find_packages(include=['peptide_fragmentor', 'peptide_fragmentor.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/fufezan-lab/peptide_fragments',
    version=version,
    zip_safe=False,
)
