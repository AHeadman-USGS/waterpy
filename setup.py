#!/usr/bin/env python
# -*- coding: utf-8 -*-

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("CHANGELOG.rst") as changelog_file:
    changelog = changelog_file.read()

with open("LICENSE") as license_file:
    license = license_file.read()


requirements = [
    "click",
    "jinja2",
    "matplotlib",
    "mpld3",
    "numpy",
    "pandas",
    "scipy",
]

test_requirements = [
    "pytest",
]


setup(
    name="waterpy",
    version="0.2",
    description=("A rainfall-runoff model that predicts the amount "
                 "of water flow in rivers."),
    long_description=readme + "\n\n" + changelog,
    author="Alexander Headman",
    author_email="aheadman@usgs.gov",
    url="",
    packages=[
        "waterpy",
    ],
    package_dir={"waterpy": "waterpy"},
    package_data={"waterpy": ["templates/*.html"]},
    entry_points={"console_scripts": ["waterpy = waterpy.cli:main"]},
    include_package_data=True,
    install_requires=requirements,
    license=license,
    zip_safe=False,
    keywords="waterpy",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: CC0 1.0 Universal (CC0 1.0) Public Domain Dedication",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    test_suite="tests",
    tests_require=test_requirements
)
