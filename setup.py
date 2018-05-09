import os
from glob import glob
from setuptools import setup, find_packages


package = 'jetstream'
version_file = os.path.join(os.path.dirname(__file__), 'VERSION')
with open(version_file, 'r') as fp:
    __version__ = fp.readline().strip()


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name=package,
    version=__version__,
    author="Ryan Richholt",
    author_email="ryan@tgen.org",
    url="https://github.com/tgen/jetstream",
    description="NGS analysis pipeline at TGen.",
    long_description=read('README.md'),
    long_description_content_type="text/markdown",
    keywords="ngs pipeline automation",
    packages=find_packages(),
    include_package_data=True,
    python_requires='>=3',
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
    ],
    install_requires=[
        'networkx',
        'ruamel.yaml',
        'ulid-py',
        'tempstore',
        'requests',
        'jinja2'
    ],
    extras_require={
        'dev': [
            'sphinx-argparse',
            'sphinx-autobuild',
            'sphinx-rtd-theme',
            'bumpversion'
        ]
    },
    scripts=[s for s in glob('scripts/**', recursive=True) if os.path.isfile(s)],
    package_data={
        'jetstream': ['templates/*', 'etc/*']
    },
    entry_points={
        'console_scripts': [
            'jetstream=jetstream.cli.jetstream:main',
        ],
    }
)
