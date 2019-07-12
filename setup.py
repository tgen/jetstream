import os
from setuptools import setup, find_packages

package = 'jetstream'
version_file = os.path.join(os.path.dirname(__file__), 'VERSION')
with open(version_file, 'r') as fp:
    __version__ = fp.readline().strip()


def read(fname):
    with open(os.path.join(os.path.dirname(__file__), fname)) as fp:
        return fp.read()

setup(
    name=package,
    version=__version__,
    author='Ryan Richholt',
    author_email='ryan@tgen.org',
    url='https://github.com/tgen/jetstream',
    description='NGS analysis pipeline at TGen.',
    long_description=read('README.md'),
    long_description_content_type='text/markdown',
    keywords='ngs pipeline automation',
    packages=find_packages(exclude=('test',)),
    include_package_data=True,
    python_requires='>=3',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Utilities',
    ],
    install_requires=[
        'networkx',
        'pyyaml',
        'ulid-py',
        'jinja2',
        'filelock',
        'confuse'
    ],
    package_data={
        'jetstream': ['config_default.yaml']
    },
    entry_points={
        'console_scripts': [
            'jetstream=jetstream.cli:main',
        ],
    }
)
